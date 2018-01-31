c PHONON module of the NMscatt program package
c
c Author: Franci Merzel
c         

      Program phonon_calc
      implicit none
      integer nmax,nelm,i,j,k,imd,ii,jj,i0,i1,i2,lwork,INFO,nat,nn,mxat
      integer j1,ik1,ik2,ik3
      parameter(mxat=3600) !(mxat=10000)
      parameter(nmax=3*mxat)
      integer lst(0:mxat)
      real*8 xx,ss,eigval,x,y,z,qq(3),dotql,pi,fct(3),ref
      real*8 brva,brvb,brvc,vlc,rcpa(3),rcpb(3),rcpc(3),pfc
      real*8 tr1,tr2,vv1,vv2,sqr_mass,amltp,dist,ad,cutoff,dpx,dpy,dpz
      real*8 RWORK,sgn,amx,avrg,conv,re,hh,dmin,distc,dx,dy,dz
      real*4 wr,wi
      double complex ic,dd,WORK
      logical lf_nz,lf_co,flg,lf_binr,deut
      character*4  ch
      character*4  arg
      character*7  arg1
      character*11 name,name_crd

      common/ blckr / eigval(nmax),RWORK(3*nmax-2)
      common/ blckc / dd(nmax,nmax),WORK(2*nmax-1)
      common/ blck1 / sqr_mass(nmax),x(mxat),y(mxat),z(mxat)
      common/ blck2 / dist(mxat,mxat),hh(nmax,nmax)
      common/ blck3 / wr(nmax),wi(nmax)
      common/ vct   / brva(3),brvb(3),brvc(3)

c ####################################################################
c #### DATA given at the command line ################################
c suffix for the eigen-value/vector output files
      call getarg(1,arg)
      read(arg,*)ch

c fractional coordinates for k-point
      call getarg(2,arg1)
      read(arg1,*)fct(1)
      call getarg(3,arg1)
      read(arg1,*)fct(2)
      call getarg(4,arg1)
      read(arg1,*)fct(3)

      write(*,'(a,3f8.5)')'k-vect (rtnl-BZ):',(fct(k),k=1,3)

c ####  SIMULATION DATA required to specify  ##########################
c ####  for each system according to the simulation protocol  #########
c 
      cutoff   = 7.5          ! cutoff radius

c Bravais lattice vectors
      data brva/ 15.00,  0.00,  0.00/
      data brvb/  0.00, 15.00,  0.00/
      data brvc/  0.00,  0.00, 20.25/

      lf_binr  = .FALSE.      ! logical flag (formatted/binary) - (0/1)
      name     = 'hessian'

      name_crd = 'coord'      ! CHARMM coordinate file 

      deut     =  .TRUE.      ! deuterated case

c #####################################################################

      pi = 4.d0*atan(1.d0)
      ic = cmplx(0.0,1.0)

c conversion factor for freq: sqrt(kcal/mol/amu/A^2) --> cm-1
      conv = 33.356*sqrt(418.415)/(2.0*pi)


c --------------------------------------------------------------------

      call mix(brva,brvb,brvc,vlc)
      pfc = 2.*pi/vlc

c reciprocal vectors
      call cross(brvb,brvc,pfc,rcpa)
      call cross(brvc,brva,pfc,rcpb)
      call cross(brva,brvb,pfc,rcpc)
      
      write(*,'(a,3f12.5)')'rcp_v1',rcpa(1),rcpa(2),rcpa(3)
      write(*,'(a,3f12.5)')'rcp_v2',rcpb(1),rcpb(2),rcpb(3)
      write(*,'(a,3f12.5)')'rcp_v3',rcpc(1),rcpc(2),rcpc(3)
      write(*,*)
c  q - vector (3D)
      do k=1,3
        qq(k) = fct(1)*rcpa(k) + fct(2)*rcpb(k) + fct(3)*rcpc(k)
      end do

      open(3,file='eig_val_'//ch,form='formatted',status ='unknown')
      open(4,file='eig_vec_'//ch,form='unformatted',status ='unknown')
      write(3,'(a,3f8.4)')'# q_i = fct_i*pi/L_i; fct_i:',
     1     fct(1),fct(2),fct(3)

c  assigning atom types and atom masses
      call assig_mass(nat,name_crd,deut)

c filling up the distance matrix
      do i=1,nat
        do j=i+1,nat
          ad = sqrt((x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2)
          dist(i,j) = ad
          dist(j,i) = ad
        end do
      end do

c ===== reading the upper triangle of the Hessian matrix =====
c 
      if (lf_binr) then
c hessian is the binary - unformatted file
        open(2, file=name,
     1          form ='unformatted', status ='unknown')
        read(2) nelm

        if (nelm/3.ne.nat) then
          write(*,*)'<coord> and <hessian> mismatch!!',nelm/3,nat
          stop
        end if
        ii = 0
        jj = 0
        i0 = 0
 11     continue
          read(2,err=1000,end=2000) i1,i2,xx
          ii = ii + 1
          hh(i1,i2) = xx
          if (abs(xx).gt.1.0e-6) then
            jj = jj + 1
          end if
          if (i0.ne.i1) then
            if (mod(i1,100).eq.0) write(*,'(a,i7,a)')
     1                                    'Passed:',i1,' lines'
            i0 = i1
          end if
        goto 11
1000    write(*,*)'Error While reading!'
        stop
2000    continue
        close(2)
      else
c hessian is the ASCII - formatted file
        open(2, file=name,
     1          form ='formatted', status ='unknown')
        read(2,*) nelm

        if (nelm/3.ne.nat) then
          write(*,*)'<coord> and <hessian> mismatch!!',nelm/3,nat
          stop
        end if
        ii = 0
        jj = 0
        i0 = 0
 22     continue
          read(2,*,err=3000,end=4000) i1,i2,xx
          ii = ii + 1
          hh(i1,i2) = xx
          if (abs(xx).gt.1.0e-6) then
            jj = jj + 1
          end if
          if (i0.ne.i1) then
            if (mod(i1,100).eq.0) write(*,'(a,i7,a)')
     1                                    'Passed:',i1,' lines'
            i0 = i1
          end if
        goto 22
3000    write(*,*)'Error While reading!'
        stop
4000    continue
        close(2)
      end if

      write(*,*)
      write(*,'(a,f10.2,a)')'Matrix contains  :',
     1      float(jj)/float(ii)*100.0,'% non-zero elements'
      write(*,*)
c ================================

c constructing dynamical matrix 
      do i=1,nelm
        i1 = (i-1)/3 + 1
        do j=i,nelm
          j1 = (j-1)/3 + 1
          if (abs(hh(i,j)).gt.0.0001) then
            if (dist(i1,j1).lt.cutoff) then
              dotql = 0.d0
            else
              dmin = cutoff
              ii = 0
              do ik1 = -1,1
                do ik2 = -1,1
                  do ik3 = -1,1
                    if ((abs(ik1)+abs(ik2)+abs(ik3)).ne.0) then
                      dx = float(ik1)*brva(1) + float(ik2)*brvb(1) +
     1                     float(ik3)*brvc(1)
                      dy = float(ik1)*brva(2) + float(ik2)*brvb(2) +
     1                     float(ik3)*brvc(2) 
                      dz = float(ik1)*brva(3) + float(ik2)*brvb(3) +
     1                     float(ik3)*brvc(3)
                      distc = sqrt((x(i1)-(x(j1)+dx))**2 + 
     1                             (y(i1)-(y(j1)+dy))**2 + 
     1                             (z(i1)-(z(j1)+dz))**2)
                      if (distc.lt.cutoff) then
                        ii = ii + 1
                        if (distc.lt.dmin) then
                          dmin = distc
                          dpx = dx
                          dpy = dy
                          dpz = dz
                        end if
                      end if
                    end if
                  end do
                end do
              end do
              dotql = qq(1)*dpx + qq(2)*dpy + qq(3)*dpz
            end if
            dd(i,j) = cmplx(hh(i,j),0.d0)*exp(-ic*dotql)
          else
            dd(i,j) = cmplx(0.d0,0.d0)
          end if
        end do
      end do
c =========================================================

      tr1 = 0.0d0
c multiplying by the mass matrix
      do i=1,nelm
        do j=i,nelm
          dd(i,j) = dd(i,j)/sqr_mass(i)/sqr_mass(j)
        end do
        tr1 = tr1 + re(dd(i,i))
      end do
      write(*,*)'trace:',tr1
      write(*,*)

      write(*,*)'Dimension of the matrix to diagonalize:',nelm
      lwork = 2*nelm-1
      nn = nmax
c      CALL CHEEV('V','U',nelm,dd,nn,eigval,WORK,lwork,RWORK,INFO)
      CALL ZHEEV('V','U',nelm,dd,nn,eigval,WORK,lwork,RWORK,INFO)
      write(*,*)'Diagonalizaton finished!'
      write(*,*)
      ss = 0.0
      write(3,'(a)')
     1'#         kcal/mol/amu/A^2      cm-1'
      do i=1,nelm
        if (eigval(i).gt.0.0d0) then
          write(3,'(i7,e15.3,f17.7)')i,eigval(i),sqrt(eigval(i))*conv
        else if (abs(eigval(i)).lt.1.0d-6) then
          write(3,'(i7,e15.3,f17.7)')i,eigval(i),0.0d0
        else
          write(3,'(i7,e15.3)')i,eigval(i)
        end if
        ss = ss + eigval(i)
      end do
      tr2 = ss
      close(3)
      write(*,'(a,2e17.8,5x,f10.3)')'comparison of the traces',tr1,tr2,
     1      tr1-tr2
      
      do i=1,nelm
        do j=1,nelm
          wr(j) = re(dd(i,j))
          wi(j) = imag(dd(i,j))
        end do
        write(4)(wr(j),wi(j),j=1,nelm)
      end do
      close(4)
      end

      real*8 function re(cmpl)
      double complex cmpl
      re = 0.5*(cmpl+conjg(cmpl))
      return
      end

      subroutine assig_mass(nat,name_crd,deut)
      implicit none
      integer nmax,nat,i,ii,ir,j,mxat
      parameter(mxat=3600) !(mxat=10000)
      parameter(nmax=3*mxat)
      real*8 sqr_mass,x,y,z,amass
      character*4 aty
      character*6 ars
      character*11 name_crd

      logical wt,oh,hn,ht,hr,r1,r2,r3,r4,r5,r6,r7,r8,exch,deut

      common/ blck1 / sqr_mass(nmax),x(mxat),y(mxat),z(mxat)

      open(9, file=name_crd, form   ='formatted',
     1                       status ='unknown')
c  suppose 4 lines off-set
c      read(9,*)
      read(9,*)
      read(9,*)
      read(9,*)

      read(9,*)nat
      do i=1,nat
        read(9,'(2i5,a6,a4,3f10.5)')ii,ir,ars,aty,x(i),y(i),z(i)
c sorting out atom types        

        if (aty(1:1).eq.'H') then  ! hydrogen

c the following block checks for exchangable hydrogens 
c in case of deuteration
          exch = .false. 
          if (deut) then
            wt = (ars(2:5).eq.'TIP3')   ! selecting water hydrogens
            oh = (aty(2:4).eq."2' ")    ! OH hydrogens
            hn = (aty(2:3).eq.'N ')     ! NH hydrogens
            ht = (aty(2:2).eq.'T')      ! n-terminal hydrogens

            r1 = (ars(2:4).eq.'ARG').and.
     1 ((aty(2:2).eq.'E').or.(aty(2:2).eq.'H'))
            r2 = (ars(2:4).eq.'ASN').and.
     1 (aty(2:2).eq.'D')
            r3 = (ars(2:4).eq.'GLN').and.
     1 (aty(2:2).eq.'E')
            r4 = (ars(2:4).eq.'LYS').and.
     1 (aty(2:2).eq.'Z')
            r5 = (ars(2:4).eq.'SER').and.
     1 (aty(2:2).eq.'G')
            r6 = (ars(2:4).eq.'THR').and.
     1 (aty(2:2).eq.'G')
            r7 = (ars(2:4).eq.'TRP').and.
     1 (aty(2:3).eq.'E1')
            r8 = (ars(2:4).eq.'TYR').and.
     1 (aty(2:2).eq.'H')
            
            hr = (r1.or.r2.or.r3.or.r4.or.r5.or.r6.or.r7.or.r8)
            exch = (wt.or.oh.or.hn.or.ht.or.hr)
          end if
          if (exch) then
c deuterium 
            amass = 2.009
          else
c hydrogen
            amass =  1.008
          end if
        else if (aty(1:3).eq.'CLA') then
          amass = 35.45
        else if (aty(1:1).eq.'C') then
          amass = 12.011
        else if (aty(1:1).eq.'N') then
          amass = 14.007
        else if (aty(1:1).eq.'O') then
          amass = 15.9994
        else if (aty(1:1).eq.'P') then
          amass = 30.974
        else if (aty(1:1).eq.'L') then
          amass = 6.941
        else if (aty(1:1).eq.'S') then
          amass = 32.060
        else if (aty(1:1).eq.'F') then
          amass = 55.847
        else
          write(*,*)'Atom type not defined!   coord - line:',i
          stop
        end if
        do j=2,0,-1
          sqr_mass(3*i-j) = sqrt(amass)
        end do
      end do
      close(9)
      return
      end

      subroutine cross(a,b,pfc,c)
      implicit double precision(a-h,o-z)
      dimension a(3),b(3),c(3)
      c(1) = pfc*(a(2)*b(3) - b(2)*a(3))
      c(2) = pfc*(a(3)*b(1) - b(3)*a(1))
      c(3) = pfc*(a(1)*b(2) - b(1)*a(2)) 
      return
      end

      subroutine mix(a,b,c,v)
      implicit double precision(a-h,o-z)
      dimension a(3),b(3),c(3),d(3)
      call cross(a,b,1.d0,d)
      v = c(1)*d(1) + c(2)*d(2) + c(3)*d(3)
      return
      end


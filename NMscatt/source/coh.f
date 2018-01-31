c COH module (coherent scattering)
c of the NMscatt program package
c Author: Franci Merzel
c         
      Program coherent_scatt
      implicit none
      integer nmax,mxat,ii,i0,i1,i2,nelm,nat,i,j,mfr,ifr,npt,ls,fbz,mz
      integer neigh,int,mip,nip,k,k1,k2,j1,j2,l,nkp,nfr,mqp,ixf,ibz
      parameter(mxat=3600,npt=3000,mip=10,mfr=20,mqp=10)
      parameter(nmax=3*mxat)
      integer jev(0:mfr),nat_o,ntfr,nk,mk,nbz1,nbz2
      real*4 mass(mxat),cs(mxat),x(mxat),y(mxat),z(mxat),wgh(nmax),eps
      real*4 fct(0:mfr,3),pi,ev(0:mfr,nmax),skpp,skpn,fmtp,qn(0:mfr,3)
      real*4 ai,bi,xm,xr,xi(mip),wi(mip),xc,wc,ss,alph,bet,ff,q0n
      real*4 wint,aa,bb,kbt,ang,qp(0:mfr,3),q0p,phi,omg(npt),cc,ss3
      real*4 hbar,dwfp(mxat,0:mfr),dotqe,ss1,ss_n,re,wmx(3),kap(3),debye
      real*4 sqw(npt,-mfr:mfr),distf(npt),qavs(mxat,npt),ang0,ggg(3)
      real*4 frqmx,ww,xy,ssn,www,qq,rrr,ss_n1,sq(3),rrp(mxat),temp0
      real*4 enf,akf,aki,alp,frqmx1,sqw1(npt,-mfr:mfr),offst,sev(3)
      real*4 ffp,ffp1,csp(mxat),csn(mxat),dwfn(mxat,0:mfr),ggg0(3)
      real*4 brva(3),brvb(3),brvc(3),vlc,rcpa(3),rcpb(3),rcpc(3),pfc
      real*4 ffn,ffn1
      real*4 evd,ccd,ggd,xx(nmax),yy(nmax)

      complex dotp,dotn,rr(mxat,3,3),bmat(mxat,3,3),ssc,vct(3)
      complex pom,pp1,pp2,cdsfp,cdsfp1,cdsfn,cdsfn1,img

      complex evct(nmax,nmax),vc1,vc2,dotd,p1,p2

      character*1 ch,at(mxat),sg(mxat),sct
      character*2 ch2,cbz
      character*4 arg      
      character*11 name_crd

      logical deut

      pi  = 4.0*atan(1.0)

c ####################################################################
c #### DATA given at the command line ################################
      call getarg(1,arg)    ! max k-point suffix n according to files
      read(arg,*)nfr        ! eig_val/vec_n, n = 0..nfr 0: Gamma-point
                            ! for evaluation of DW factors

      call getarg(2,arg)    ! k-point suffix indicating q sampling
      read(arg,*)nkp        ! range n=0..nkp, which fixes a q-direction;
                            ! each k_n must lie along the chosen q direction

      call getarg(3,arg)    ! scattering type (x-ray/neutrons)
      read(arg,*)sct

      call getarg(4,arg)    ! n corresponding to k-point defining
      read(arg,*)ixf        ! the linear regime (0-n) of acoustic branch

c experimental data 
      temp0  = 10.0         ! temperature [K]
      deut   = .false.      ! logical flag for deuteration

      frqmx1 = 2600.0       ! freq. range for S(q,w): 0-frqmx1

c Bravais lattice vectors
      data brva/ 15.00,  0.00,  0.00/
      data brvb/  0.00, 15.00,  0.00/
      data brvc/  0.00,  0.00, 20.25/

      call mix(brva,brvb,brvc,vlc)
      pfc = 2.*pi/vlc

c reciprocal vectors
      call cross(brvb,brvc,pfc,rcpa)
      call cross(brvc,brva,pfc,rcpb)
      call cross(brva,brvb,pfc,rcpc)

c reciprocal lattice vector G assigning a
c direction for sampling the momentum transfer vector q
      do k=1,3
        ggg0(k) = rcpc(k)
      end do

      nbz1 = 0             ! nbz1..nbz2: Brillouin zones defining
      nbz2 = 1             ! the sampling interval of the momentum
                           ! transfer vector q
c numerical constants
      nip = 9              ! Gaussian integration points <=mip

      name_crd = 'coord'   ! CHARMM coordinate file 
                           ! (for assigning atom types)
c ####################################################################
      
      write(*,'(a,3f12.5)')'rcp_v1',rcpa(1),rcpa(2),rcpa(3)
      write(*,'(a,3f12.5)')'rcp_v2',rcpb(1),rcpb(2),rcpb(3)
      write(*,'(a,3f12.5)')'rcp_v3',rcpc(1),rcpc(2),rcpc(3)

c cm-1 units used for energy terms!
      kbt  = temp0*0.695028
      hbar   = 33.4664       ! hbar^2/at.m.u [cm-1 A2]

      if (sct.eq.'x') then
        write(*,*)'INELASTIC COHERENT X-ray SCATTERING'
      else
        write(*,*)'INELASTIC COHERENT neutron SCATTERING'
      end if
      write(*,'(a,f8.1)')'TEMPERATURE:',temp0
      write(*,*)
      write(*,*)'Reciprocal lattice vectors:'
      write(*,*)
      write(*,'(a,3f12.5)')'rcp_v1',rcpa(1),rcpa(2),rcpa(3)
      write(*,'(a,3f12.5)')'rcp_v2',rcpb(1),rcpb(2),rcpb(3)
      write(*,'(a,3f12.5)')'rcp_v3',rcpc(1),rcpc(2),rcpc(3)
      write(*,*)
      img = cmplx(0.0,1.0)
      call assign_atoms(nat,mass,cs,sg,at,x,y,z,deut,sct,name_crd)

      call gauss(nip,nip,xi,wi)

      nelm = 3*nat
      write(*,*)'nelm:',nelm

      do ifr=0,nfr
        if (ifr.lt.10) then
          write(ch,'(i1)')ifr
          open(3,file='eig_val_'//ch,form='formatted',
     1         status ='unknown')
        else
          write(ch2,'(i2)')ifr
          open(3,file='eig_val_'//ch2,form='formatted',
     1         status ='unknown')
        end if
        read(3,'(28x,3f8.4)')(fct(ifr,k),k=1,3)
        write(*,'(a,3f7.3)')'fct:',(fct(ifr,k),k=1,3)

        if ((fct(ifr,1).gt.0.5).or.(fct(ifr,2).gt.0.5).or.
     1      (fct(ifr,3).gt.0.5)) then
          write(*,*)'k-point outside of the 1st Brillouin zone!'
          stop
        end if
        read(3,*)
        do j=1,10
          read(3,*)ii,xy
          if (xy.le.0.0) then
            write(*,*)'Negative frequencies detected in frame:',ifr
            stop
          end if
        end do
        close(3)
      end do

      frqmx = 0.0
      do ifr=0,nfr
        if (ifr.lt.10) then
          write(ch,'(i1)')ifr
          open(3,file='eig_val_'//ch,form='formatted',
     1         status ='unknown')
        else
          write(ch2,'(i2)')ifr
          open(3,file='eig_val_'//ch2,form='formatted',
     1         status ='unknown')
        end if

        read(3,*)
        read(3,*)
        do j=1,nelm
          read(3,*)ii,xy,ev(ifr,j)
        end do
        if (ev(ifr,nelm).gt.frqmx) frqmx = ev(ifr,nelm)
        close(3)
        if (ifr.eq.ixf) then
          write(*,*)'acoustic frequencies at Xa-point:'
          do j=1,3         ! acoustic frequencies
            wmx(j) = ev(ifr,j)
            write(*,'(i2,f10.5)')j,wmx(j)
          end do
        end if
      end do
      frqmx = frqmx + 100.0

      if (frqmx1.gt.frqmx) frqmx1 = frqmx     

      wint = frqmx/float(npt)
      do k=1,npt
        omg(k) = wint*(float(k)-0.5)
      end do

      do i=1,nat
        do j=1,3
          do k=1,3
            bmat(i,j,k) = cmplx(0.0,0.0)
            rr(i,j,k) = cmplx(0.0,0.0)
          end do
        end do
      end do

c evaluating Debye function
      do i=1,3
        kap(i) = debye(wmx(i)/kbt,nip,xi,wi)
      end do

 100  continue
      do ifr=0,nfr ! loop over phonon k-points

        if (ifr.lt.10) then
          write(ch,'(i1)')ifr
          open(2,file='eig_vec_'//ch,form='unformatted',
     1         status ='unknown')
        else
          write(ch2,'(i2)')ifr
          open(2,file='eig_vec_'//ch2,form='unformatted',
     1         status ='unknown')
        end if

        write(*,*)'Reading eigen-vectors....'
        write(*,*)
        ii = 0
        i0 = 0
        do ii=1,nelm
          read(2)(xx(j),yy(j),j=1,nelm)
          do j=1,nelm
            evct(ii,j) = cmplx(xx(j),yy(j))
          end do  
          if (i0.ne.ii) then
            if (mod(ii,100).eq.0) write(*,'(a,i7,a,i2,a,i2)')'Passed:',
     1         ii,' lines   -  ',ifr+1,' / ',nfr+1
            i0 = ii
          end if
        end do
        close(2)

        if (nelm/3.ne.nat) then
          write(*,*)'<coord> and <eig-vect> mismatch!!'
          stop
        end if

        do j=1,nelm
          if (mod(j,500).eq.0) write(*,'(i5,a,i5)')j,' / ',nelm
          evd = ev(ifr,j)
          ccd = 0.5d0*evd/kbt

          if (ifr.le.ixf.and.j.le.3) then
            if (ifr.lt.ixf) then
              ggd = 0.0
            else
              ggd = 3.0*kap(j)/wmx(j)
            end if
          else
            if (ccd.gt.15.0) then
              ggd = 1.0d0/evd
            else
              ggd = (exp(ccd)+exp(-ccd))/(exp(ccd)-exp(-ccd))/evd
            end if
          end if

          do i=1,nat
            do k1=1,3
              do k2=1,3
                j1 = 3*(i-1) + k1
                j2 = 3*(i-1) + k2
                vc1 = evct(j1,j)
                vc2 = conjg(evct(j2,j))
                dotd = vc1*vc2
                p1 = dotd
                p2 = dotd*ggd
                pp1 = p1
                pp2 = p2
                rr(i,k1,k2)   = rr(i,k1,k2)   + pp1/
     1                          float(nfr+1)/float(nelm)
                bmat(i,k1,k2) = bmat(i,k1,k2) + pp2/
     1                          float(nfr+1)/float(nelm)
              end do
            end do
          end do
        end do
      end do

c normalization of partial g(w) --> 1/r
      rrr = float(nelm)
      do i=1,nat
        do j=1,3
          rrp(i) = rrp(i) + re(rr(i,j,j))/3.0
        end do
        rrp(i) = 1.0/rrp(i)
        if (abs(rrr-rrp(i)).gt.0.2) 
     1    write(*,'(i6,3f12.3)')i,rrr,rrp(i),rrr-rrp(i)
      end do

c ---------------------------------------------------------------------
      write(*,*)
      write(*,'(a,f6.2,a,f6.2)')'Default values: T = ',
     1      temp0

c B matrix
      do i=1,nat
        do k1=1,3
          do k2=k1,3
            pom = 0.5*hbar*bmat(i,k1,k2)*rrr/mass(i)
            bmat(i,k1,k2) = pom
            if (k1.ne.k2) bmat(i,k2,k1) = conjg(bmat(i,k1,k2))
          end do
        end do
      end do

c w-range is modified here w:(o..frqmx1)
      wint = frqmx1/float(npt)
      do k=1,npt
        omg(k) = wint*(float(k)-0.5)
      end do

c      write(*,*)'Smoothing factor for frequency peaks:'
c      write(*,*)'bet > 1.0: sharper, bet < 1.0: smoother;    bet = ?'
c      read*,bet
      bet = 0.25
      eps = 1.e-4
      alph = (4.0/(3.0*wint))**2*bet

      ww = sqrt(log(sqrt(pi/alph)/eps)/alph)
      neigh = nint(ww/wint)
      if (neigh.lt.0) neigh = 1

c*******************************************************
c******************  S(q,w) part  **********************
c*******************************************************

      do ibz=nbz1,nbz2     ! loop over BZ
        if (ibz.lt.10) then
          write(cbz,'(i1)')ibz
        else
          write(cbz,'(i2)')ibz
        end if
        do k=1,3
          ggg(k) = float(ibz)*ggg0(k)  ! needs to be chosen along 
        end do                         ! the dispersion direction
      do ifr=0,nkp ! loop over phonon k-points

c momentum transfer vector
        do k=1,3
          qp(ifr,k) = fct(ifr,1)*rcpa(k) + fct(ifr,2)*rcpb(k) + 
     1                fct(ifr,3)*rcpc(k) + ggg(k)
          qn(ifr,k) =-fct(ifr,1)*rcpa(k) - fct(ifr,2)*rcpb(k) - 
     1                fct(ifr,3)*rcpc(k) + ggg(k)
        end do
c Debye-Waller for the corresponding q
        do j=1,nat
          do k=1,3
            vct(k) = cmplx(0.0,0.0)
            do l=1,3
              vct(k) = vct(k) + bmat(j,k,l)*qp(ifr,l)
            end do
          end do
          ssc = cmplx(0.0,0.0)
          do k=1,3
            ssc = ssc + vct(k)*qp(ifr,k)
          end do
          dwfp(j,ifr) = 0.5*re(ssc)

          do k=1,3
            vct(k) = cmplx(0.0,0.0)
            do l=1,3
              vct(k) = vct(k) + bmat(j,k,l)*qn(ifr,l)
            end do
          end do
          ssc = cmplx(0.0,0.0)
          do k=1,3
            ssc = ssc + vct(k)*qn(ifr,k)
          end do
          dwfn(j,ifr) = 0.5*re(ssc)  !*4.0*pi**2
        end do

c reset variables
        do k=1,npt
          sqw(k,  ifr+1) = 0.0
          sqw1(k, ifr+1) = 0.0
          sqw(k, -ifr-1) = 0.0
          sqw1(k,-ifr-1) = 0.0
        end do

        q0p = sqrt(qp(ifr,1)**2+qp(ifr,2)**2+qp(ifr,3)**2)
        q0n = sqrt(qn(ifr,1)**2+qn(ifr,2)**2+qn(ifr,3)**2)
        write(*,'(a,2f12.5)')'q0p, q0n:',q0p,q0n
        do i=1,nat
          csp(i) = cs(i)
          csn(i) = cs(i)
        end do
c re-assigning of scattering lengths for X-rays
        if (sct.eq.'x') then
          write(*,*)
          call xsctlng(nat,q0p,csp)
          call xsctlng(nat,q0n,csn)
cc test output (first 10 atoms)
c          do i=1,10
c             write(*,'(i3,a,2x,2f7.3,2f10.4)')i,at(i),q0p,q0n,
c     1                                        csp(i),csn(i)
c          end do
c          write(*,*)
        end if

        if (nfr.ne.0) then
        if (ifr.lt.10) then
          write(ch,'(i1)')ifr
          open(2,file='eig_vec_'//ch,form='unformatted',
     1         status ='unknown')
        else
          write(ch2,'(i2)')ifr
          open(2,file='eig_vec_'//ch2,form='unformatted',
     1         status ='unknown')
        end if

        write(*,*)'Reading eigen-vectors....'
        write(*,*)
        ii = 0
        i0 = 0
        do ii=1,nelm
          read(2)(xx(j),yy(j),j=1,nelm)
          do j=1,nelm
            evct(ii,j) = cmplx(xx(j),yy(j))
          end do  
          if (i0.ne.ii) then
            if (mod(ii,100).eq.0) write(*,'(a,i7,a,i2,a,i2)')'Passed:',
     1         ii,' lines   -  ',ifr+1,' / ',nfr+1
            i0 = ii
          end if
        end do
        close(2)
        end if

        do k=1,npt
          if (mod(k,100).eq.0) 
     1    write(*,'(a,i5,a,i5,a,i3,a,i3,a,i3)')'S(q,w): ',k,'/',npt,
     2            '  k:',ifr,'  BZ: ',ibz,'/',nbz2
          ai = omg(k)-0.5*wint
          bi = omg(k)+0.5*wint
          xm = 0.5*(ai+bi)
          xr = 0.5*(bi-ai)

          ffp  = 0.0
          ffp1 = 0.0
          ffn  = 0.0
          ffn1 = 0.0
          do int=1,nip
            xc = xm + xr*xi(int)
            wc = wi(int)*xr
            do j=1,nelm
              if (ev(ifr,j).le.frqmx1) then
                cc = ev(ifr,j)/kbt
                aa = abs(ev(ifr,j)-omg(k))
                bb = (float(neigh)+0.5)*wint
                if (aa.lt.bb) then
                  qq = sqrt(alph/pi)*exp(-alph*(ev(ifr,j)-xc)**2)*
     1                 0.5*hbar/(1.0-exp(-cc))/ev(ifr,j)

c coherent dynamical structor factor
                  cdsfp  = cmplx(0.0,0.0)
                  cdsfp1 = cmplx(0.0,0.0)
                  cdsfn  = cmplx(0.0,0.0)
                  cdsfn1 = cmplx(0.0,0.0)
                  do i=1,nat
                    dotp = cmplx(0.0,0.0)
                    dotn = cmplx(0.0,0.0)
                    do l=1,3
                      j1 = 3*(i-1) + l
                      dotp = dotp + conjg(evct(j1,j))*qp(ifr,l)
                      dotn = dotn + evct(j1,j)*qn(ifr,l)
                    end do
                    skpp = x(i)*qp(ifr,1)+y(i)*qp(ifr,2)+z(i)*qp(ifr,3)
                    skpn = x(i)*qn(ifr,1)+y(i)*qn(ifr,2)+z(i)*qn(ifr,3)

                    cdsfp = cdsfp + csp(i)*dotp/sqrt(mass(i))*
     1                              exp(-dwfp(i,ifr) + img*skpp)
                    cdsfn = cdsfn + csn(i)*dotn/sqrt(mass(i))*
     1                              exp(-dwfn(i,ifr) + img*skpn)
                    if (sg(i).eq.'B') then
                      cdsfp1 = cdsfp1 + csp(i)*dotp/sqrt(mass(i))*
     1                                  exp(-dwfp(i,ifr) + img*skpp)
                      cdsfn1 = cdsfn1 + csn(i)*dotn/sqrt(mass(i))*
     1                                  exp(-dwfn(i,ifr) + img*skpn)
                    end if
                  end do
                  ffp  = ffp  + qq*wc * re(cdsfp *conjg(cdsfp))
                  ffp1 = ffp1 + qq*wc * re(cdsfp1*conjg(cdsfp1))
                  ffn  = ffn  + qq*wc * re(cdsfn *conjg(cdsfn))
                  ffn1 = ffn1 + qq*wc * re(cdsfn1*conjg(cdsfn1))
                end if
              end if
            end do
          end do
          sqw(k ,ifr+1)  = sqw(k ,ifr+1)  + ffp/wint
          sqw1(k,ifr+1)  = sqw1(k,ifr+1)  + ffp1/wint

          sqw(k ,-ifr-1) = sqw(k ,-ifr-1) + ffn/wint
          sqw1(k,-ifr-1) = sqw1(k,-ifr-1) + ffn1/wint
        end do
      end do !ifr

c*******************************************************
      if (ibz.lt.10) then
        open(8,file='Sqw_coh'//cbz(1:1),form='formatted',
     1         status ='unknown')
        open(9,file='Sqw_coh_r'//cbz(1:1),form='formatted',
     1         status ='unknown')
      else
        open(8,file='Sqw_coh'//cbz(1:2),form='formatted',
     1         status ='unknown')
        open(9,file='Sqw_coh_r'//cbz(1:2),form='formatted',
     1         status ='unknown')
      end if

      write(8,'(a,f8.3,2i3)')'# temp, n k-pts, n q-pts',
     1                    temp0,(nfr+1),2*nkp+1
      write(9,'(a,f8.3,2i3)')'# temp, n k-pts, n q-pts',
     1                    temp0,(nfr+1),2*nkp+1
      nk = (nkp+1)/3
      mk = mod(nkp+1,3)

      do i=1,nk
        write(8,'(a,3(5x,3f6.3,a))')'#',((qn(nkp-ii+1,k),k=1,3),';',
     1        ii=3*(i-1)+1,3*(i-1)+3)
        write(9,'(a,3(5x,3f6.3,a))')'#',((qn(nkp-ii+1,k),k=1,3),';',
     1        ii=3*(i-1)+1,3*(i-1)+3)
      end do
      write(8,'(a,3(5x,3f6.3,a))')'#',((qn(nkp-ii+1,k),k=1,3),';',
     1        ii=3*nk+1,3*nk+mk)
      write(9,'(a,3(5x,3f6.3,a))')'#',((qn(nkp-ii+1,k),k=1,3),';',
     1        ii=3*nk+1,3*nk+mk)

      do i=1,nk
        write(8,'(a,3(5x,3f6.3,a))')'#',((qp(ii-1,k),k=1,3),';',
     1        ii=3*(i-1)+1,3*(i-1)+3)
        write(9,'(a,3(5x,3f6.3,a))')'#',((qp(ii-1,k),k=1,3),';',
     1        ii=3*(i-1)+1,3*(i-1)+3)
      end do
      write(8,'(a,3(5x,3f6.3,a))')'#',((qp(ii-1,k),k=1,3),';',
     1        ii=3*nk+1,3*nk+mk)
      write(9,'(a,3(5x,3f6.3,a))')'#',((qp(ii-1,k),k=1,3),';',
     1        ii=3*nk+1,3*nk+mk)

      fmtp = 100.0  ! constant for scaling output values 
      do k=1,npt
        write(8,'(f10.3,21f12.3)')omg(k),
     1       (fmtp*sqw(k,ifr), ifr=-nkp-1,-1),
     2       (fmtp*sqw(k,ifr), ifr=1,nkp+1)
        write(9,'(f10.3,21f12.3)')omg(k),
     1       (fmtp*sqw1(k,ifr),ifr=-nkp-1,-1),
     2       (fmtp*sqw1(k,ifr),ifr=1,nkp+1)
      end do
      close(8)
      close(9)
      end do  ! ibz
      end

      subroutine assign_atoms(nat,mass,cs,sg,at,x,y,z,deut,sct,name_crd)
      implicit none
      integer nmax,nat,i,ii,ir,mxat,k
      parameter(mxat=3600)
      parameter(nmax=3*mxat)
      real*4 mass(mxat),cs(mxat),x(mxat),y(mxat),z(mxat),amass,acs,acsx
      character*4 aty,seg
      character*6 ars
      character*1 at(mxat),sg(mxat),sct
      character*11 name_crd

      logical wt,oh,hn,ht,hr,r1,r2,r3,r4,r5,r6,r7,r8,exch,deut

      open(9, file=name_crd, form='formatted', status ='unknown')
c  suppose 3 lines off-set
      do k=1,3
        read(9,*)
      end do

      read(9,*)nat
      do i=1,nat
        read(9,'(2i5,a6,a4,3f10.5,1x,a4)')ii,ir,ars,aty,x(i),y(i),z(i),
     1                                    seg
c sorting atom types        
        if (aty(1:1).eq.'H') then
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
c            acs = 2.04   !incoh 
            acs = 7.64   !scatt
          else
c hydrogen
            amass =  1.008
c            acs = 79.90  !incoh
            acs = 81.66  !scatt            
          end if
          acsx = 1.0
        else if (aty(1:1).eq.'C') then
          if (aty(2:3).eq.'LA') then
            amass = 35.45
c            acs = 5.2  !incoh
            acs = 16.7  !scatt
            acsx = 6.0
          else
            amass = 12.011
c            acs = 0.001  !incoh
            acs = 5.555  !scatt
            acsx = 17.0
          end if
        else if (aty(1:1).eq.'N') then
          amass = 14.007
c          acs = 0.49   !incoh
          acs = 11.50  !scatt
          acsx = 7.0
        else if (aty(1:1).eq.'O') then
          amass = 15.9994
c          acs = 0.00   !incoh
          acs = 4.24   !scatt
          acsx = 8.0
        else if (aty(1:1).eq.'P') then
          amass = 30.974
c          acs = 0.01   !incoh
          acs = 3.31   !scatt
          acsx = 15.0
        else if (aty(1:1).eq.'S') then
          amass = 32.07
c          acs = 0.01   !incoh
          acs = 1.03   !scatt
          acsx = 16.0
        else if (aty(1:1).eq.'F') then
          amass = 55.85
c          acs = 0.39   !incoh
          acs = 11.83   !scatt
          acsx = 26.0
        else 
          write(*,*)'Atom type not defined!   coord - line:',i
          stop
        end if
        mass(i) = amass
        if (sct.eq.'x') then
          cs(i) = acsx
        else
          cs(i) = acs
        end if
        at(i) = aty(1:1)
        sg(i) = seg(1:1)
c assuming the initial of the segment name of water being 'W'
        if (sg(i).ne.'W') sg(i)='B' ! select biomolecule
      end do
      close(9)
      return
      end

      subroutine gauss(nn,n,x,w)
c  weights and coordinates for Gausian quadrature of level n
      implicit none
      integer i,j,m,n,nn
      real*4 eps,p1,p2,p3,pp,x1,x2,xl,xm,z,z1,x,w
      parameter(eps=3.e-6)
      dimension x(nn),w(nn)
      x1=-1.0
      x2= 1.0
      m=(n+1)/2
      xm=0.5*(x2+x1)
      xl=0.5*(x2-x1)
      do i=1,m
        z=cos(3.141592654*(i-0.25)/(n+0.5))
1       continue
          p1=1.0
          p2=0.0
          do j=1,n
            p3=p2
            p2=p1
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
          end do
          pp=n*(z*p1-p2)/(z*z-1.0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.eps)go to 1
        x(i)=xm-xl*z
        x(n-i+1)=xm+xl*z
        w(i)=2.0*xl/((1.0-z*z)*pp*pp)
        w(n-i+1)=w(i)
      end do
      return
      end

      real*4 function re(cmplx)
      complex cmplx
      re = 0.5*(cmplx+conjg(cmplx))
      return
      end

      real*4 function debye(xx,nip,xi,wi)
      implicit none
      integer mip,nip,i
      parameter(mip=10)
      real*4 xx,xi,wi,ss,xh,xc,wc,ff
      dimension xi(mip),wi(mip)

      ss = 0.0
      xh = 0.5*xx
      do i=1,nip
        xc = xh + xh*xi(i)
        wc = wi(i)*xh
        ff = xc/(exp(xc)-1.0)/xx**2
        ss = ss + wc*ff
      end do
      debye = ss + 0.25
      return
      end

c electron density for isolated atoms according to the Slater rules
      real*4 function den(z,r)
      implicit none
      integer i,j,k,z
      real*4 a,pi,r,c,fak,coef(18,3,5)
      data (((coef(i,j,k),k=1,5),j=1,3),i=1,18)
     */ 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
     *  1.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  1.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  1.6875, 0.0000, 0.0000, 0.0000, 0.0000, 
     *  2.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  1.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  2.6906, 1.2792, 0.0000, 0.0000, 0.0000, 
     *  2.0000, 1.0000, 0.0000, 0.0000, 0.0000,
     *  1.0000, 2.0000, 0.0000, 0.0000, 0.0000,
     *  3.6848, 1.9120, 0.0000, 0.0000, 0.0000, 
     *  2.0000, 2.0000, 0.0000, 0.0000, 0.0000,
     *  1.0000, 2.0000, 0.0000, 0.0000, 0.0000,
     *  4.6795, 2.5762, 2.4214, 0.0000, 0.0000, 
     *  2.0000, 2.0000, 1.0000, 0.0000, 0.0000,
     *  1.0000, 2.0000, 2.0000, 0.0000, 0.0000,
     *  5.6727, 3.2166, 3.1358, 0.0000, 0.0000, 
     *  2.0000, 2.0000, 2.0000, 0.0000, 0.0000,
     *  1.0000, 2.0000, 2.0000, 0.0000, 0.0000,
     *  6.6651, 3.8474, 3.8340, 0.0000, 0.0000, 
     *  2.0000, 2.0000, 3.0000, 0.0000, 0.0000,
     *  1.0000, 2.0000, 2.0000, 0.0000, 0.0000,
     *  7.6579, 4.4916, 4.4532, 0.0000, 0.0000, 
     *  2.0000, 2.0000, 4.0000, 0.0000, 0.0000,
     *  1.0000, 2.0000, 2.0000, 0.0000, 0.0000,
     *  8.6501, 5.1276, 5.1000, 0.0000, 0.0000, 
     *  2.0000, 2.0000, 5.0000, 0.0000, 0.0000,
     *  1.0000, 2.0000, 2.0000, 0.0000, 0.0000,
     *  9.6421, 5.7584, 5.7584, 0.0000, 0.0000, 
     *  2.0000, 2.0000, 6.0000, 0.0000, 0.0000,
     *  1.0000, 2.0000, 2.0000, 0.0000, 0.0000,
     * 10.6259, 6.5714, 6.8018, 2.5074, 0.0000, 
     *  2.0000, 2.0000, 6.0000, 1.0000, 0.0000,
     *  1.0000, 2.0000, 2.0000, 3.0000, 0.0000,
     * 11.6089, 7.3920, 7.8258, 3.3075, 0.0000, 
     *  2.0000, 2.0000, 6.0000, 2.0000, 0.0000,
     *  1.0000, 2.0000, 2.0000, 3.0000, 0.0000,
     * 12.5910, 8.2136, 8.9634, 4.1172, 4.0656, 
     *  2.0000, 2.0000, 6.0000, 2.0000, 1.0000,
     *  1.0000, 2.0000, 2.0000, 3.0000, 3.0000,
     * 13.5745, 9.0200, 9.9450, 4.9032, 4.2852, 
     *  2.0000, 2.0000, 6.0000, 2.0000, 2.0000,
     *  1.0000, 2.0000, 2.0000, 3.0000, 3.0000,
     * 14.5578, 9.8250,10.9612, 5.6418, 4.8864, 
     *  2.0000, 2.0000, 6.0000, 2.0000, 3.0000,
     *  1.0000, 2.0000, 2.0000, 3.0000, 3.0000,
     * 15.5409,10.6288,11.9770, 6.3669, 5.4819, 
     *  2.0000, 2.0000, 6.0000, 2.0000, 4.0000,
     *  1.0000, 2.0000, 2.0000, 3.0000, 3.0000,
     * 16.5239,11.4304,12.9932, 7.0683, 6.1161, 
     *  2.0000, 2.0000, 6.0000, 2.0000, 5.0000,
     *  1.0000, 2.0000, 2.0000, 3.0000, 3.0000,
     * 17.5075,12.2304,14.0082, 7.7568, 6.7641, 
     *  2.0000, 2.0000, 6.0000, 2.0000, 6.0000,
     *  1.0000, 2.0000, 2.0000, 3.0000, 3.0000/
      a = 1.0 
      pi=4.*atan(1.0)
      den=0.0
      do j=1,5
        if (coef(z,1,j).ne.0.0) then
          c = ((2.0*coef(z,1,j)/coef(z,3,j))**(coef(z,3,j)+0.5))/
     1        sqrt(fak(int(2.*coef(z,3,j))))
          den = den+coef(z,2,j)*(c*r**(coef(z,3,j)-1.0)*
     1          exp(-coef(z,1,j)*r/coef(z,3,j)/a))**2
        end if
      end do
      den = den/4.0/pi
      return
      end

      real*4 function fak(n)
      implicit none
      integer i,n
      fak = 1.0
      do i=2,n
        fak = fak*dfloat(i)
      end do
      end

      real*4 function func(q,iz,x)
      implicit none
      integer iz
      real*4 den,q,pi,x
      pi = 4.*atan(1.0)
      func = 4.0*pi*x**2*den(iz,x)*sin(q*x)/(q*x)
      return
      end

      subroutine radinteg(q,scl)
      implicit none
      integer nm,mz,nz,iz,nn,mm,nn1,nn2,j
      parameter(nm=50,mz=18)
      real*4 q,xi1(nm),wi1(nm),xi2(nm),wi2(nm),a1,a2,a3,ss
      real*4 scl(mz),xm,xr,dx,func

      mm = nm
      a1 =  0.0
      a2 = 10.0
      a3 = 30.0
      nn1 = 50
      nn2 = 20
      call gauss(mm,nn1,xi1,wi1)
      call gauss(mm,nn2,xi2,wi2)
      do iz=1,mz
        xm =0.5*(a1+a2)
        xr =0.5*(a2-a1)
        ss = 0.0
        do j=1,nn1
          dx = xr*xi1(j)
          ss = ss + wi1(j)*func(q,iz,xm+dx)
        end do
        scl(iz) = ss*xr

        xm =0.5*(a2+a3)
        xr =0.5*(a3-a2)
        ss = 0.0
        do j=1,nn2
          dx = xr*xi2(j)
          ss = ss + wi2(j)*func(q,iz,xm+dx)
        end do
        scl(iz) = scl(iz) + ss*xr
      end do
      return
      end

      subroutine xsctlng(nat,q,cs)
      implicit none
      integer mxat,mz,nat,i,j
      parameter(mxat=3600,mz=18)
      real*4 q,scl(mz),cs(mxat)

      call radinteg(q,scl)
      do i=1,nat
        do j=1,mz
          if (nint(cs(i)).eq.j) cs(i) = scl(j)
        end do
      end do
      return
      end

      subroutine cross(a,b,pfc,c)
      implicit real*4(a-h,o-z)
      dimension a(3),b(3),c(3)
      c(1) = pfc*(a(2)*b(3) - b(2)*a(3))
      c(2) = pfc*(a(3)*b(1) - b(3)*a(1))
      c(3) = pfc*(a(1)*b(2) - b(1)*a(2)) 
      return
      end

      subroutine mix(a,b,c,v)
      implicit real*4(a-h,o-z)
      dimension a(3),b(3),c(3),d(3)
      call cross(a,b,1.0,d)
      v = c(1)*d(1) + c(2)*d(2) + c(3)*d(3)
      return
      end


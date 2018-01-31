c BEAD module (vibrational modes analysis)
c of the NMscatt program package
c Author: Franci Merzel
c
      Program bead_display
      implicit none
      integer nmax,mxat,ii,im,nelm,nat,i,ix,iy,iz,ism,nat1,j,k,mres,nres
      parameter(mxat=3600,mres=1000)
      parameter(nmax=3*mxat)
      integer bead(mxat),nabd(mres),ir,ib
      real*4 ev,x(mxat),y(mxat),z(mxat),xx(2*nmax),scl
      real*4 amass(mxat),scl0,cmsb(3,mres),disp(3,mres),ambd(mres)

      character*1  ch1,rsd(mxat),bd(mres)
      character*2  ch,at(mxat)
      character*5  arg
      character*11 name_crd

      common/par/ev(nmax,nmax)

c ####################################################################
c #### DATA given at the command line ################################
      call getarg(1,arg) ! 5 sequential normal modes to display: 
      read(arg,*)ism     ! starting with ism: (ism, ism+1, ism+2...)
c
      ch = '0'           ! suffix of the eigen-vector file used for 
                         ! vizualization of displacements (Gamma-point)

      nat1 = 245         ! number of protein atoms
      nres = 21          ! number of residues in the protein part

      scl0 = 100.0       ! scaling factor for displacement vector

      name_crd = 'coord' ! CHARMM coordinate file 
                         ! (for assigning atom types)
c ####################################################################

      call assig_atyp(nat,at,bead,rsd,amass,x,y,z,name_crd)

      write(*,*)'nat:',nat

      nelm = 3*nat
      open(2,file='eig_vec_'//ch,form='unformatted',status ='unknown')
      write(*,*)'Reading eigen-vectors....'
      write(*,*)
      do ii=1,nelm
        read(2)(xx(j),j=1,2*nelm)
c taking into account only real part
        do j=1,nelm
          ev(ii,j) = xx(2*j-1)
        end do
        if (mod(ii,100).eq.0) write(*,'(a,i7,a)')'Passed:',ii,' lines'
      end do
      close(2)

      do im=ism,ism+4
        write(ch1,'(i1)')im-ism+1
        write(*,'(a,i5,a,i5,a,i2)')
     1          'Writing mode no. ',im,' / ',nelm,' to file ',im-ism+1

        open(3,file='atm'//ch1//'.xyz',form='formatted',
     1                                 status ='unknown')
        write(3,*)nat1
        write(3,*)
        do i=1,nat1
          ix = 3*i-2
          iy = 3*i-1
          iz = 3*i-0
          scl = scl0/sqrt(amass(i))
          write(3,'(a,3f12.7,a,3f12.7,3x,a,i4)')at(i)(1:2),
     1        x(i),y(i),z(i),'  atom_vector ',
     2        ev(ix,im)*scl,ev(iy,im)*scl,ev(iz,im)*scl
        end do

        open(7,file='bead'//ch1//'.xyz',form='formatted',
     1                                 status ='unknown')
        write(7,*)nres
        write(7,*)

        do ir=1,nres
          nabd(ir) = 0
          ambd(ir) = 0.0
          do k=1,3
            cmsb(k,ir) = 0.0
            disp(k,ir) = 0.0
          end do
        end do

        do ir=1,nres
          do i=1,nat1
            if (bead(i).eq.ir) then
              bd(ir) = rsd(i)
              nabd(ir) = nabd(ir) + 1 
              ambd(ir) = ambd(ir) + amass(i)
c center of mass
              cmsb(1,ir) = cmsb(1,ir) + amass(i)*x(i)
              cmsb(2,ir) = cmsb(2,ir) + amass(i)*y(i)
              cmsb(3,ir) = cmsb(3,ir) + amass(i)*z(i)
c displacement
              scl = scl0/sqrt(amass(i))
              disp(1,ir) = disp(1,ir) + ev(3*i-2,im)*scl
              disp(2,ir) = disp(2,ir) + ev(3*i-1,im)*scl
              disp(3,ir) = disp(3,ir) + ev(3*i-0,im)*scl
            end if
          end do
        end do

c normalization
        do ir=1,nres
          do k=1,3
            cmsb(k,ir) = cmsb(k,ir)/ambd(ir)
            disp(k,ir) = disp(k,ir)/float(nabd(ir))
          end do
        end do

        do ir=1,nres
          ib = ((ir-1)/7) + 1
          if (ib.eq.1) bd(ir) = 'A'
          if (ib.eq.2) bd(ir) = 'B'
          if (ib.eq.3) bd(ir) = 'C'
          write(7,'(a,2x,3f12.7,a,3f12.7)')bd(ir),
     1    cmsb(1,ir),cmsb(2,ir),cmsb(3,ir),'  atom_vector ',
     2    disp(1,ir),disp(2,ir),disp(3,ir)
        end do
      end do
      write(*,*)
      write(*,*)'Total number of modes:',nelm
      close(7)
      close(3)
      end

      subroutine assig_atyp(nat,at,bead,rsd,amass,x,y,z,name_crd)
      implicit none
      integer nat,i,ii,ir,j,k,mxat
      parameter(mxat=3600)
      integer bead(mxat)
      real*4 x(mxat),y(mxat),z(mxat)
      real*4 amass(mxat)
      character rsd(mxat)
      character*2 at(mxat)
      character*4 aty
      character*6 ars
      character*11 name_crd
      logical wt,ha,hb,hc,hd,he,hf,hg,hh,hi,hj,ht,exch

      open(9, file=name_crd, form ='formatted',status ='unknown')
c  suppose 4 lines off-set
      do k=1,3
        read(9,*)
      end do

      read(9,*)nat
      do i=1,nat
        read(9,'(2i5,a6,a4,3f10.5)')ii,ir,ars,aty,x(i),y(i),z(i)
        bead(i) = ir
        rsd(i) = ars(2:2)
c sorting atom types 
        if (aty(1:1).eq.'H') then
          wt = (ars(2:5).eq.'TIP3')   ! selecting water hydrogens
 
          ht = (aty(2:2).eq.'T')

          ha = (aty(2:3).eq.'N ')

          hb = (ars(2:4).eq.'ARG').and.
     1 ((aty(2:2).eq.'E').or.(aty(2:2).eq.'H'))

          hc = (ars(2:4).eq.'ASN').and.(aty(2:3).eq.'D2')

          hd = (ars(2:4).eq.'GLN').and.(aty(2:3).eq.'E2')

          he = (ars(2:4).eq.'HIS').and.
     1 ((aty(2:3).eq.'D1').or.(aty(2:3).eq.'E2'))

          hf = (ars(2:4).eq.'LYS').and.(aty(2:2).eq.'Z')

          hg = ((ars(2:4).eq.'SER').or.(ars(2:4).eq.'THR')).and.
     1          (aty(2:3).eq.'G1')

          hh = (ars(2:4).eq.'TRP').and.(aty(2:3).eq.'E1')

          hi = (ars(2:4).eq.'TYR').and.(aty(2:2).eq.'H')

          exch = (ha.or.hb.or.hc.or.hd.or.he.or.
     1            hf.or.hg.or.hh.or.hi.or.ht)
c deuterium 
          if (wt.or.exch) then
            amass(i) = 2.009
            at(i)(1:2) = 'D '
          else
            amass(i) =  1.008
            at(i)(1:2) = 'H '
          end if
        else if (aty(1:1).eq.'C') then
          amass(i) = 12.011
          at(i)(1:2) = 'C '
        else if (aty(1:1).eq.'N') then
          amass(i) = 14.007
          at(i)(1:2) = 'N '
        else if (aty(1:1).eq.'O') then
          amass(i) = 15.9994
          at(i)(1:2) = 'O '
        else if (aty(1:1).eq.'P') then
          amass(i) = 30.974
          at(i)(1:2) = 'P '
        else if (aty(1:1).eq.'S') then
          amass(i) = 32.07
          at(i)(1:2) = 'S '
        else if (aty(1:1).eq.'F') then
          amass(i) = 55.85
          at(i)(1:2) = 'FE'
        else
          write(*,*)'Atom type not defined!   coord - line:',i
          stop
        end if
      end do
      close(9)
      return
      end

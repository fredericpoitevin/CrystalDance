c INCOH module (incoherent scattering)
c of the NMscatt program package
c Author: Franci Merzel
c         

      Program incoherent_scatt
      implicit none
      integer nmax,mxat,ii,i0,i1,i2,nelm,nat,i,j,mfr,ifr,npt,ls
      integer neigh,int,mip,nip,k,k1,k2,j1,j2,l,nqp,nfr,mqp,ixf
      parameter(mxat=3600,npt=3000,mip=10,mfr=19,mqp=50)
      parameter(nmax=3*mxat)
      integer jev(0:mfr),ibf,nat_o
      real*4 mass(mxat),cs(mxat),wgh(nmax),eps
      real*4 fct(3),pi,ev(0:mfr,nmax),ssr,dd,spcn(3,npt)
      real*4 ai,bi,xm,xr,xi(mip),wi(mip),xc,wc,ss0,alph,bet
      real*4 wint,aa,bb,kbt,q(mqp,3),q0,phi,omg(npt),cc,ss2,ff,ff0 
      real*4 hbar,dwf(mxat,mqp),dotqe,ss1,ss_n0,re,wmx(3),kap(3),debye
      real*4 sqw0(npt),distf(npt),qavs(mxat,npt),sqa(mxat),temp0,rp
      real*4 frqmx,ww,xy,ss(3),www,qq,rrr,ss_n1,sq(3),rrp(mxat),ss_n2
      real*4 enf,akf,aki,alp,frqmx1,sqw1(npt),offst,sev(3),sqw2(npt)
      real*4 dos_t(npt),dos_b(npt),snt
      real*4 evd,ccd,ggd,xx(nmax),yy(nmax)

      complex dot,rr(mxat,3,3),bmat(mxat,3,3),ssc,vct(3),ssa(mxat,3,3)
      complex pom,pp1,pp2
      complex evct(nmax,nmax),vc1,vc2,dotd,p1,p2

      character*1  ch,at(mxat),sg(mxat)
      character*2  ch2
      character*4  arg      
      character*11 name_crd
      
      logical deut,lf_dos

      pi  = 4.0*atan(1.0)
c ####################################################################
c #### DATA given at the command line ################################
      call getarg(1,arg)    ! max k-point suffix n according to files
      read(arg,*)nfr        ! eig_val/vec_n, n = 0..nfr 0: Gamma-point

      call getarg(2,arg)    ! number of q vector orientations
      read(arg,*)nqp        ! for powder averaging

      call getarg(3,arg)    ! n corresponding to k-point defining
      read(arg,*)ixf        ! the linear regime (0-n) of acoustic branch

      lf_dos = .true.       ! logical flag for calculation of DOS
                            ! density-of-states
      frqmx1 = 3000.0       ! freq. range for S(q,w): 0-frqmx1

c experimental data 
      temp0  = 10.0         ! temperature [K]
      deut   = .true.       ! logical flag for deuteration

c instrument geometry (TOSCA)
      enf = 3.0             ! fixed neutron energy of the scattered 
                            ! beam in [meV]
      akf = sqrt(enf/2.09)  
      alp = 135.0 *pi/180.0 ! fixed angle between incident and s
                            ! cattered beam 

c numerical constants
      nip = 9               ! Gaussian integration points <=mip

      name_crd = 'coord'    ! CHARMM coordinate file 
                            ! (for assigning atom types)
c ####################################################################

      call gauss(nip,xi,wi)

c cm-1 units used for frequencies
      kbt  = temp0*0.695028
      hbar   = 33.4664      ! hbar^2/at.m.u [cm-1 A2]

      call assign_atoms(nat,mass,cs,sg,at,deut,name_crd)

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
        read(3,'(28x,3f8.4)')(fct(k),k=1,3)
        write(*,'(a,3f7.3)')'fct:',(fct(k),k=1,3)

        if (fct(1).gt.0.5.or.fct(2).gt.0.5.or.fct(3).gt.0.5) then
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
          if (j.le.10) write(*,'(i4,2f12.6)')j,evd,ggd
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
     1    write(*,'(a,i6,3f12.3)')'rrr',i,rrr,rrp(i),rrr-rrp(i)
      end do

c ---------------------------------------------------------------------
      write(*,*)
      write(*,'(a)')'Powder averaging!'
      write(*,*)
      write(*,'(a,f6.2,a,f6.2)')'Selected temperature T[K] = ',
     1      temp0

      call sphere(mqp,nqp,q)

c B-matrix
      do i=1,nat
        do k1=1,3
          do k2=k1,3
            pom = 0.5*hbar*bmat(i,k1,k2)*rrr/mass(i)
            bmat(i,k1,k2) = pom
            if (k1.ne.k2) bmat(i,k2,k1) = conjg(bmat(i,k1,k2))
          end do
        end do
      end do

c B-matrix output - optional (to check DW factors)
      do j=1,3
        sq(j) = 0.0
      end do
      nat_o = 100  ! no. of atoms
      do i=1,nat_o
        write(*,'(i5,2x,a,4f15.8)')i,at(i),re(bmat(i,1,1)),
     1      re(bmat(i,2,2)),re(bmat(i,3,3)),imag(bmat(i,1,1))
        do j=1,3
          sq(j) = sq(j) + bmat(i,j,j)
        end do
      end do
      write(*,*)
      write(*,'(a,3e16.7)')'sum:',sq(1)/float(nat_o),sq(2)/float(nat_o),
     1                            sq(3)/float(nat_o)

c Debye - Waller factor
      do i=1,nqp
        do j=1,nat
          do k=1,3
            vct(k) = cmplx(0.0,0.0)
            do l=1,3
              vct(k) = vct(k) + bmat(j,k,l)*q(i,l)
            end do
          end do
          ssc = cmplx(0.0,0.0)
          do k=1,3
            ssc = ssc + vct(k)*q(i,k)
          end do
          dwf(j,i) = 0.5*re(ssc)
        end do
      end do

c*******************************************************
c******************  S(q,w) part  **********************
c*******************************************************

c w-range is modified here w:(o..frqmx1)
      
c reset variables
      wint = frqmx1/float(npt)
      do k=1,npt
        omg(k) = wint*(float(k)-0.5)
        sqw0(k)  = 0.0
        sqw1(k)  = 0.0
        sqw2(k)  = 0.0
        dos_t(k) = 0.0
        dos_b(k) = 0.0
      end do

c      write(*,*)'Smoothing factor for frequency peaks:'
c      write(*,*)'bet > 1.0: sharper, bet < 1.0: smoother;    bet = ?'
c      read*,bet
      bet = 2.5 
      eps = 1.e-4
      alph = (4.0/(3.0*wint))**2*bet

      ww = sqrt(log(sqrt(pi/alph)/eps)/alph)
      neigh = nint(ww/wint)
      if (neigh.lt.0) neigh = 1

      do ifr=0,nfr ! loop over phonon k-points

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

        do i=1,nat
          do j=1,npt
            qavs(i,j) = 0.0
          end do
        end do        

c   integration over frequency range -->
c   atomic dynamical structure factor qavs(i,j)

        do ii=1,nqp
        do k=1,npt
          do i=1,nat
            sqa(i) = 0.0
            do k1=1,3
              do k2=1,3
                ssa(i,k1,k2) = cmplx(0.0,0.0)
              end do
            end do
          end do
          if (mod(k,100).eq.0) 
     1    write(*,'(a,i5,a,i5,a,i3,a,i3)')'S(q,w): ',k,'/',npt,'   q:',
     1       ii,'  k:',ifr

          ai = omg(k)-0.5*wint
          bi = omg(k)+0.5*wint
          xm = 0.5*(ai+bi)
          xr = 0.5*(bi-ai)
          do int=1,nip
            xc = xm + xr*xi(int)
            wc = wi(int)*xr
            do j=1,nelm
              if (ev(ifr,j).le.frqmx1) then
                cc = ev(ifr,j)/kbt
                aa = abs(ev(ifr,j)-omg(k))
                bb = (float(neigh)+0.5)*wint
                if (aa.lt.bb) then
                  ff0 = sqrt(alph/pi)*exp(-alph*(ev(ifr,j)-xc)**2)
                  ff  = ff0/(1.0-exp(-cc))/ev(ifr,j)
                  do i=1,nat
                    pom = cmplx(0.0,0.0)
                    do l=1,3
                      j1 = 3*(i-1) + l
                      pom = pom + evct(j1,j)*q(ii,l)
c DOS
                      if (lf_dos.and.ii.eq.1) then
                        do k2=1,3
                          j2 = 3*(i-1) + k2
                          dot = evct(j1,j)*conjg(evct(j2,j))
                          ssa(i,l,k2) = ssa(i,l,k2) +
     1                                  ff0*wc*dot
                        end do
                      end if
                    end do
                    dotqe = re(pom*conjg(pom))
                    sqa(i) = sqa(i) + ff*wc*dotqe
                  end do
                end if
              end if
            end do
          end do

c triangular rule
          aki = sqrt((enf+omg(k)/8.006)/2.09)
          q0 = sqrt(aki**2+akf**2-2.*aki*akf*cos(alp))

          do i=1,nat
            qavs(i,k) = qavs(i,k) + 0.5*cs(i)*hbar*sqa(i)*q0*q0/mass(i)*
     1      exp(-2.0*dwf(i,ii)*q0*q0)/float(nqp)/float(nfr+1)/wint
          end do


c DOS     
          if (lf_dos) then
            do i=1,nat
              www = 0.0
              do k1=1,3
                www = www + re(ssa(i,k1,k1))/float(nfr+1)/wint
              end do
              dos_t(k) = dos_t(k) + www
              if (sg(i).eq.'B') dos_b(k) = dos_b(k) + www
            end do
          end if
        end do
        end do

c sum over atomic contributions
        do i=1,npt
          ss0 = 0.0
          ss1 = 0.0
          ss2 = 0.0

          do j=1,nat
c total system
            ss0 = ss0 + qavs(j,i)
c system without water = biomolecule only
            if (sg(j).eq.'B') ss1 = ss1 + qavs(j,i)
c water only
            if (sg(j).eq.'W') ss2 = ss2 + qavs(j,i)
          end do

c dynamical structure factor
          sqw0(i) = sqw0(i) + ss0/float(nfr+1)*float(npt)
          sqw1(i) = sqw1(i) + ss1/float(nfr+1)*float(npt)
          sqw2(i) = sqw2(i) + ss2/float(nfr+1)*float(npt)
        end do 
      end do      

c*******************************************************
c DOS
      if (lf_dos) then
        open(7,file='DOS.dat', form='formatted',status ='unknown')
        write(7,'(a)')'# omg dos_tot, dos_biomol'
        snt = 0.0
        do k=1,npt
          snt = snt + dos_t(k)
        end do
        do k=1,npt
          write(7,'(f10.3,5f12.7)')omg(k),dos_t(k)/snt,dos_b(k)/snt
        end do
        close(7)
      end if

c SQW
      open(8,file='Sqw_inch',form='formatted',status ='unknown')
      write(8,'(a,f8.3,2i3)')'# temp, n k-pts, n q-pts',
     1                    temp0,(nfr+1),nqp
c normalization for plotting
      offst = 0.0 !25.0
      ss_n0 = 0.0
      ss_n1 = 0.0
      ss_n2 = 0.0
      ii = 0
      do k=1,npt
        if (omg(k).gt.offst) then
          ss_n0 = ss_n0 + sqw0(k)
          ss_n1 = ss_n1 + sqw1(k)
          ss_n2 = ss_n2 + sqw2(k)
          ii    = ii    + 1
        end if
      end do
      ss_n0 = ss_n0/float(ii)
      ss_n1 = ss_n1/float(ii)
      ss_n2 = ss_n2/float(ii)
      write(8,'(a,a)')
     1 '# columns: 1: freq., S(q,w) 2: tot. sys, 3: biomol, 4: water',
     2 ' 5-7: normalized 2-4'
      write(8,'(a,3f10.3)')'# normalization factors:',ss_n0,ss_n1,ss_n2
      do k=1,npt
        write(8,'(f10.3,6f15.5)')omg(k),sqw0(k),sqw1(k),sqw2(k),
     1        sqw0(k)/ss_n0,sqw1(k)/ss_n1,sqw2(k)/ss_n2
      end do
      close(8)

c smothing of the spectrum: resolution = rp% of the energy transfer
      rp = 2.0

      rp = rp/100.
      open(9,file='Sqw_inch_smth',form='formatted',status ='unknown')
      write(9,'(a)')
     1 '# columns: 1: freq., S(q,w) 2: tot. sys, 3: biomol, 4: water'
      wint = omg(npt)/float(npt)
      alph = (2.0/wint)**2
      do k=1,npt
        do j=1,3
          ss(j) = 0.0
        end do
        ssr = 0.0
        alph = (2.0/(wint*(1.0 + rp*omg(k)/omg(1))))**2
        do l=1,npt
          dd = float(k-l)*wint
          ww = sqrt(alph/pi)*exp(-alph*dd**2)
          ssr = ssr + ww
          ss(1) = ss(1) + sqw0(l)*ww
          ss(2) = ss(2) + sqw1(l)*ww
          ss(3) = ss(3) + sqw2(l)*ww
        end do
        do j=1,3
          spcn(j,k) = ss(j)/ssr
        end do
      end do
      do k=1,npt
        write(9,'(f10.3,4f15.5)')omg(k),
     1        spcn(1,k),spcn(2,k),spcn(3,k)
      end do
      close(9)
      end

      subroutine assign_atoms(nat,mass,cs,sg,at,deut,name_crd)
      implicit none
      integer nmax,nat,i,ii,ir,mxat,k
      parameter(mxat=3600)
      parameter(nmax=3*mxat)
      real*4 mass(mxat),cs(mxat),x(mxat),y(mxat),z(mxat),amass,acs

      character*4  aty,seg
      character*6  ars
      character*1  at(mxat),sg(mxat)
      character*11 name_crd

      logical wt,oh,hn,ht,hr,r1,r2,r3,r4,r5,r6,r7,r8,exch,deut

      open(9, file=name_crd, form='formatted', status ='unknown')
c  suppose 4 lines off-set
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

        else if (aty(1:1).eq.'C') then
          if (aty(2:3).eq.'LA') then
            amass = 35.45
c            acs = 5.2  !incoh
            acs = 16.7  !scatt            
          else
            amass = 12.011
c            acs = 0.001  !incoh
            acs = 5.555  !scatt
          end if
        else if (aty(1:1).eq.'N') then
          amass = 14.007
c          acs = 0.49   !incoh
          acs = 11.50  !scatt
        else if (aty(1:1).eq.'O') then
          amass = 15.9994
c          acs = 0.00   !incoh
          acs = 4.24   !scatt
        else if (aty(1:1).eq.'P') then
          amass = 30.974
c          acs = 0.01   !incoh
          acs = 3.31   !scatt
        else if (aty(1:1).eq.'S') then
          amass = 32.07
c          acs = 0.01   !incoh
          acs = 1.03   !scatt
        else if (aty(1:1).eq.'F') then
          amass = 55.85
c          acs = 0.39   !incoh
          acs = 11.83   !scatt
        else 
          write(*,*)'Atom type not defined!   coord - line:',i
          stop
        end if
        mass(i) = amass
        cs(i) = acs
        at(i) = aty(1:1)
        sg(i) = seg(1:1)
c assuming the initial of the segment name of water being 'W'
        if (sg(i).ne.'W') sg(i)='B' ! select biomolecule
      end do
      close(9)
      return
      end

      subroutine gauss(n,x,w)
c  weights and coordinates for Gausian quadrature of level n
c  from: Flannery et al., Numerical Recepies
      implicit none
      integer mip,i,j,m,n
      parameter (mip=10)
      real*4 eps,p1,p2,p3,pp,x1,x2,xl,xm,z,z1,x,w
      parameter(eps=3.e-6)
      dimension x(mip),w(mip)
      x1=-1.0
      x2= 1.0
      m=(n+1)/2
      xm=0.5*(x2+x1)
      xl=0.5*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654*(i-0.25)/(n+0.5))
1       continue
          p1=1.0
          p2=0.0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.0)
          z1=z
          z=z1-p1/pp
          write(*,*)abs(z-z1)
        if(abs(z-z1).gt.eps)go to 1
        x(i)=xm-xl*z
        x(n-i+1)=xm+xl*z
        w(i)=2.0*xl/((1.0-z*z)*pp*pp)
        w(n-i+1)=w(i)
12    continue
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
      parameter (mip=10)
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

      subroutine sphere(mx,nx,qq)
      implicit none
      integer i,mx,nx,idum,j,k
      real*4 w(3),rnd,x,q,th,fi,rr,pi,phi,dd0,dd,qq(mx,3)
      logical flg

      idum = 1357
      rr = 1.0
      pi = 4.0*atan(1.0)

      phi = 2.0/sqrt(float(nx))
      dd0 = rr*phi
      do i=1,nx
 100    continue
        x = rnd(idum)
        q = 1.0-2.0*x
        th = atan(sqrt(1.0-q*q)/q)
        if (th.lt.0.0) th = pi+th
        x = rnd(idum)
        fi = 2.0*pi*x
        
        w(1) = rr*sin(th)*cos(fi)
        w(2) = rr*sin(th)*sin(fi)
        w(3) = rr*cos(th)
        flg = .true.
        do j=1,i-1
          dd = (qq(j,1)-w(1))**2+(qq(j,2)-w(2))**2+(qq(j,3)-w(3))**2
          dd = sqrt(dd)
          if (dd.lt.dd0) flg = .false.
        end do
        if (flg.or.i.eq.1) then
          qq(i,1) = w(1)
          qq(i,2) = w(2)
          qq(i,3) = w(3)
        else
          goto 100
        end if
      end do
      return
      end

      REAL*4 FUNCTION RND(ISEED)
      implicit none
      INTEGER ISEED
      REAL*8 DSEED,DIVIS,DENOM,MULTIP,RNDP
C
      DATA  DIVIS/2147483647.D0/
      DATA  DENOM /2147483711.D0/
      DATA  MULTIP/16807.D0/
C
      IF(ISEED.LE.1) ISEED=314159
      DSEED=MULTIP*ISEED
      DSEED=MOD(DSEED,DIVIS)
      RNDP=DSEED/DENOM
      ISEED=DSEED
      RND = RNDP
C
      RETURN
      END

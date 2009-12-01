      subroutine read_orb_loc
c Written by Cyrus Umrigar
c Reads in either analytic or numerical orbitals for finite systems (molecules, Je spheres, quantum dots)
      use all_tools_mod
      use mpi_mod
      use const_mod
      use contr2_mod
      use contr3_mod
      use contr_ylm_mod
      use dim_mod
      use basis1_mod, only : zex
      use basis2_mod
      use numbas_mod
      use coefs_mod
      use wfsec_mod
      use distance_mod
      use orbital_num_mod
      implicit real*8(a-h,o-z)

      dimension orb(nelec,norb),dorb(3,nelec,norb),ddorb(nelec,norb),r(2)

c Check if orbitals_num exists.  If it does not exist, we will create it from analytic orbs.
c Do not confuse analytic orbs with analytic basis fns.
      if(inum_orb.ne.0) then
        num_orb_exist=0
        open(4,file='orbitals_num',form='formatted',status='old',err=10)
        num_orb_exist=1
      endif

   10 if(inum_orb.ne.0 .and. num_orb_exist.eq.1) then ! numerical orbs
        call my_second(1,'read_o')
        call read_orb_loc_num
        write(6,'(''read in numerical orbitals'')')
        close(4)
        call my_second(2,'read_o')
       else ! analytical orbitals (basis functions can be either analytical or numerical)
        if(numr.ge.-3 .and. numr .le. 1) then
          call read_orb_loc_ana
         else
          stop 'numr must be between -3 and 1'
        endif
c       if(maxval(zex).gt.0.d0) call distinct_radial_bas
        call distinct_radial_bas
        write(6,'(''Done reading in LCAO coefs. of orbitals'')')
      endif

c Check that irecursion_ylm=1 if l of basis function >=5
      do 20 ibasis=1,nbasis
        if(l_bas(ibasis).ge.5 .and. irecursion_ylm.eq.0) then
          write(6,'(''if basis functions with l>=5 are used, set irecursion_ylm=1 in read_input'')')
          stop 'if basis functions with l>=5 are used, set irecursion_ylm=1 in read_input'
        endif
   20 continue

c The rest of this routine is for numerical orbs and so we could have one big if statement
c Generate grid and orbitals on grid if using numerical orbitals that do not already exist on disk
      if(inum_orb.ne.0 .and. num_orb_exist.eq.0) then
        read(5,*) ngrid_orbx,ngrid_orby,sizex,sizey
        write(6,'(''ngrid_orbx,ngrid_orby,sizex,sizey='',2i5,2f8.3)') ngrid_orbx,ngrid_orby,sizex,sizey
!JT        if(ngrid_orbx.gt.MGRID_ORB) stop 'ngrid_orbx > MGRID_ORB in read_orb_loc'
!JT        if(ngrid_orby.gt.MGRID_ORB) stop 'ngrid_orby > MGRID_ORB in read_orb_loc'
      endif

c Set grid info.  Grid extends from -sizex to sizex, etc.
      if(inum_orb.ne.0) then
        call alloc ('xorb_grid',  xorb_grid, ngrid_orbx)
        call alloc ('yorb_grid',  yorb_grid, ngrid_orby)
        hx=2*sizex/(ngrid_orbx-1)
        hy=2*sizey/(ngrid_orby-1)
        hxi=1/hx
        hyi=1/hy
        do 30 ix=1,ngrid_orbx
   30     xorb_grid(ix)=-sizex+(ix-1)*hx
        do 35 iy=1,ngrid_orby
   35     yorb_grid(iy)=-sizey+(iy-1)*hy
      endif

c Calculate orbitals on grid and spline them
      if(inum_orb.ne.0 .and. num_orb_exist.eq.0) then
  
        call alloc ('orb_num',  orb_num, 4, ngrid_orbx, ngrid_orby, norb)

        do 76 ix=1,ngrid_orbx
          r(1)=-sizex+(ix-1)*hx
          do 76 iy=1,ngrid_orby
            r(2)=-sizey+(iy-1)*hy

c Calculate e-N inter-particle distances
            call distancese(1,r)

            iwf=1
            call orbitals_loc_ana(1,rvec_en,r_en,orb,dorb,ddorb)
c           call orbitals_loc_ana_grade(1,rvec_en,r_en,orb,dorb,ddorb)

c Warning: Note that we need to save only the orbitals at all grid pts. and
c we need to save the 1st derivatives only on the boundary.
c However, we are presently saving everything everywhere.
            do 76 iorb=1,norb
              orb_num(1,ix,iy,iorb)=orb(1,iorb)
              orb_num(4,ix,iy,iorb)=ddorb(1,iorb)
              do 76 k=1,ndim
   76           orb_num(k+1,ix,iy,iorb)=dorb(k,1,iorb)

        call my_second(2,'orb_nu')

c Write orbitals_num if it does not exist and inum_orb.ne.0
        if(index(mode,'mpi').eq.0 .or. idtask.eq.0) then
c         open(4,file='orbitals_num',form='unformatted',status='new')
          open(4,file='orbitals_num',form='formatted',status='new')
          write(4,'(2i5,2f8.4,i4)') ngrid_orbx,ngrid_orby,sizex,sizey,norb
          do 90 iorb=1,norb
   90       write(4,'(1p10d22.14)') ((orb_num(1,ix,iy,iorb),ix=1,ngrid_orbx),iy=1,ngrid_orby)
          close(4)
        endif

        call my_second(2,'writeo')

c       write(6,'(''orbnum='',/,(5f12.8))') ((orb_num(1,ix,iy,1),iy=0,ngrid_orby-1),ix=0,ngrid_orbx-1)

      endif

c Construct the spline
      if(inum_orb.ne.0) call spline_orb(num_orb_exist)

      return
      end
c-----------------------------------------------------------------------

      subroutine read_orb_loc_ana
c Written by Cyrus Umrigar
c Reads in localized orbitals, in one or more of
c 1) a slater basis
c 2) a gaussian basis
c 3) a gauss-slater basis
c In addition there can be asymptotic function (nsa,npa,nda) at the end
c For the p and d basis functions, I maintain the old order to avoid
c having to change old inputs.
c The order I read in the p functions is x,y,z, which is m=1,-1,0 and
c the order that I read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, which is m=0,2,-2,1,-1,
c in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.
      use all_tools_mod
      use const_mod
      use dim_mod
      use atom_mod
      use basis1_mod
      use basis2_mod
      use numbas_mod
      use coefs_mod
      use pseudo_mod
      implicit real*8(a-h,o-z)

c Note: the order that we read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.

      allocate(n1s(nctype))
      allocate(n2s(nctype),n2p(-1:1,nctype))
      allocate(n3s(nctype),n3p(-1:1,nctype),n3d(-2:2,nctype))
      allocate(n4s(nctype),n4p(-1:1,nctype),n4d(-2:2,nctype),n4f(-3:3,nctype))
      allocate(n5s(nctype),n5p(-1:1,nctype),n5d(-2:2,nctype),n5f(-3:3,nctype),n5g(-4:4,nctype))
      allocate(n6d(-2:2,nctype),n6f(-3:3,nctype),n6g(-4:4,nctype),n6h(-5:5,nctype))
c     allocate(n7g(-4:4,nctype),n7h(-5:5,nctype),n7i(-6:6,nctype))
c     allocate(n8i(-6:6,nctype),n8j(-7:7,nctype))
c     allocate(n9k(-8:8,nctype))
c     allocate(n10l(-9:9,nctype))
c     allocate(n11m(-10:10,nctype))
c     allocate(n12n(-11:11,nctype))
c     allocate(n13o(-12:12,nctype))
      allocate(nsa(nctype),npa(-1:1,nctype),nda(-2:2,nctype))
      call alloc ('nbasis_ctype', nbasis_ctype, nctype)

c     initialize arrays
      n1s=0
      n2s=0; n2p=0
      n3s=0; n3p=0; n3d=0
      n4s=0; n4p=0; n4d=0; n4f=0
      n5s=0; n5p=0; n5d=0; n5f=0; n5g=0
      n6d=0; n6f=0; n6g=0; n6h=0
c     n7g=0; n7h=0; n7i=0;
c     n8i=0; n8j=0
c     n9k=0
c     n10l=0
c     n11m=0
c     n12n=0
c     n13o=0
      nsa=0; npa=0; nda=0

      nbas_tot=0
      do 50 ict=1,nctype

        if(ndim.eq.2) then
          if(numr.le.0) then
            read(5,*) n1s(ict),n2p(1,ict),n2p(-1,ict),n3d(2,ict),n3d(-2,ict)
     &      ,n4f(3,ict),n4f(-3,ict),n5g(4,ict),n5g(-4,ict),n6h(5,ict),n6h(-5,ict)
           else ! numerical radial basis functions
            read(5,*) (m_bas(ib),ib=1,nbasis)
            write(6,'(/,''m for basis functions:'',/,100i5)') (m_bas(ib),ib=1,nbasis)
          endif
         elseif(ndim.eq.3) then
          if(numr.le.0) then
            read(5,*) n1s(ict)
     &      ,n2s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &      ,n3s(ict),n3p(1,ict),n3p(-1,ict),n3p(0,ict)
     &      ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &      ,n4s(ict),n4p(1,ict),n4p(-1,ict),n4p(0,ict)
     &      ,n4d(0,ict),(n4d(m,ict),n4d(-m,ict),m=2,1,-1)
     &      ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &      ,n5s(ict),n5p(1,ict),n5p(-1,ict),n5p(0,ict)
     &      ,n5d(0,ict),(n5d(m,ict),n5d(-m,ict),m=2,1,-1)
     &      ,(n5f(m,ict),n5f(-m,ict),m=3,1,-1),n5f(0,ict)
     &      ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &      ,nsa(ict),npa(1,ict),npa(-1,ict),npa(0,ict)
     &      ,nda(0,ict),(nda(m,ict),nda(-m,ict),m=2,1,-1)
            write(6,'(''n1s etc:'',90i2)') n1s(ict)
     &      ,n2s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &      ,n3s(ict),n3p(1,ict),n3p(-1,ict),n3p(0,ict)
     &      ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &      ,n4s(ict),n4p(1,ict),n4p(-1,ict),n4p(0,ict)
     &      ,n4d(0,ict),(n4d(m,ict),n4d(-m,ict),m=2,1,-1)
     &      ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &      ,n5s(ict),n5p(1,ict),n5p(-1,ict),n5p(0,ict)
     &      ,n5d(0,ict),(n5d(m,ict),n5d(-m,ict),m=2,1,-1)
     &      ,(n5f(m,ict),n5f(-m,ict),m=3,1,-1),n5f(0,ict)
     &      ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &      ,nsa(ict),npa(1,ict),npa(-1,ict),npa(0,ict)
     &      ,nda(0,ict),(nda(m,ict),nda(-m,ict),m=2,1,-1)
           else ! numerical radial basis functions
            read(5,*) n1s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &      ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &      ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &      ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &      ,(n6h(m,ict),n6h(-m,ict),m=5,1,-1),n6h(0,ict)
            write(6,'(''ns etc:'',90i2)') n1s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &      ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &      ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &      ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &      ,(n6h(m,ict),n6h(-m,ict),m=5,1,-1),n6h(0,ict)
          endif
        endif

c Check that the total number of basis functions for each center type is correct
        nbas_typ=           iabs(n1s(ict))+iabs(n2s(ict))+iabs(n3s(ict))+iabs(n4s(ict))+iabs(n5s(ict))
        nbas_tot=nbas_tot + iabs(n1s(ict))+iabs(n2s(ict))+iabs(n3s(ict))+iabs(n4s(ict))+iabs(n5s(ict))
        do 21 m=-1,1
          nbas_typ=nbas_typ + iabs(n2p(m,ict)) + iabs(n3p(m,ict)) + iabs(n4p(m,ict)) + iabs(n5p(m,ict))
   21     nbas_tot=nbas_tot + iabs(n2p(m,ict)) + iabs(n3p(m,ict)) + iabs(n4p(m,ict)) + iabs(n5p(m,ict))
        do 22 m=-2,2
          nbas_typ=nbas_typ + iabs(n3d(m,ict)) + iabs(n4d(m,ict)) + iabs(n5d(m,ict))
   22     nbas_tot=nbas_tot + iabs(n3d(m,ict)) + iabs(n4d(m,ict)) + iabs(n5d(m,ict))
        do 23 m=-3,3
          nbas_typ=nbas_typ + iabs(n4f(m,ict)) + iabs(n5f(m,ict))
   23     nbas_tot=nbas_tot + iabs(n4f(m,ict)) + iabs(n5f(m,ict))
        do 24 m=-4,4
          nbas_typ=nbas_typ + iabs(n5g(m,ict))
   24     nbas_tot=nbas_tot + iabs(n5g(m,ict))
        do 25 m=-5,5
          nbas_typ=nbas_typ + iabs(n6h(m,ict))
   25     nbas_tot=nbas_tot + iabs(n6h(m,ict))

        nbasis_ctype(ict)=nbas_typ

c Read in which radial basis function is used by each basis function
c We are now reading this in even if the radial basis functions are analytical to allow for a mixed analyltical-numerical basis
        mbasis_ctype = maxval(nbasis_ctype)
        if(ndim.eq.3 .or. (ndim.eq.2.and.numr.eq.1)) then
          call alloc ('iwrwf', iwrwf, mbasis_ctype ,nctype)
          write(6,'(''Reading iwrwf'')')
          read(5,*) (iwrwf(ibct,ict),ibct=1,nbasis_ctype(ict))
          write(6,'(''Center'',i5,'' uses radial bas. fns:'',(100i3))') ict,(iwrwf(ibct,ict),ibct=1,nbasis_ctype(ict))
        endif
        nbas_tot=nbas_tot+(ncent_ctype(ict)-1)*nbasis_ctype(ict)
 
   50 continue

      if(nbas_tot.ne.nbasis) then
        write(6,'(''nbas_tot,nbasis='',9i4)') nbas_tot,nbasis
        stop 'nbas_tot not equal to nbasis'
      endif

c Calculate betaq for asymptotic basis functions
      betaq=0
      do 55 ic=1,ncent
   55   betaq=betaq+znuc(iwctype(ic))
      betaq=betaq-nelec+1

      write(6,'(/,''Basis functions:'')')
      write(6,'(''center type'',t12,(12i3))') (i,i=1,nctype)

      if(maxval(abs(n1s)).ne.0) write(6,'(''1s'',t12,(12i3))') (n1s(i),i=1,nctype)

      if(maxval(abs(n2s)).ne.0) write(6,'(''2s'',t12,(12i3))') (n2s(i),i=1,nctype)
      if(maxval(abs(n2p(1,:))).ne.0)  write(6,'(''2px'',t12,(12i3))') (n2p(1,i),i=1,nctype)
      if(maxval(abs(n2p(-1,:))).ne.0) write(6,'(''2py'',t12,(12i3))') (n2p(-1,i),i=1,nctype)
      if(maxval(abs(n2p(0,:))).ne.0) write(6,'(''2pz'',t12,(12i3))') (n2p(0,i),i=1,nctype)

      if(maxval(abs(n3s)).ne.0) write(6,'(''3s'',t12,(12i3))') (n3s(i),i=1,nctype)
      if(maxval(abs(n3p(1,:))).ne.0)  write(6,'(''3px'',t12,(12i3))') (n3p(1,i),i=1,nctype)
      if(maxval(abs(n3p(-1,:))).ne.0) write(6,'(''3py'',t12,(12i3))') (n3p(-1,i),i=1,nctype)
      if(maxval(abs(n3p(0,:))).ne.0) write(6,'(''3pz'',t12,(12i3))') (n3p(0,i),i=1,nctype)
      if(maxval(abs(n3d(0,:))).ne.0) write(6,'(''3dzr'',t12,(12i3))') (n3d(0,i),i=1,nctype)
      if(maxval(abs(n3d(2,:))).ne.0) write(6,'(''3dx2'',t12,(12i3))') (n3d(2,i),i=1,nctype)
      if(maxval(abs(n3d(-2,:))).ne.0) write(6,'(''3dxy'',t12,(12i3))') (n3d(-2,i),i=1,nctype)
      if(maxval(abs(n3d(1,:))).ne.0) write(6,'(''3dxz'',t12,(12i3))') (n3d(1,i),i=1,nctype)
      if(maxval(abs(n3d(-1,:))).ne.0) write(6,'(''3dyz'',t12,(12i3))') (n3d(-1,i),i=1,nctype)

      if(maxval(abs(n4s)).ne.0) write(6,'(''4s'',t12,(12i3))') (n4s(i),i=1,nctype)
      if(maxval(abs(n4p(1,:))).ne.0)  write(6,'(''4px'',t12,(12i3))') (n4p(1,i),i=1,nctype)
      if(maxval(abs(n4p(-1,:))).ne.0) write(6,'(''4py'',t12,(12i3))') (n4p(-1,i),i=1,nctype)
      if(maxval(abs(n4p(0,:))).ne.0) write(6,'(''4pz'',t12,(12i3))') (n4p(0,i),i=1,nctype)
      if(maxval(abs(n4d(0,:))).ne.0) write(6,'(''4dzr'',t12,(12i3))') (n4d(0,i),i=1,nctype)
      if(maxval(abs(n4d(2,:))).ne.0) write(6,'(''4dx2'',t12,(12i3))') (n4d(2,i),i=1,nctype)
      if(maxval(abs(n4d(-2,:))).ne.0) write(6,'(''4dxy'',t12,(12i3))') (n4d(-2,i),i=1,nctype)
      if(maxval(abs(n4d(1,:))).ne.0) write(6,'(''4dxz'',t12,(12i3))') (n4d(1,i),i=1,nctype)
      if(maxval(abs(n4d(-1,:))).ne.0) write(6,'(''4dyz'',t12,(12i3))') (n4d(-1,i),i=1,nctype)
      do 60 m=-3,3
   60   if(maxval(abs(n4f(m,:))).ne.0) write(6,'(''4f('',i2,'')'',t12,(12i3))') m,(n4f(m,i),i=1,nctype)

      if(maxval(abs(n5s)).ne.0) write(6,'(''4s'',t12,(12i3))') (n5s(i),i=1,nctype)
      if(maxval(abs(n5p(1,:))).ne.0)  write(6,'(''5px'',t12,(12i3))') (n5p(1,i),i=1,nctype)
      if(maxval(abs(n5p(-1,:))).ne.0) write(6,'(''5py'',t12,(12i3))') (n5p(-1,i),i=1,nctype)
      if(maxval(abs(n5p(0,:))).ne.0) write(6,'(''5pz'',t12,(12i3))') (n5p(0,i),i=1,nctype)
      if(maxval(abs(n5d(0,:))).ne.0) write(6,'(''5dzr'',t12,(12i3))') (n5d(0,i),i=1,nctype)
      if(maxval(abs(n5d(2,:))).ne.0) write(6,'(''5dx2'',t12,(12i3))') (n5d(2,i),i=1,nctype)
      if(maxval(abs(n5d(-2,:))).ne.0) write(6,'(''5dxy'',t12,(12i3))') (n5d(-2,i),i=1,nctype)
      if(maxval(abs(n5d(1,:))).ne.0) write(6,'(''5dxz'',t12,(12i3))') (n5d(1,i),i=1,nctype)
      if(maxval(abs(n5d(-1,:))).ne.0) write(6,'(''5dyz'',t12,(12i3))') (n5d(-1,i),i=1,nctype)
      do 70 m=-3,3
   70   if(maxval(abs(n5f(m,:))).ne.0) write(6,'(''5f('',i2,'')'',t12,(12i3))') m,(n5f(m,i),i=1,nctype)
      do 80 m=-4,4
   80   if(maxval(abs(n5g(m,:))).ne.0) write(6,'(''5g('',i2,'')'',t12,(12i3))') m,(n5g(m,i),i=1,nctype)

      if(maxval(abs(nsa)).ne.0) write(6,'(''sa'',t12,(12i3))') (nsa(i),i=1,nctype)
      if(maxval(abs(npa(1,:))).ne.0) write(6,'(''pxa'',t12,(12i3))') (npa(1,i),i=1,nctype)
      if(maxval(abs(npa(-1,:))).ne.0) write(6,'(''pya'',t12,(12i3))') (npa(-1,i),i=1,nctype)
      if(maxval(abs(npa(0,:))).ne.0) write(6,'(''pza'',t12,(12i3))') (npa(0,i),i=1,nctype)
      if(maxval(abs(nda(0,:))).ne.0) write(6,'(''dzra'',t12,(12i3))') (nda(0,i),i=1,nctype)
      if(maxval(abs(nda(2,:))).ne.0) write(6,'(''dx2a'',t12,(12i3))') (nda(2,i),i=1,nctype)
      if(maxval(abs(nda(-2,:))).ne.0) write(6,'(''dxya'',t12,(12i3))') (nda(-2,i),i=1,nctype)
      if(maxval(abs(nda(1,:))).ne.0) write(6,'(''dxza'',t12,(12i3))') (nda(1,i),i=1,nctype)
      if(maxval(abs(nda(-1,:))).ne.0) write(6,'(''dyza'',t12,(12i3))') (nda(-1,i),i=1,nctype)

      write(6,'(/,''charge'',t12,(12i3))')(nint(znuc(i)),i=1,nctype)
      write(6,*)

c Asign character labels, n_bas, l_bas, m_bas, centers and center types to each of the basis functions
c The basis fns. in the LCAO coefs. are either in the atomic filling order
c or in order of increasing l.  The former is the order I used to use the
c latter is the order from GAMESS when doing all-electron calculations.
c Either order will work when numr.eq.1.
      if(numr.eq.0 .or. numr.eq.1) then
        call orb_loc_ana_original_order
       elseif(numr.eq.-1 .or. numr.eq.-2 .or. numr.eq.-3) then
        call orb_loc_ana_gamess_order
       else
        stop 'numr must be between -3 and 1'
      endif

      call object_modified ('nbasis_ctype')
      call object_modified ('mbasis_ctype')

      write(6,'(''orbital coefficients'')')
      do 260 i=1,norb
        read(5,*) (coef(j,i,1),j=1,nbasis)
  260   write(6,'(12f10.6)') (coef(j,i,1),j=1,nbasis)

      write(6,'(''screening constants'')')
      read(5,*) (zex(i,1),i=1,nbasis)
      write(6,'(12f10.6)') (zex(i,1),i=1,nbasis)
      write(6,*)
      if(minval(zex(:,1)).lt.0.d0) then
        write(6,'(''Basis exponents must be >= 0.  If 0 then numerical radial basis fn is used'')')
        stop 'Basis exponents must be >= 0.  If 0 then numerical radial basis fn is used'
      endif

      do 262 i=1,nbasis
        if(zex(i,1).lt.0.d0) then
          write(6,'(''read_orb_loc: basis exponents zex must be >=0'')')
          stop 'read_orb_loc: basis exponents zex must be >=0'
        endif
  262 continue

      return
      end
c-----------------------------------------------------------------------

      subroutine orb_loc_ana_original_order
c Written by Cyrus Umrigar
c Reads in localized orbitals, in one or more of
c 1) a slater basis
c 2) a gaussian basis
c 3) a gauss-slater basis
c In addition there can be asymptotic function (nsa,npa,nda) at the end
c For the p and d basis functions, I maintain the old order to avoid
c having to change old inputs.
c The order I read in the p functions is x,y,z, which is m=1,-1,0 and
c the order that I read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, which is m=0,2,-2,1,-1,
c in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.
      use all_tools_mod
      use atom_mod
      use coefs_mod
      use basis1_mod
      use basis2_mod
      use lbas_mod
      implicit real*8(a-h,o-z)

      parameter(nprime=5)

      character*6 l1s(nprime),l2s(nprime),l2p(nprime,-1:1)
     &,l3s(nprime),l3p(nprime,-1:1),l3d(nprime,-2:2)
     &,l4s(nprime),l4p(nprime,-1:1),l4d(nprime,-2:2)
     &,lsa,lpa(-1:1),lda(-2:2)
      character*2 lmbas
      character*4 lcent,l4f(nprime),l5g(nprime),l6h(nprime)

      data l1s /'    1s','   1s''','   1s"','  1s"''','  1s""'/
      data l2s /'    2s','   2s''','   2s"','  2s"''','  2s""'/
      data l3s /'    3s','   3s''','   3s"','  3s"''','  3s""'/
      data l4s /'    4s','   4s''','   4s"','  4s"''','  4s""'/
      data l2p /'   2py','  2py''','  2py"',' 2py"''',' 2py""'
     &         ,'   2pz','  2pz''','  2pz"',' 2pz"''',' 2pz""'
     &         ,'   2px','  2px''','  2px"',' 2px"''',' 2px""'/
      data l3p /'   3py','  3py''','  3py"',' 3py"''',' 3py""'
     &         ,'   3pz','  3pz''','  3pz"',' 3pz"''',' 3pz""'
     &         ,'   3px','  3px''','  3px"',' 3px"''',' 3px""'/
      data l4p /'   4py','  4py''','  4py"',' 4py"''',' 4py""'
     &         ,'   4pz','  4pz''','  4pz"',' 4pz"''',' 4pz""'
     &         ,'   4px','  4px''','  4px"',' 4px"''',' 4px""'/
      data l3d /'  3dxy',' 3dxy''',' 3dxy"','3dxy"''','3dxy""'
     &         ,'  3dyz',' 3dyz''',' 3dyz"','3dyz"''','3dyz""'
     &         ,'  3dzr',' 3dzr''',' 3dzr"','3dzr"''','3dzr""'
     &         ,'  3dxz',' 3dxz''',' 3dxz"','3dxz"''','3dxz""'
     &         ,'  3dx2',' 3dx2''',' 3dx2"','3dx2"''','3dx2""'/
      data l4d /'  4dxy',' 4dxy''',' 4dxy"','4dxy"''','4dxy""'
     &         ,'  4dyz',' 4dyz''',' 4dyz"','4dyz"''','4dyz""'
     &         ,'  4dzr',' 4dzr''',' 4dzr"','4dzr"''','4dzr""'
     &         ,'  4dxz',' 4dxz''',' 4dxz"','4dxz"''','4dxz""'
     &         ,'  4dx2',' 4dx2''',' 4dx2"','4dx2"''','4dx2""'/
      data l4f /' f','f''','f"','f"''','f""'/
      data l5g /' g','g''','g"','g"''','g""'/
      data l6h /' h','h''','h"','h"''','h""'/
      data lsa /'    sa'/
      data lpa /'   pya','   pza','   pxa'/
      data lda /'   dxy','   dyz','   dzr','   dxz','   dx2'/

c the following sets up character strings that identify basis
c functions on centers for the printout of coefficients and
c screening constants. Character strings are of form
c 2px'(3) for the second 2px function on the third center

      call alloc ('n_bas', n_bas, nbasis)
      call alloc ('l_bas', l_bas, nbasis)
      call alloc ('m_bas', m_bas, nbasis)
      call alloc ('icenter_basis', icenter_basis, nbasis)
      call alloc ('ictype_basis', ictype_basis, nbasis)
      call alloc ('lbasis', lbasis, nbasis)

      ib=0
      do 250 ic=1,ncent
        i=iwctype(ic)

        write(lcent,'(''('',i2,'')'')') ic
        if(ic.lt.10)  write(lcent,'(''(0'',i1,'')'')') ic

        do 90 j=1,iabs(n1s(i))
c         if(iabs(n1s(i)).gt.nprime) stop 'number of 1s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(1,n1s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
   90     lbasis(ib)=l1s(min(j,nprime))//lcent

        do 100 j=1,iabs(n2s(i))
c         if(iabs(n2s(i)).gt.nprime) stop 'number of 2s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(2,n2s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  100     lbasis(ib)=l2s(min(j,nprime))//lcent

c Here we are preseving the input order, i.e., m=1,-1,0 (see loop over 108)
        do 110 m=1,0,-1
          do 105 j=1,iabs(n2p(m,i))
c           if(iabs(n2p(m,i)).gt.nprime) stop 'number of 2p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(2,n2p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  105       lbasis(ib)=l2p(min(j,nprime),m)//lcent
          if(m.eq.0) goto 110
          do 108 j=1,iabs(n2p(-m,i))
c           if(iabs(n2p(-m,i)).gt.nprime) stop 'number of 2p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(2,n2p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  108       lbasis(ib)=l2p(min(j,nprime),-m)//lcent
  110   continue

        do 120 j=1,iabs(n3s(i))
c         if(iabs(n3s(i)).gt.nprime) stop 'number of 3s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(3,n3s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  120     lbasis(ib)=l3s(min(j,nprime))//lcent

        do 130 m=1,0,-1
          do 125 j=1,iabs(n3p(m,i))
c           if(iabs(n3p(m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  125       lbasis(ib)=l3p(min(j,nprime),m)//lcent
          if(m.eq.0) goto 130
          do 128 j=1,iabs(n3p(-m,i))
c           if(iabs(n3p(-m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  128       lbasis(ib)=l3p(min(j,nprime),-m)//lcent
  130   continue

c Here we are preseving the input order, i.e., m=0,2,-2,1,-1 (see loop over 140)
          m=0
          do 132 j=1,iabs(n3d(m,i))
c           if(iabs(n3d(m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  132       lbasis(ib)=l3d(min(j,nprime),m)//lcent
        do 140 m=2,1,-1
          do 135 j=1,iabs(n3d(m,i))
c           if(iabs(n3d(m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  135       lbasis(ib)=l3d(min(j,nprime),m)//lcent
          do 140 j=1,iabs(n3d(-m,i))
c           if(iabs(n3d(-m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(-m,i))
            l_bas(ib)=2
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  140       lbasis(ib)=l3d(min(j,nprime),-m)//lcent

        do 150 j=1,iabs(n4s(i))
c         if(iabs(n4s(i)).gt.nprime) stop 'number of 4s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(4,n4s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
          ib=ib+1
  150     lbasis(ib)=l4s(min(j,nprime))//lcent

        do 160 m=1,0,-1
          do 155 j=1,iabs(n4p(m,i))
c           if(iabs(n4p(m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  155       lbasis(ib)=l4p(min(j,nprime),m)//lcent
          if(m.eq.0) goto 160
          do 158 j=1,iabs(n4p(-m,i))
c           if(iabs(n4p(-m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  158       lbasis(ib)=l4p(min(j,nprime),-m)//lcent
  160   continue

          m=0
          do 162 j=1,iabs(n4d(m,i))
            if(iabs(n4d(m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  162       lbasis(ib)=l4d(min(j,nprime),m)//lcent
        do 170 m=2,1,-1
          do 165 j=1,iabs(n4d(m,i))
            if(iabs(n4d(m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  165       lbasis(ib)=l4d(min(j,nprime),m)//lcent
          do 170 j=1,iabs(n4d(-m,i))
            if(iabs(n4d(-m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4d(-m,i))
            l_bas(ib)=2
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  170       lbasis(ib)=l4d(min(j,nprime),-m)//lcent

c Ordering is now consistent with input ordering
        do 185 m=3,0,-1
          do 180 j=1,iabs(n4f(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 4f basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4f(m,i))
            l_bas(ib)=3
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
c           write(lmbas,'(''('',i2,'')'')') m
            write(lmbas,'(i2)') m
  180       lbasis(ib)=l4f(min(j,nprime))//lmbas//lcent
          if(m.ne.0) then
          do 182 j=1,iabs(n4f(-m,i))
c           if(iabs(n4f(-m,i)).gt.nprime) stop 'number of 4f basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4f(-m,i))
            l_bas(ib)=3
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
c           write(lmbas,'(''('',i2,'')'')') -m
            write(lmbas,'(i2)') -m
  182       lbasis(ib)=l4f(min(j,nprime))//lmbas//lcent
          endif
  185   continue

        do 195 m=4,0,-1
          do 190 j=1,iabs(n5g(m,i))
c           if(iabs(n5g(m,i)).gt.nprime) stop 'number of 5g basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5g(m,i))
            l_bas(ib)=4
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') m
  190       lbasis(ib)=l5g(min(j,nprime))//lmbas//lcent
          if(m.ne.0) then
          do 192 j=1,iabs(n5g(-m,i))
c           if(iabs(n5g(-m,i)).gt.nprime) stop 'number of 5g basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5g(-m,i))
            l_bas(ib)=4
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') -m
  192       lbasis(ib)=l5g(min(j,nprime))//lmbas//lcent
          endif
  195   continue

        do 205 m=5,0,-1
          do 200 j=1,iabs(n6h(m,i))
c           if(iabs(n6h(m,i)).gt.nprime) stop 'number of 6h basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(6,n6h(m,i))
            l_bas(ib)=5
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') m
  200       lbasis(ib)=l6h(min(j,nprime))//lmbas//lcent
          if(m.ne.0) then
          do 202 j=1,iabs(n6h(-m,i))
c           if(iabs(n6h(-m,i)).gt.nprime) stop 'number of 6h basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(6,n6h(-m,i))
            l_bas(ib)=5
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') -m
  202       lbasis(ib)=l6h(min(j,nprime))//lmbas//lcent
          endif
  205   continue

        if(iabs(nsa(i)).ge.1) then
          ib=ib+1
          n_bas(ib)=0
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
          lbasis(ib)=lsa//lcent
        endif

        do 220 m=1,-1,-1
          if(iabs(npa(m,i)).ge.1) then
            ib=ib+1
            n_bas(ib)=0
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            lbasis(ib)=lpa(m)//lcent
          endif
  220   continue

c Here we are preseving the input order, i.e., m=0,2,-2,1,-1 (see loop over 232)
          m=0
          do 225 j=1,iabs(nda(m,i))
c           if(iabs(nda(m,i)).gt.nprime) stop 'number of da basis fns > nprime'
            ib=ib+1
            n_bas(ib)=0
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  225       lbasis(ib)=lda(m)//lcent
        do 235 m=2,1,-1
          do 230 j=1,iabs(nda(m,i))
c           if(iabs(nda(m,i)).gt.nprime) stop 'number of da basis fns > nprime'
            ib=ib+1
            n_bas(ib)=0
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  230       lbasis(ib)=lda(m)//lcent
          do 232 j=1,iabs(nda(-m,i))
c           if(iabs(nda(-m,i)).gt.nprime) stop 'number of da basis fns > nprime'
            ib=ib+1
            n_bas(ib)=0
            l_bas(ib)=2
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  232       lbasis(ib)=lda(-m)//lcent
  235     continue

  250 continue

!     for now, ML_BAS needs to be at least 4 because of the code in basis_fns.f
      ML_BAS = max(4, maxval(l_bas))

      write(6,'(/,(12a10))') (lbasis(j),j=1,nbasis)

      return
      end
c-----------------------------------------------------------------------

      subroutine orb_loc_ana_gamess_order
c Written by Cyrus Umrigar
c Reads in localized orbitals, in one or more of
c 1) a slater basis
c 2) a gaussian basis
c 3) a gauss-slater basis
c In addition there can be asymptotic function (nsa,npa,nda) at the end
c For the p and d basis functions, I maintain the old order to avoid
c having to change old inputs.
c The order I read in the p functions is x,y,z, which is m=1,-1,0 and
c the order that I read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, which is m=0,2,-2,1,-1,
c in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.
c This routine differs from read_orb_loc_ana in that the order of the fns.
c is changed from
c 1s, 2s, 2px, 2py, ...          to
c 1s, 2s, 3s, 4s,  2px, 2py, 2pz, 3px, ... ., 4pz,  3d, 4f, 5g, 6h
      use all_tools_mod
      use atom_mod
      use coefs_mod
      use basis1_mod
      use basis2_mod
      use lbas_mod
      implicit real*8(a-h,o-z)

      parameter(nprime=5)

      character*6 l1s(nprime),l2s(nprime),l2p(nprime,-1:1)
     &,l3s(nprime),l3p(nprime,-1:1),l3d(nprime,-2:2)
     &,l4s(nprime),l4p(nprime,-1:1),l4d(nprime,-2:2)
     &,l5s(nprime),l5p(nprime,-1:1),l5d(nprime,-2:2)
     &,lsa,lpa(-1:1),lda(-2:2)
      character*2 lmbas
      character*4 lcent,l4f(nprime),l5f(nprime),l5g(nprime),l6h(nprime)

      data l1s /'    1s','   1s''','   1s"','  1s"''','  1s""'/
      data l2s /'    2s','   2s''','   2s"','  2s"''','  2s""'/
      data l3s /'    3s','   3s''','   3s"','  3s"''','  3s""'/
      data l4s /'    4s','   4s''','   4s"','  4s"''','  4s""'/
      data l5s /'    5s','   5s''','   5s"','  5s"''','  5s""'/
      data l2p /'   2py','  2py''','  2py"',' 2py"''',' 2py""'
     &         ,'   2pz','  2pz''','  2pz"',' 2pz"''',' 2pz""'
     &         ,'   2px','  2px''','  2px"',' 2px"''',' 2px""'/
      data l3p /'   3py','  3py''','  3py"',' 3py"''',' 3py""'
     &         ,'   3pz','  3pz''','  3pz"',' 3pz"''',' 3pz""'
     &         ,'   3px','  3px''','  3px"',' 3px"''',' 3px""'/
      data l4p /'   4py','  4py''','  4py"',' 4py"''',' 4py""'
     &         ,'   4pz','  4pz''','  4pz"',' 4pz"''',' 4pz""'
     &         ,'   4px','  4px''','  4px"',' 4px"''',' 4px""'/
      data l5p /'   5py','  5py''','  5py"',' 5py"''',' 5py""'
     &         ,'   5pz','  5pz''','  5pz"',' 5pz"''',' 5pz""'
     &         ,'   5px','  5px''','  5px"',' 5px"''',' 5px""'/
      data l3d /'  3dxy',' 3dxy''',' 3dxy"','3dxy"''','3dxy""'
     &         ,'  3dyz',' 3dyz''',' 3dyz"','3dyz"''','3dyz""'
     &         ,'  3dzr',' 3dzr''',' 3dzr"','3dzr"''','3dzr""'
     &         ,'  3dxz',' 3dxz''',' 3dxz"','3dxz"''','3dxz""'
     &         ,'  3dx2',' 3dx2''',' 3dx2"','3dx2"''','3dx2""'/
      data l4d /'  4dxy',' 4dxy''',' 4dxy"','4dxy"''','4dxy""'
     &         ,'  4dyz',' 4dyz''',' 4dyz"','4dyz"''','4dyz""'
     &         ,'  4dzr',' 4dzr''',' 4dzr"','4dzr"''','4dzr""'
     &         ,'  4dxz',' 4dxz''',' 4dxz"','4dxz"''','4dxz""'
     &         ,'  4dx2',' 4dx2''',' 4dx2"','4dx2"''','4dx2""'/
      data l5d /'  5dxy',' 5dxy''',' 5dxy"','5dxy"''','5dxy""'
     &         ,'  5dyz',' 5dyz''',' 5dyz"','5dyz"''','5dyz""'
     &         ,'  5dzr',' 5dzr''',' 5dzr"','5dzr"''','5dzr""'
     &         ,'  5dxz',' 5dxz''',' 5dxz"','5dxz"''','5dxz""'
     &         ,'  5dx2',' 5dx2''',' 5dx2"','5dx2"''','5dx2""'/
      data l4f /' f','f''','f"','f"''','f""'/
      data l5f /' f','f''','f"','f"''','f""'/
      data l5g /' g','g''','g"','g"''','g""'/
      data l6h /' h','h''','h"','h"''','h""'/
      data lsa /'    sa'/
      data lpa /'   pya','   pza','   pxa'/
      data lda /'   dxy','   dyz','   dzr','   dxz','   dx2'/

c the following sets up character strings that identify basis
c functions on centers for the printout of coefficients and
c screening constants. Character strings are of form
c 2px'(3) for the second 2px function on the third center

      call alloc ('n_bas', n_bas, nbasis)
      call alloc ('l_bas', l_bas, nbasis)
      call alloc ('m_bas', m_bas, nbasis)
      call alloc ('icenter_basis', icenter_basis, nbasis)
      call alloc ('ictype_basis', ictype_basis, nbasis)
      call alloc ('lbasis', lbasis, nbasis)

      ib=0
      do 250 ic=1,ncent
        i=iwctype(ic)

        write(lcent,'(''('',i2,'')'')') ic
        if(ic.lt.10)  write(lcent,'(''(0'',i1,'')'')') ic

        do 90 j=1,iabs(n1s(i))
c         if(iabs(n1s(i)).gt.nprime) stop 'number of 1s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(1,n1s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
   90     lbasis(ib)=l1s(min(j,nprime))//lcent

        do 100 j=1,iabs(n2s(i))
c         if(iabs(n2s(i)).gt.nprime) stop 'number of 2s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(2,n2s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  100     lbasis(ib)=l2s(min(j,nprime))//lcent

        do 120 j=1,iabs(n3s(i))
c         if(iabs(n3s(i)).gt.nprime) stop 'number of 3s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(3,n3s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  120     lbasis(ib)=l3s(min(j,nprime))//lcent

        do 150 j=1,iabs(n4s(i))
c         if(iabs(n4s(i)).gt.nprime) stop 'number of 4s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(4,n4s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  150     lbasis(ib)=l4s(min(j,nprime))//lcent

        do 151 j=1,iabs(n5s(i))
c         if(iabs(n5s(i)).gt.nprime) stop 'number of 5s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(5,n5s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  151     lbasis(ib)=l5s(min(j,nprime))//lcent

c Here we are preseving the input order, i.e., m=1,-1,0 (see loop over 108)
        do 110 m=1,0,-1
          do 105 j=1,iabs(n2p(m,i))
c           if(iabs(n2p(m,i)).gt.nprime) stop 'number of 2p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(2,n2p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  105       lbasis(ib)=l2p(min(j,nprime),m)//lcent
          do 125 j=1,iabs(n3p(m,i))
c           if(iabs(n3p(m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  125       lbasis(ib)=l3p(min(j,nprime),m)//lcent
          do 155 j=1,iabs(n4p(m,i))
c           if(iabs(n4p(m,i)).gt.nprime) stop 'number of 4p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  155       lbasis(ib)=l4p(min(j,nprime),m)//lcent
          do 156 j=1,iabs(n5p(m,i))
c           if(iabs(n5p(m,i)).gt.nprime) stop 'number of 5p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  156       lbasis(ib)=l5p(min(j,nprime),m)//lcent
          if(m.eq.0) goto 110
          do 108 j=1,iabs(n2p(-m,i))
c           if(iabs(n2p(-m,i)).gt.nprime) stop 'number of 2p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(2,n2p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  108       lbasis(ib)=l2p(min(j,nprime),-m)//lcent
          do 128 j=1,iabs(n3p(-m,i))
c           if(iabs(n3p(-m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  128       lbasis(ib)=l3p(min(j,nprime),-m)//lcent
          do 158 j=1,iabs(n4p(-m,i))
c           if(iabs(n4p(-m,i)).gt.nprime) stop 'number of 4p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  158       lbasis(ib)=l4p(min(j,nprime),-m)//lcent
          do 159 j=1,iabs(n5p(-m,i))
c           if(iabs(n5p(-m,i)).gt.nprime) stop 'number of 5p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  159       lbasis(ib)=l5p(min(j,nprime),-m)//lcent
  110   continue

c Here we are preseving the input order, i.e., m=0,2,-2,1,-1 (see loop over 140)
          m=0
          do 132 j=1,iabs(n3d(m,i))
c           if(iabs(n3d(m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  132       lbasis(ib)=l3d(min(j,nprime),m)//lcent
          do 162 j=1,iabs(n4d(m,i))
c           if(iabs(n4d(m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  162       lbasis(ib)=l4d(min(j,nprime),m)//lcent
          do 163 j=1,iabs(n5d(m,i))
c           if(iabs(n5d(m,i)).gt.nprime) stop 'number of 5d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  163       lbasis(ib)=l5d(min(j,nprime),m)//lcent
        do 175 m=2,1,-1
          do 135 j=1,iabs(n3d(m,i))
c           if(iabs(n3d(m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  135       lbasis(ib)=l3d(min(j,nprime),m)//lcent
          do 165 j=1,iabs(n4d(m,i))
c           if(iabs(n4d(m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  165       lbasis(ib)=l4d(min(j,nprime),m)//lcent
          do 166 j=1,iabs(n5d(m,i))
c           if(iabs(n5d(m,i)).gt.nprime) stop 'number of 5d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  166       lbasis(ib)=l5d(min(j,nprime),m)//lcent
          do 140 j=1,iabs(n3d(-m,i))
c           if(iabs(n3d(-m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(-m,i))
            l_bas(ib)=2
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  140       lbasis(ib)=l3d(min(j,nprime),-m)//lcent
         do 170 j=1,iabs(n4d(-m,i))
c           if(iabs(n4d(-m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
           ib=ib+1
           n_bas(ib)=sign(4,n4d(-m,i))
           l_bas(ib)=2
           m_bas(ib)=-m
           icenter_basis(ib)=ic
           ictype_basis(ib)=i
 170       lbasis(ib)=l4d(min(j,nprime),-m)//lcent
         do 171 j=1,iabs(n5d(-m,i))
c           if(iabs(n5d(-m,i)).gt.nprime) stop 'number of 5d basis fns > nprime'
           ib=ib+1
           n_bas(ib)=sign(5,n5d(-m,i))
           l_bas(ib)=2
           m_bas(ib)=-m
           icenter_basis(ib)=ic
           ictype_basis(ib)=i
 171       lbasis(ib)=l5d(min(j,nprime),-m)//lcent
  175     continue

c Ordering is now consistent with input ordering
        do 185 m=3,0,-1
          do 180 j=1,iabs(n4f(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 4f basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4f(m,i))
            l_bas(ib)=3
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
c           write(lmbas,'(''('',i2,'')'')') m
            write(lmbas,'(i2)') m
  180       lbasis(ib)=l4f(min(j,nprime))//lmbas//lcent
          do 181 j=1,iabs(n5f(m,i))
c           if(iabs(n5f(m,i)).gt.nprime) stop 'number of 5f basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5f(m,i))
            l_bas(ib)=3
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
c           write(lmbas,'(''('',i2,'')'')') m
            write(lmbas,'(i2)') m
  181       lbasis(ib)=l5f(min(j,nprime))//lmbas//lcent
          if (m.ne.0) then ! JT
          do 182 j=1,iabs(n4f(-m,i))
c           if(iabs(n4f(-m,i)).gt.nprime) stop 'number of 4f basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4f(-m,i))
            l_bas(ib)=3
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
c           write(lmbas,'(''('',i2,'')'')') -m
            write(lmbas,'(i2)') -m
  182       lbasis(ib)=l4f(min(j,nprime))//lmbas//lcent
          do 183 j=1,iabs(n5f(-m,i))
c           if(iabs(n5f(-m,i)).gt.nprime) stop 'number of 5f basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5f(-m,i))
            l_bas(ib)=3
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
c           write(lmbas,'(''('',i2,'')'')') -m
            write(lmbas,'(i2)') -m
  183       lbasis(ib)=l5f(min(j,nprime))//lmbas//lcent
          endif ! JT
  185   continue

        do 195 m=4,0,-1
          do 190 j=1,iabs(n5g(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 5g basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5g(m,i))
            l_bas(ib)=4
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') m
  190       lbasis(ib)=l5g(min(j,nprime))//lmbas//lcent
          if(m.ne.0) then
          do 192 j=1,iabs(n5g(-m,i))
c           if(iabs(n4f(-m,i)).gt.nprime) stop 'number of 5g basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5g(-m,i))
            l_bas(ib)=4
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') -m
  192       lbasis(ib)=l5g(min(j,nprime))//lmbas//lcent
          endif
  195   continue

        do 205 m=5,0,-1
          do 200 j=1,iabs(n6h(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 6h basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(6,n6h(m,i))
            l_bas(ib)=5
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') m
  200       lbasis(ib)=l6h(min(j,nprime))//lmbas//lcent
          if(m.ne.0) then
          do 202 j=1,iabs(n6h(-m,i))
c           if(iabs(n4f(-m,i)).gt.nprime) stop 'number of 6h basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(6,n6h(-m,i))
            l_bas(ib)=5
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') -m
  202       lbasis(ib)=l6h(min(j,nprime))//lmbas//lcent
          endif
  205   continue

        if(iabs(nsa(i)).ge.1) then
          ib=ib+1
          n_bas(ib)=0
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
          lbasis(ib)=lsa//lcent
        endif

        do 220 m=1,-1,-1
          if(iabs(npa(m,i)).ge.1) then
            ib=ib+1
            n_bas(ib)=0
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            lbasis(ib)=lpa(m)//lcent
          endif
  220   continue

c Here we are preseving the input order, i.e., m=0,2,-2,1,-1 (see loop over 232)
          m=0
          do 225 j=1,iabs(nda(m,i))
c           if(iabs(nda(m,i)).gt.nprime) stop 'number of da basis fns > nprime'
            ib=ib+1
            n_bas(ib)=0
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  225       lbasis(ib)=lda(m)//lcent
        do 235 m=2,1,-1
          do 230 j=1,iabs(nda(m,i))
c           if(iabs(nda(m,i)).gt.nprime) stop 'number of da basis fns > nprime'
            ib=ib+1
            n_bas(ib)=0
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  230       lbasis(ib)=lda(m)//lcent
          do 232 j=1,iabs(nda(-m,i))
c           if(iabs(nda(-m,i)).gt.nprime) stop 'number of da basis fns > nprime'
            ib=ib+1
            n_bas(ib)=0
            l_bas(ib)=2
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  232       lbasis(ib)=lda(-m)//lcent
  235     continue

  250 continue

      ML_BAS = max(4, maxval(l_bas))

      write(6,'(/,(12a10))') (lbasis(j),j=1,nbasis)

      return
      end
c-----------------------------------------------------------------------

      subroutine distinct_radial_bas
c Written by Cyrus Umrigar
c Determine distinct radial basis functions
      use all_tools_mod
      use const_mod, only : ipr
      use atom_mod
      use coefs_mod
      use basis1_mod
      use numbas_mod
      use basis2_mod
      use wfsec_mod
      use contr3_mod
      use forcepar_mod
      implicit real*8(a-h,o-z)

c If an analytical basis is used find out how many distinct radial orbitals there are
c and store n and zex for each of these.
c This is done so that only the unique analytic radial wavefns. are evaluated in basis_fns.f and basis_fnse.f
c zex2 is used in func when optimizing the exponents.
c ib  goes up to the total number of basis functions
c ib2 goes up to the number of basis functions on each center type
c irb goes up to the number of radial basis functions on each center type
c Warning: If this is an optimization run (either fit or new optimization) that optimizes exponents
c then it is possible that exponents that start out being equal may not end that way.
c In particular this happens because GAMESS creates an s function out
c of d functions and these start off with the same exponent.
c To be correct and efficient, one could check for equality of just n_bas and zex when exponents are not being optimized
c and check in addition for equality of l_bas when exponents are being optimized.
c The problem is that one does not know this until later in the input.  So, at present we always check for equality of all 3 things.
c Now that we are allowing a mixed analytical and numerical basis, we use this routine to find the distinct analytical
c basis functions and put them first and the numerical ones will be put later.
c The decision of which are analytical and which are numerical is based on whether zex is nonzero or zero.

!JT      write(6,'(''beg of distinct_rad_bas'')')
!JT      do 10 ict=1,nctype
!JT        write(6,'(''ict='',20i5)') ict
!JT   10   write(6,'(''distinct: iwrwf='',20i5)') (iwrwf(ibct,ict),ibct=1,nbasis_ctype(ict))
 
!JT      write(6,'(''beg of distinct_radial_bas'')')
      call alloc ('iwrwf2', iwrwf2, nbasis)
      call alloc ('nrbas_analytical', nrbas_analytical, nctype)
!     allocate iwrwf here in the case when it is not allocated and read in before (new style format input)
      call alloc ('iwrwf', iwrwf, mbasis_ctype ,nctype)

c Figure out the number of different analytical radial basis functions
c Note that we are checking not only for the equality of of n_bas and zex but also of l_bas.  The latter is not necessary if this is
c not an optimization run, but since we do not know at this stage if it is or not, we do the safe thing.
c Note also that in principle different centers of the same center type could have different basis exponents for the analytical basis functions
c (although the numerical ones must be the same).  However, the dimensioning of iwrwf and zex2 does not allow this.
c So, although we are looping over all centers below, we could loop over just the center types and copy the information to centers of the same center type.
      if(ipr.ge.2) write(6,'(''n_bas(ib),l_bas(ib)'',100(2i2,x))') (n_bas(ib),l_bas(ib),ib=1,nbasis)
      ib=0
      kbct=0
      do 290 ic=1,ncent
        ict=iwctype(ic)
        nrbas_analytical(ict)=0
        do 280 ibct=1,nbasis_ctype(ict)
          ib=ib+1
          if(ib.ne.ibct+kbct) stop 'ib .ne. ibct+kbct'
          if(zex(ib,1).gt.0.d0) then
            do 270 jbct=1,ibct-1  ! check if any of the earlier basis functions on the same atom uses the same analytical radial basis fn.
              if(zex(jbct+kbct,1).gt.0.d0) then
                if(ipr.ge.2) then
                  write(6,'(''ibct,jbct,l_bas(ib),l_bas(jbct+kbct),n_bas(ib),n_bas(jbct+kbct),zex(ib,1),zex(jbct+kbct,1)'',
     &            6i5,9f15.8)')
     &            ibct,jbct,l_bas(ib),l_bas(jbct+kbct),n_bas(ib),n_bas(jbct+kbct),zex(ib,1),zex(jbct+kbct,1)
                endif
c               if((index(mode,'fit').eq.0 .and. n_bas(ib).eq.n_bas(jbct+kbct) .and. zex(ib,1).eq.zex(jbct+kbct,1)) .or.
                if((l_bas(ib).eq.l_bas(jbct+kbct) .and. n_bas(ib).eq.n_bas(jbct+kbct) .and. zex(ib,1).eq.zex(jbct+kbct,1))) then
                  irb=iwrwf2(jbct+kbct)
                  goto 275
                endif
              endif
  270       continue
            nrbas_analytical(ict)=nrbas_analytical(ict)+1
            irb=nrbas_analytical(ict)
!JT         if(irb.gt.MRWF) stop 'nbas > MRWF'
!           MRWF is initilized to 0 in numbas_mod
!           it cannot be initialized in this routine because
!           this routine is called several times for the correlated sampling part
!           of optimization runs, and this would destroy zex2
            MRWF = max (MRWF, nrbas_analytical(ict))
            call alloc ('n_bas2', n_bas2, MRWF, nctype)
            call alloc ('zex2', zex2, MRWF, nctype, nwf)
            n_bas2(irb,ict)=n_bas(ib)
            zex2(irb,ict,iwf)=zex(ib,iwf) ! JT: 1 -> iwf
  275       iwrwf2(ib)=irb
            iwrwf(ibct,ict)=irb
          endif
  280     continue
  290   kbct=kbct+nbasis_ctype(ict)

      if(ipr.ge.2) then
        write(6,'(''iwrwf2='',40i3)') (iwrwf2(ib),ib=1,nbasis)
        do 300 ict=1,nctype
          write(6,'(''ict,nrbas_analytical(ict)='',20i5)') ict,nrbas_analytical(ict)
          write(6,'(''iwrwf='',20i5)') (iwrwf(ibct,ict),ibct=1,nbasis_ctype(ict))
          write(6,'(''nbas2='',20i5)') (n_bas2(irb,ict),irb=1,nrbas_analytical(ict))
  300     write(6,'(''zex2='',20f6.2)') (zex2(irb,ict,iwf),irb=1,nrbas_analytical(ict))
      endif

      call object_modified ('zex2')

      return
      end
c-----------------------------------------------------------------------

      subroutine copy_zex_zex2
c Written by Cyrus Umrigar, 21 Oct 09
c Copy zex to zex2
c zex and zex2 are the same, but indexed differently.
      use atom_mod,   only : nctype, ncent_ctype
      use basis1_mod, only : zex
      use basis2_mod, only : zex2, nbasis_ctype
      use numbas_mod, only : iwrwf, nrbas_analytical
      use wfsec_mod,  only : iwf
      implicit real*8(a-h,o-z)

      ib=0
      do 20 ict=1,nctype
        do 10 ibct=1,nbasis_ctype(ict)
          ib=ib+1
          irb=iwrwf(ibct,ict)
          if(irb.le.nrbas_analytical(ict)) then
            write(6,'(''zex2('',2i3,'')=zex('',i3,'')'')') irb,ict,ib
            zex2(irb,ict,iwf)=zex(ib,iwf)
          endif
   10   continue
   20   ib=ib+(ncent_ctype(ict)-1)*nbasis_ctype(ict)

      return
      end
c-----------------------------------------------------------------------

      subroutine read_orb_loc_num
c Reads in numerical orbitals and V_ext on a grid.
c At present V_ext that is read in is not used.
      use all_tools_mod
      use coefs_mod
      use orbital_num_mod
      implicit real*8(a-h,o-z)

c     dimension v_ext(MGRID_ORB,MGRID_ORB)

      open(4,file='orbitals_num',form='formatted',status='old')

      read(4,*) ngrid_orbx,ngrid_orby,sizex,sizey,norba
!JT      if(ngrid_orbx.gt.MGRID_ORB .or. ngrid_orby.gt.MGRID_ORB) then
!JT          write(6,'(a,i10,a,i10,a)') 'ngrid_orbx=',ngrid_orbx,' ngrid_orby=',ngrid_orby,' > MGRID_ORB'
!JT          stop 'ngrid_orbx,ngrid_orby for numerical orbitals must be < MGRID_ORB'
!JT      endif
      if(mod(ngrid_orbx,2).ne.1) stop 'ngrid_orbx must be odd in read_orb_loc_num'
      if(mod(ngrid_orby,2).ne.1) stop 'ngrid_orby must be odd in read_orb_loc_num'
      if(norba.lt.norb) stop 'norba not large enough in read_orb_loc_num'

c     read(4,*) ((v_ext(i,j),i=1,ngrid_orbx),j=1,ngrid_orby)
      call alloc ('orb_num', orb_num, 4, ngrid_orbx, ngrid_orby, norba)
      do 10 iorb=1,norba
   10   read(4,*) ((orb_num(1,i,j,iorb),i=1,ngrid_orbx),j=1,ngrid_orby)

c     write(6,'(''orb'',i2,(65d10.2))') (((orb_num(1,i,j,iorb),i=1,ngrid_orbx),j=1,ngrid_orby),iorb=1,2)

      close(4)

      return
      end
c-----------------------------------------------------------------------

      subroutine spline_orb(num_orb_exist)
c Spline numerical orbitals
      use coefs_mod
      use orbital_num_mod
      implicit real*8(a-h,o-z)
      parameter(MWORK=21)

      dimension bcxmin(ngrid_orby),bcxmax(ngrid_orby),bcymin(ngrid_orbx),bcymax(ngrid_orbx),
     &wk(MWORK*ngrid_orbx*ngrid_orby)

      nwk=MWORK*ngrid_orbx*ngrid_orby

c ibc* sets the boundary condition.  If we are reading in orbitals on a grid
c then nothing is known about derivatives so use "not a knot" b.c.
c If we are calculating orbitals on a grid calculate 1st derivatives and use them.
      if(num_orb_exist.eq.1) then
        ibcxmin=0
        ibcxmax=0
        ibcymin=0
        ibcymax=0
       else
        ibcxmin=1
        ibcxmax=1
        ibcymin=1
        ibcymax=1
      endif

c     do 30 ix=1,ngrid_orbx
c  30   xorb_grid(ix)=-sizex+(ix-1)*hx
c     do 35 iy=1,ngrid_orby
c  35   yorb_grid(iy)=-sizey+(iy-1)*hy

c Takes data in orb_num(1,ix,iy,iorb) as input and calculates bicubic spline
c coeff's in orb_num(2-4,ix,iy,iorb)
      do 40 iorb=1,norb
        if(num_orb_exist.eq.1) then
          do 10 iy=1,ngrid_orby
            bcxmin(iy)=orb_num(2,1,iy,iorb)
   10       bcxmax(iy)=orb_num(2,ngrid_orbx,iy,iorb)
          do 20 ix=1,ngrid_orbx
            bcymin(ix)=orb_num(3,ix,1,iorb)
   20       bcymax(ix)=orb_num(3,ix,ngrid_orbx,iorb)
        endif
        call r8mkbicubw(xorb_grid,ngrid_orbx,yorb_grid,ngrid_orby,orb_num(1,1,1,iorb),ngrid_orbx,
     &    ibcxmin,bcxmin,ibcxmax,bcxmax,
     &    ibcymin,bcymin,ibcymax,bcymax,
     &    wk,nwk,ilinx,iliny,ier)
   40   if(ier.ne.0) stop 'error in r8mkbicubw'

c     write(6,'(''xorb_grid(1),'',6i9,20d12.4)') ngrid_orbx,ngrid_orby,nwk,ilinx,iliny,ier,xorb_grid(1),yorb_grid(2),((orb_num(i,2,3,iorb),i=1,4),iorb=1,2)

c ict(i) is a control array, used by splint_orb.
c ict:  1 = compute, 0 = do not compute.
c i=1 f
c   2 df/dx
c   3 df/dy
c   4 d2f/dx2
c   5 d2f/dy2
c   6 d2f/dxdy
      do 50 n=1,5
   50   ict(n)=1
      ict(6)=0

      write(6,'(''orb_num1b'',10d12.4)') (orb_num(1,1,1,iorb),iorb=1,norb)

      return
      end
c-----------------------------------------------------------------------

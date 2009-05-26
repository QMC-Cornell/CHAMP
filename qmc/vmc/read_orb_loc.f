      subroutine read_orb_loc
c Written by Cyrus Umrigar
c Reads in either analytic or localized orbitals

      use mpi_mod
      use orbitals_mod, only: orb_tot_nb

      use coefs_mod
      use const_mod
      use dim_mod
      use numbas_mod
      use basis2_mod
      use contr2_mod
      use wfsec_mod
      use contr3_mod
      use distance_mod
      use contr_ylm_mod
      implicit real*8(a-h,o-z)


      common /orbital_num/ orb_num(4,MGRID_ORB,MGRID_ORB,MORB_OCC),xorb_grid(MGRID_ORB),yorb_grid(MGRID_ORB)
     &,sizex,sizey,hx,hy,hxi,hyi,ngrid_orbx,ngrid_orby,ict(6)

      dimension orb(nelec,orb_tot_nb),dorb(3,nelec,orb_tot_nb),ddorb(nelec,orb_tot_nb),r(2)

c Check if orbitals_num exists.  If it does not exist, we will create it from analytic orbs.
c Do not confuse analytic orbs with analytic basis fns.
      if(inum_orb.ne.0) then
        num_orb_exist=0
        open(4,file='orbitals_num',form='formatted',status='old',err=10)
        num_orb_exist=1
      endif

   10 if(inum_orb.ne.0 .and. num_orb_exist.eq.1) then
        call my_second(1,'read_o')
        call read_orb_loc_num
        write(6,'(''read in numerical orbitals'')')
        close(4)
        call my_second(2,'read_o')
       else
        if(numr.eq.0 .or. numr.eq.1 .or. numr.eq.-1 .or. numr.eq.-2 .or. numr.eq.-3 .or. numr.eq.-4) then
          call read_orb_loc_ana
         else
          stop 'numr must be between -4 and 1'
        endif
        call distinct_radial_bas
        write(6,'(''read in analytical orbitals'')')
      endif

c Check that irecursion_ylm=1 if l of basis function >=4
      do 20 ibasis=1,nbasis
        if(l_bas(ibasis).gt.ML_BAS) then
          write(6,'(''(l_bas(ibasis) > ML_BAS'')')
          stop 'l_bas(ibasis) > ML_BAS'
        endif
        if(l_bas(ibasis).ge.5 .and. irecursion_ylm.eq.0) then
          write(6,'(''if basis functions with l>=5 are used, set irecursion_ylm=1 in read_input'')')
          stop 'if basis functions with l>=5 are used, set irecursion_ylm=1 in read_input'
        endif
   20 continue

      if(inum_orb.ne.0 .and. num_orb_exist.eq.0) then
        read(5,*) ngrid_orbx,ngrid_orby,sizex,sizey
        write(6,'(''ngrid_orbx,ngrid_orby,sizex,sizey='',2i5,2f8.3)') ngrid_orbx,ngrid_orby,sizex,sizey
        if(ngrid_orbx.gt.MGRID_ORB) stop 'ngrid_orbx > MGRID_ORB in read_orb_loc'
        if(ngrid_orby.gt.MGRID_ORB) stop 'ngrid_orby > MGRID_ORB in read_orb_loc'
      endif

c Set grid info.  Grid extends from -sizex to sizex, etc.
      if(inum_orb.ne.0) then
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
c 2) an asymptotic basis
c 3) a gaussian basis
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
      use const_mod
      use dim_mod
      use numbas_mod
      use basis2_mod
      use pseudo_mod
      implicit real*8(a-h,o-z)



c     character*6 l1s(nprime),l2s(nprime),l2p(nprime,-1:1)
c    &,l3s(nprime),l3p(nprime,-1:1),l3d(nprime,-2:2)
cc   &,l4s(nprime),l4p(nprime,-1:1),l4d(nprime,-2:2)
c    &,l4s(nprime),l4p(nprime,-1:1)
c    &,lsa,lpa(-1:1),lda(-2:2)
c     character*2 lmbas
c     character*4 lcent,l4f(nprime),l5g(nprime),l6h(nprime)
c     character*10 lbasis(MBASIS)

c     data l1s /'    1s','   1s''','   1s"','  1s"''','  1s""'/
c     data l2s /'    2s','   2s''','   2s"','  2s"''','  2s""'/
c     data l3s /'    3s','   3s''','   3s"','  3s"''','  3s""'/
c     data l4s /'    4s','   4s''','   4s"','  4s"''','  4s""'/
c     data l2p /'   2py','  2py''','  2py"',' 2py"''',' 2py""'
c    &         ,'   2pz','  2pz''','  2pz"',' 2pz"''',' 2pz""'
c    &         ,'   2px','  2px''','  2px"',' 2px"''',' 2px""'/
c     data l3p /'   3py','  3py''','  3py"',' 3py"''',' 3py""'
c    &         ,'   3pz','  3pz''','  3pz"',' 3pz"''',' 3pz""'
c    &         ,'   3px','  3px''','  3px"',' 3px"''',' 3px""'/
c     data l4p /'   4py','  4py''','  4py"',' 4py"''',' 4py""'
c    &         ,'   4pz','  4pz''','  4pz"',' 4pz"''',' 4pz""'
c    &         ,'   4px','  4px''','  4px"',' 4px"''',' 4px""'/
c     data l3d /'  3dxy',' 3dxy''',' 3dxy"','3dxy"''','3dxy""'
c    &         ,'  3dyz',' 3dyz''',' 3dyz"','3dyz"''','3dyz""'
c    &         ,'  3dzr',' 3dzr''',' 3dzr"','3dzr"''','3dzr""'
c    &         ,'  3dxz',' 3dxz''',' 3dxz"','3dxz"''','3dxz""'
c    &         ,'  3dx2',' 3dx2''',' 3dx2"','3dx2"''','3dx2""'/
c     data l4d /'  4dxy',' 4dxy''',' 4dxy"','4dxy"''','4dxy""'
c    &         ,'  4dyz',' 4dyz''',' 4dyz"','4dyz"''','4dyz""'
c    &         ,'  4dzr',' 4dzr''',' 4dzr"','4dzr"''','4dzr""'
c    &         ,'  4dxz',' 4dxz''',' 4dxz"','4dxz"''','4dxz""'
c    &         ,'  4dx2',' 4dx2''',' 4dx2"','4dx2"''','4dx2""'/
c     data l4f /' f','f''','f"','f"''','f""'/
c     data l5g /' g','g''','g"','g"''','g""'/
c     data l6h /' h','h''','h"','h"''','h""'/
c     data lsa /'    sa'/
c     data lpa /'   pya','   pza','   pxa'/
c     data lda /'   dxy','   dyz','   dzr','   dxz','   dx2'/

c Note: the order that we read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.

!     allocations
      allocate(n1s(nctype))
      allocate(n2s(nctype),n2p(-1:1,nctype))
      allocate(n3s(nctype),n3p(-1:1,nctype),n3d(-2:2,nctype))
      allocate(n4s(nctype),n4p(-1:1,nctype),n4d(-2:2,nctype),n4f(-3:3,nctype))
      allocate(n5s(nctype),n5p(-1:1,nctype),n5d(-2:2,nctype),n5f(-3:3,nctype),n5g(-4:4,nctype))
      allocate(n6d(-2:2,nctype),n6f(-3:3,nctype),n6g(-4:4,nctype),n6h(-5:5,nctype))
      allocate(n7g(-4:4,nctype),n7h(-5:5,nctype),n7i(-6:6,nctype))
      allocate(n8i(-6:6,nctype),n8j(-7:7,nctype))
      allocate(n9k(-8:8,nctype))
      allocate(n10l(-9:9,nctype))
      allocate(n11m(-10:10,nctype))
      allocate(n12n(-11:11,nctype))
      allocate(n13o(-12:12,nctype))
      allocate(nsa(nctype),npa(-1:1,nctype),nda(-2:2,nctype))

!     initializations
      n1s=0
      n2s=0; n2p=0
      n3s=0; n3p=0; n3d=0
      n4s=0; n4p=0; n4d=0; n4f=0
      n5s=0; n5p=0; n5d=0; n5f=0; n5g=0
      n6d=0; n6f=0; n6g=0; n6h=0
      n7g=0; n7h=0; n7i=0;
      n8i=0; n8j=0
      n9k=0
      n10l=0
      n11m=0
      n12n=0
      n13o=0
      nsa=0; npa=0; nda=0

      if(numr.le.0) then
c Analytic radial basis
        do 10 ict=1,nctype
          if(ndim.eq.2) then
            read(5,*) n1s(ict),n2p(1,ict),n2p(-1,ict),n3d(2,ict),n3d(-2,ict)
     &      ,n4f(3,ict),n4f(-3,ict),n5g(4,ict),n5g(-4,ict),n6h(5,ict),n6h(-5,ict)
            n2s(ict)=0
            n2p(0,ict)=0
            n3d(0,ict)=0
            n2s(ict)=0
            n3s(ict)=0
            n4s(ict)=0
            n2p(0,ict)=0
            n3d(0,ict)=0
            do 2 m=-1,1
              n3d(m,ict)=0
    2         n4d(m,ict)=0
            do 3 m=-2,2
    3         n4f(m,ict)=0
            do 4 m=-3,3
    4         n5g(m,ict)=0
            do 5 m=-4,4
    5         n6h(m,ict)=0
           else
! JT: if numr = -2, allows up to 4f functions
           if(numr.eq.-2) then
            read(5,*) n1s(ict)
     &      ,n2s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &      ,n3s(ict),n3p(1,ict),n3p(-1,ict),n3p(0,ict)
     &      ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &      ,n4s(ict),n4p(1,ict),n4p(-1,ict),n4p(0,ict)
     &      ,n4d(0,ict),(n4d(m,ict),n4d(-m,ict),m=2,1,-1)
     &      ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &      ,nsa(ict),npa(1,ict),npa(-1,ict),npa(0,ict)
     &      ,nda(0,ict),(nda(m,ict),nda(-m,ict),m=2,1,-1)
! JT: if numr = -3, allows up to 5g functions
           elseif(numr.eq.-3 .or. numr.eq.-4) then
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
           else
            read(5,*) n1s(ict)
     &      ,n2s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &      ,n3s(ict),n3p(1,ict),n3p(-1,ict),n3p(0,ict)
     &      ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &      ,n4s(ict),n4p(1,ict),n4p(-1,ict),n4p(0,ict)
c    &      ,n4d(0,ict),(n4d(m,ict),n4d(-m,ict),m=2,1,-1)
     &      ,nsa(ict),npa(1,ict),npa(-1,ict),npa(0,ict)
     &      ,nda(0,ict),(nda(m,ict),nda(-m,ict),m=2,1,-1)
           endif
          endif
   10   continue

        nbas_tot=0
        do 24 ic=1,ncent
          ict=iwctype(ic)
          nbas_typ=           iabs(n1s(ict))+iabs(n2s(ict))+iabs(n3s(ict))+iabs(n4s(ict))+iabs(n5s(ict))
          nbas_tot=nbas_tot + iabs(n1s(ict))+iabs(n2s(ict))+iabs(n3s(ict))+iabs(n4s(ict))+iabs(n5s(ict))
          do 21 m=-1,1
            nbas_typ=nbas_typ + iabs(n2p(m,ict)) + iabs(n3p(m,ict)) + iabs(n4p(m,ict)) + iabs(n5p(m,ict))
   21       nbas_tot=nbas_tot + iabs(n2p(m,ict)) + iabs(n3p(m,ict)) + iabs(n4p(m,ict)) + iabs(n5p(m,ict))
          do 22 m=-2,2
            nbas_typ=nbas_typ + iabs(n3d(m,ict)) + iabs(n4d(m,ict)) + iabs(n5d(m,ict))
   22       nbas_tot=nbas_tot + iabs(n3d(m,ict)) + iabs(n4d(m,ict)) + iabs(n5d(m,ict))
          do 23 m=-3,3
            nbas_typ=nbas_typ + iabs(n4f(m,ict)) + iabs(n5f(m,ict))
   23       nbas_tot=nbas_tot + iabs(n4f(m,ict)) + iabs(n5f(m,ict))
          do 123 m=-4,4
            nbas_typ=nbas_typ + iabs(n5g(m,ict))
  123       nbas_tot=nbas_tot + iabs(n5g(m,ict))

   24     nbasis_ctype(ict)=nbas_typ
        if (numr.ne.-4) then ! JT: do not stop if numr=-4, basis information will be read in later
         if(nbas_tot.ne.nbasis) then
          write(6,'(''nbas_tot,nbasis='',9i4)') nbas_tot,nbasis
          stop 'nbas_tot not equal to nbasis'
         endif
        endif

        betaq=0
        do 25 ic=1,ncent
   25     betaq=betaq+znuc(iwctype(ic))
        betaq=betaq-nelec+1

       else

c Numerical radial basis

c For the moment read basis info in new format for nloc=-1 only, but
c eventually change to this always.  For the moment we assume that it is
c single dot (ncent=nctype=1) and set nbasis_ctype=nbasis.
c (Do same for nloc = -5, since this is almost the same as nloc = -1 -ACM)
        if((nloc.eq.-1).or.(nloc.eq.-5)) then

          nbasis_ctype(1)=nbasis
          write(6,'(''Reading iwrwf'')')
          read(5,*) ((iwrwf(ib,ict),ib=1,nbasis),ict=1,nctype)
          write(6,'(''Center'',i5,'' uses radial bas. fns:'',(100i3))')
     &    (ict,(iwrwf(ib,ict),ib=1,nbasis),ict=1,nctype)
          if(ndim.eq.3) then
            read(5,*) (l_bas(ib),ib=1,nbasis)
            write(6,'(/,''l for basis functions:'',/,100i5)') (l_bas(ib),ib=1,nbasis)
          endif
          read(5,*) (m_bas(ib),ib=1,nbasis)
          write(6,'(/,''m for basis functions:'',/,100i5)') (m_bas(ib),ib=1,nbasis)
          do 30 ib=1,nbasis
            if(abs(m_bas(ib)).gt.ML_BAS) then
              write(6,'(''abs(m_bas(ib)) > ML_BAS, m_bas(ib),ML_BAS='',9i3)') m_bas(ib),ML_BAS
              stop 'abs(m_bas(ib)) > ML_BAS'
            endif
            if(l_bas(ib).gt.ML_BAS) then
              write(6,'(''l_bas(ib) > ML_BAS, l_bas(ib),ML_BAS='',9i3)') l_bas(ib),ML_BAS
              stop 'l_bas(ib) > ML_BAS'
            endif
   30     continue

         else

        do 49 ict=1,nctype
          write(6,'(/,''center type'',i3)') ict
          read(5,*) n1s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &    ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &    ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &    ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &    ,(n6h(m,ict),n6h(-m,ict),m=5,1,-1),n6h(0,ict)
          write(6,'(''number of basis fns. for each l:'',100i2)')
     &    n1s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &    ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &    ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &    ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &    ,(n6h(m,ict),n6h(-m,ict),m=5,1,-1),n6h(0,ict)

          nbas=iabs(n1s(ict))
          n2s(ict)=0
          n3s(ict)=0
          n4s(ict)=0
          do 41 m=-1,1
            n3p(m,ict)=0
            n4p(m,ict)=0
            if(n2p(m,ict).ne.0 .and. ML_BAS.lt.1) stop 'ML_BAS in vmc.h too small'
   41       nbas=nbas+ iabs(n2p(m,ict))
          do 42 m=-2,2
            n4d(m,ict)=0
            if(n3d(m,ict).ne.0 .and. ML_BAS.lt.2) stop 'ML_BAS in vmc.h too small'
   42       nbas=nbas+ iabs(n3d(m,ict))
          do 43 m=-3,3
            if(n4f(m,ict).ne.0 .and. ML_BAS.lt.3) stop 'ML_BAS in vmc.h too small'
   43       nbas=nbas+ iabs(n4f(m,ict))
          do 44 m=-4,4
            if(n5g(m,ict).ne.0 .and. ML_BAS.lt.4) stop 'ML_BAS in vmc.h too small'
   44       nbas=nbas+ iabs(n5g(m,ict))
          do 45 m=-5,5
            if(n6h(m,ict).ne.0 .and. ML_BAS.lt.5) stop 'ML_BAS in vmc.h too small'
   45       nbas=nbas+ iabs(n6h(m,ict))
          do m=-6,6
            if(n7i(m,ict).ne.0 .and. ML_BAS.lt.6) stop 'ML_BAS in vmc.h too small'
            nbas=nbas+ iabs(n7i(m,ict))
          enddo
          do m=-7,7
            if(n8j(m,ict).ne.0 .and. ML_BAS.lt.7) stop 'ML_BAS in vmc.h too small'
            nbas=nbas+ iabs(n8j(m,ict))
          enddo
          do m=-8,8
            if(n9k(m,ict).ne.0 .and. ML_BAS.lt.8) stop 'ML_BAS in vmc.h too small'
            nbas=nbas+ iabs(n9k(m,ict))
          enddo
          do m=-9,9
            if(n10l(m,ict).ne.0 .and. ML_BAS.lt.9) stop 'ML_BAS in vmc.h too small'
            nbas=nbas+ iabs(n10l(m,ict))
          enddo
          do m=-10,10
            if(n11m(m,ict).ne.0 .and. ML_BAS.lt.10) stop 'ML_BAS in vmc.h too small'
            nbas=nbas+ iabs(n11m(m,ict))
          enddo
          do m=-11,11
            if(n12n(m,ict).ne.0 .and. ML_BAS.lt.11) stop 'ML_BAS in vmc.h too small'
            nbas=nbas+ iabs(n12n(m,ict))
          enddo
          do m=-12,12
            if(n13o(m,ict).ne.0 .and. ML_BAS.lt.12) stop 'ML_BAS in vmc.h too small'
            nbas=nbas+ iabs(n13o(m,ict))
          enddo
          if(nbas.gt.MBASIS_CTYPE) stop 'nbas > MBASIS_CTYPE'
          nbasis_ctype(ict)=nbas

          write(6,'(''Reading iwrwf'')')
          read(5,*) (iwrwf(ib,ict),ib=1,nbas)
          write(6,'(''Center'',i5,'' uses radial bas. fns:'',(100i3))')
     &    ict,(iwrwf(ib,ict),ib=1,nbas)

          do 46 ib=1,nbas
            if(iwrwf(ib,ict).gt.nrbas(ict)) then
              write(6,'(''ict,ib,iwrwf(ib,ict),nrbas(ict)'',9i3)')
     &        ict,ib,iwrwf(ib,ict),nrbas(ict)
              stop 'iwrwf(ib,ict) > nrbas(ict)'
            endif
   46     continue

          nrwf=1
          do 48 ib=2,nbas
            do 47 ib2=1,ib
   47         if(iwrwf(ib2,ict).eq.iwrwf(ib,ict)) goto 48
            nrwf=nrwf+1
   48     continue
   49     if(nrwf.gt.MRWF) stop 'nrwf > MRWF'

        nbas_tot=0
        do 55 ic=1,ncent
          ict=iwctype(ic)
          nbas_tot=nbas_tot+ iabs(n1s(ict))
          do 51 m=-1,1
   51       nbas_tot=nbas_tot+ iabs(n2p(m,ict))
          do 52 m=-2,2
   52       nbas_tot=nbas_tot+ iabs(n3d(m,ict))
          do 53 m=-3,3
   53       nbas_tot=nbas_tot+ iabs(n4f(m,ict))
          do 54 m=-4,4
   54       nbas_tot=nbas_tot+ iabs(n5g(m,ict))
          do 55 m=-5,5
   55       nbas_tot=nbas_tot+ iabs(n6h(m,ict))
        if(nbas_tot.ne.nbasis) stop 'nbas_tot not equal to nbasis'

      endif ! end of nloc if
      endif

      if(nloc.ne.-1 .and. numr.ne.-4 .and. nloc.ne.-5) then
      write(6,'(/,''center type'',(12i4))') (i,i=1,nctype)
      write(6,'(/,''1s'',t11,(12i4))') (n1s(i),i=1,nctype)
      if(numr.le.0) write(6,'(''2s'',t11,(12i4))') (n2s(i),i=1,nctype)
      write(6,'(''2px'',t11,(12i4))') (n2p(1,i),i=1,nctype)
      write(6,'(''2py'',t11,(12i4))') (n2p(-1,i),i=1,nctype)
      write(6,'(''2pz'',t11,(12i4))') (n2p(0,i),i=1,nctype)
      if(numr.le.0) then
        write(6,'(''3s'',t11,(12i4))') (n3s(i),i=1,nctype)
        write(6,'(''3px'',t11,(12i4))') (n3p(1,i),i=1,nctype)
        write(6,'(''3py'',t11,(12i4))') (n3p(-1,i),i=1,nctype)
        write(6,'(''3pz'',t11,(12i4))') (n3p(0,i),i=1,nctype)
      endif
      write(6,'(''3dzr'',t11,(12i4))') (n3d(0,i),i=1,nctype)
      write(6,'(''3dx2'',t11,(12i4))') (n3d(2,i),i=1,nctype)
      write(6,'(''3dxy'',t11,(12i4))') (n3d(-2,i),i=1,nctype)
      write(6,'(''3dxz'',t11,(12i4))') (n3d(1,i),i=1,nctype)
      write(6,'(''3dyz'',t11,(12i4))') (n3d(-1,i),i=1,nctype)
      if(numr.le.0) then
        write(6,'(''4s'',t11,(12i4))') (n4s(i),i=1,nctype)
        write(6,'(''4px'',t11,(12i4))') (n4p(1,i),i=1,nctype)
        write(6,'(''4py'',t11,(12i4))') (n4p(-1,i),i=1,nctype)
        write(6,'(''4pz'',t11,(12i4))') (n4p(0,i),i=1,nctype)
        write(6,'(''4dzr'',t11,(12i4))') (n4d(0,i),i=1,nctype)
        write(6,'(''4dx2'',t11,(12i4))') (n4d(2,i),i=1,nctype)
        write(6,'(''4dxy'',t11,(12i4))') (n4d(-2,i),i=1,nctype)
        write(6,'(''4dxz'',t11,(12i4))') (n4d(1,i),i=1,nctype)
        write(6,'(''4dyz'',t11,(12i4))') (n4d(-1,i),i=1,nctype)
       endif
       do 60 m=-3,3
   60     write(6,'(''4f('',i2,'')'',t11,(12i4))') m,(n4f(m,i),i=1,nctype)
      if(numr.le.0) then
        write(6,'(''5s'',t11,(12i4))') (n5s(i),i=1,nctype)
        write(6,'(''5px'',t11,(12i4))') (n5p(1,i),i=1,nctype)
        write(6,'(''5py'',t11,(12i4))') (n5p(-1,i),i=1,nctype)
        write(6,'(''5pz'',t11,(12i4))') (n5p(0,i),i=1,nctype)
        write(6,'(''5dzr'',t11,(12i4))') (n5d(0,i),i=1,nctype)
        write(6,'(''5dx2'',t11,(12i4))') (n5d(2,i),i=1,nctype)
        write(6,'(''5dxy'',t11,(12i4))') (n5d(-2,i),i=1,nctype)
        write(6,'(''5dxz'',t11,(12i4))') (n5d(1,i),i=1,nctype)
        write(6,'(''5dyz'',t11,(12i4))') (n5d(-1,i),i=1,nctype)
        do 61 m=-3,3
   61     write(6,'(''5f('',i2,'')'',t11,(12i4))') m,(n5f(m,i),i=1,nctype)
       endif
       do 70 m=-4,4
   70     write(6,'(''5g('',i2,'')'',t11,(12i4))') m,(n5g(m,i),i=1,nctype)
       do 80 m=-5,5
   80     write(6,'(''6h('',i2,'')'',t11,(12i4))') m,(n6h(m,i),i=1,nctype)

      write(6,'(''sa'',t11,(12i4))') (nsa(i),i=1,nctype)
      write(6,'(''pxa'',t11,(12i4))') (npa(1,i),i=1,nctype)
      write(6,'(''pya'',t11,(12i4))') (npa(-1,i),i=1,nctype)
      write(6,'(''pza'',t11,(12i4))') (npa(0,i),i=1,nctype)
      write(6,'(''dzra'',t11,(12i4))') (nda(0,i),i=1,nctype)
      write(6,'(''dx2a'',t11,(12i4))') (nda(2,i),i=1,nctype)
      write(6,'(''dxya'',t11,(12i4))') (nda(-2,i),i=1,nctype)
      write(6,'(''dxza'',t11,(12i4))') (nda(1,i),i=1,nctype)
      write(6,'(''dyza'',t11,(12i4))') (nda(-1,i),i=1,nctype)

      write(6,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=1,nctype)
      write(6,*)

c The basis fns. in the LCAO coefs. are either in the atomic filling order
c or in order of increasing l.  The former is the order I used to use the
c latter is the order from GAMESS when doing all-electron calculations.
       if(numr.eq.0 .or. numr.eq.1) then
         call orb_loc_ana_original_order
        elseif(numr.eq.-1 .or. numr.eq.-2 .or. numr.eq.-3) then
         call orb_loc_ana_gamess_order
        else
         stop 'numr must be between -3 and 1'
       endif

      endif ! end of if(nloc.ne.-1)

      call object_modified ('nbasis_ctype') !JT

      write(6,*)
      write(6,'(''orbital coefficients'')')
      do 260 i=1,norb
        read(5,*) (coef(j,i,1),j=1,nbasis)
  260   write(6,'(12f10.6)') (coef(j,i,1),j=1,nbasis)

c Warning: Foll. "if" commented out till I modify change_input to work without the zex
c     if(inumr.ge.1) then
        write(6,'(''screening constants'')')
        read(5,*) (zex(i,1),i=1,nbasis)
        write(6,'(12f10.6)') (zex(i,1),i=1,nbasis)
        if (numr.ne.-4) then
         do 262 i=1,nbasis
  262      if(zex(i,1).le.0.d0) stop 'exponent zex should be > 0'
        endif
c     endif

      return
      end
c-----------------------------------------------------------------------

      subroutine orb_loc_ana_original_order
c Written by Cyrus Umrigar
c Reads in localized orbitals, in one or more of
c 1) a slater basis
c 2) an asymptotic basis
c 3) a gaussian basis
c For the p and d basis functions, I maintain the old order to avoid
c having to change old inputs.
c The order I read in the p functions is x,y,z, which is m=1,-1,0 and
c the order that I read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, which is m=0,2,-2,1,-1,
c in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.

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
c    &,l4s(nprime),l4p(nprime,-1:1)
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

      write(6,'(/,(12a10))') (lbasis(j),j=1,nbasis)
c     write(6,'(''orbital coefficients'')')
c     do 260 i=1,norb
c       read(5,*) (coef(j,i,1),j=1,nbasis)
c 260   write(6,'(12f10.6)') (coef(j,i,1),j=1,nbasis)

ccWarning: Foll. "if" commented out till I modify change_input to work without the zex
cc    if(inumr.ge.1) then
c       write(6,'(''screening constants'')')
c       read(5,*) (zex(i,1),i=1,nbasis)
c       write(6,'(12f10.6)') (zex(i,1),i=1,nbasis)
c       do 262 i=1,nbasis
c 262     if(zex(i,1).le.0.d0) stop 'exponent zex should be > 0'
cc    endif

      return
      end
c-----------------------------------------------------------------------

      subroutine orb_loc_ana_gamess_order
c Written by Cyrus Umrigar
c Reads in localized orbitals, in one or more of
c 1) a slater basis
c 2) an asymptotic basis
c 3) a gaussian basis
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

      write(6,'(/,(12a10))') (lbasis(j),j=1,nbasis)
c     write(6,'(''orbital coefficients'')')
c     do 260 i=1,norb
c       read(5,*) (coef(j,i,1),j=1,nbasis)
c 260   write(6,'(12f10.6)') (coef(j,i,1),j=1,nbasis)

ccWarning: Foll. "if" commented out till I modify change_input to work without the zex
cc    if(inumr.ge.1) then
c       write(6,'(''screening constants'')')
c       read(5,*) (zex(i,1),i=1,nbasis)
c       write(6,'(12f10.6)') (zex(i,1),i=1,nbasis)
c       do 262 i=1,nbasis
c 262     if(zex(i,1).le.0.d0) stop 'exponent zex should be > 0'
cc    endif

      return
      end
c-----------------------------------------------------------------------

      subroutine distinct_radial_bas
c Written by Cyrus Umrigar
c Determine distinct radial basis functions

      use all_tools_mod
      use atom_mod
      use coefs_mod
      use basis1_mod
      use numbas_mod
      use basis2_mod
      use wfsec_mod
      use contr3_mod
      implicit real*8(a-h,o-z)


c If an analytical basis is used find out how many distinct radial orbitals there are
c and store n and zex for each of these.
c This is done so that only the unique analytic radial wavefns. are evaluated in basis_fns.f and basis_fnse.f
c zex2 is used in func when optimizing the exponents.
c ib  goes up to the total number of basis functions
c ib2 goes up to the number of basis functions on each center type
c irb goes up to the number of radial basis functions on each center type
c Warning: If this is being called from fit then it is possible that exponents
c that start out being equal may not end that way.
c In particular this happens because GAMESS creates an s function out
c of d functions and these start off with the same exponent.
c So, if it is being called from fit then at present I check also that the l
c values are the same to guard against different l functions starting with same
c exponent and then becoming different.
c Warning:  This is a little inefficient and more importantly it is not fool proof.
c e.g. one may start with px,py,pz having the same exponents but then optimize
c them separately.  One could make it a bit more efficient by checking if nparme.ne.0
c but at present nparme is being set after the call to this routine.
c One could easily make it foolproof by just regarding each basis function
c as having a different radial function when nparme.ne.0, but this is even more
c inefficient during the optimization process.

      if(numr.le.0) then
!JT        write(6,'(''n_bas(ib),l_bas(ib)'',50(i2,x))') (n_bas(ib),l_bas(ib),ib=1,nbasis)
        ib=0
        ib4=0
        do 290 ic=1,ncent
          ict=iwctype(ic)
!JT          write(6,'(''n_bas(ib),l_bas(ib)'',50(i2,x))') (n_bas(i),l_bas(i),i=1,nbasis_ctype(ict))
          nrbas(ict)=0
          do 280 ib2=1,nbasis_ctype(ict)
            ib=ib+1
            if(ib.ne.ib2+ib4) stop 'ib .ne. ib2+ib4'
            do 270 ib3=1,ib2-1
              if((index(mode,'fit').eq.0 .and. n_bas(ib).eq.n_bas(ib3+ib4) .and. zex(ib,1).eq.zex(ib3+ib4,1)) .or.
     &        (l_bas(ib).eq.l_bas(ib3+ib4) .and. n_bas(ib).eq.n_bas(ib3+ib4) .and. zex(ib,1).eq.zex(ib3+ib4,1))) then
                irb=iwrwf2(ib3+ib4)
                goto 275
              endif
  270       continue
            nrbas(ict)=nrbas(ict)+1
            irb=nrbas(ict)
            if(irb.gt.MRWF) stop 'nbas > MRWF'
            n_bas2(irb,ict)=n_bas(ib)
            zex2(irb,ict,iwf)=zex(ib,iwf) ! JT: 1 -> iwf
  275       iwrwf2(ib)=irb
  280       iwrwf(ib2,ict)=irb
  290       ib4=ib4+nbasis_ctype(ict)
!JT        write(6,'(''nrbas(ict)='',40i4)') (nrbas(ict),ict=1,nctype)
!JT        write(6,'(''iwrwf2(ib)='',40i3)') (iwrwf2(ib),ib=1,nbasis)
!JT        write(6,'(''iwrwf(ib2,ict)='',40i3)') ((iwrwf(ib2,ict),ib2=1,nbasis_ctype(ict)),ict=1,nctype)
      endif

      call object_modified ('zex2')

c     if(numr.le.0) then
c       ib=0
c       ib4=0
c       do 290 ict=1,nctype
c         irb=0
c         do 290 ib2=1,nbasis_ctype(ict)
c           ib=ib+1
c           do 270 ib3=1,ib2-1
c             write(6,'(''ict,ib2,ib3,ib4,ib='',9i5)') ict,ib2,ib3,ib4,ib
c 270         if(n_bas(ib)-l_bas(ib).eq.n_bas(ib3+ib4)-l_bas(ib3+ib4)
c    &        .and. zex(ib,1).eq.zex(ib3+ib4,1)) goto 280
c           irb=irb+1
c           n_bas2(irb,ict)=n_bas(ib)
c           zex2(irb,ict,1)=zex(ib,1)
c           iwrwf2(ib)=irb
c 280       iwrwf(ib2,ict)=irb
c           if(irb.gt.MRWF) stop 'nbas > MRWF'
c           nrbas(ict)=irb
c 290       ib4=ib4+nbasis_ctype(ict)
c     endif

      return
      end
c-----------------------------------------------------------------------

      subroutine read_orb_loc_num
c Reads in numerical orbitals and V_ext on a grid.
c At present V_ext that is read in is not used.
      use coefs_mod
      implicit real*8(a-h,o-z)

      common /orbital_num/ orb_num(4,MGRID_ORB,MGRID_ORB,MORB_OCC),xorb_grid(MGRID_ORB),yorb_grid(MGRID_ORB)
     &,sizex,sizey,hx,hy,hxi,hyi,ngrid_orbx,ngrid_orby,ict(6)

c     dimension v_ext(MGRID_ORB,MGRID_ORB)

      open(4,file='orbitals_num',form='formatted',status='old')

      read(4,*) ngrid_orbx,ngrid_orby,sizex,sizey,norba
      if(ngrid_orbx.gt.MGRID_ORB .or. ngrid_orby.gt.MGRID_ORB) then
          write(6,'(a,i10,a,i10,a)') 'ngrid_orbx=',ngrid_orbx,' ngrid_orby=',ngrid_orby,' > MGRID_ORB'
          stop 'ngrid_orbx,ngrid_orby for numerical orbitals must be < MGRID_ORB'
      endif
      if(mod(ngrid_orbx,2).ne.1) stop 'ngrid_orbx must be odd in read_orb_loc_num'
      if(mod(ngrid_orby,2).ne.1) stop 'ngrid_orby must be odd in read_orb_loc_num'
      if(norba.lt.norb) stop 'norba not large enough in read_orb_loc_num'
      if(norba.gt.MORB) stop 'norba > MORB in read_orb_loc_num'

c     read(4,*) ((v_ext(i,j),i=1,ngrid_orbx),j=1,ngrid_orby)
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
      implicit real*8(a-h,o-z)
      parameter(MWORK=21)

      common /orbital_num/ orb_num(4,MGRID_ORB,MGRID_ORB,MORB_OCC),xorb_grid(MGRID_ORB),yorb_grid(MGRID_ORB)
     &,sizex,sizey,hx,hy,hxi,hyi,ngrid_orbx,ngrid_orby,ict(6)

      dimension bcxmin(MGRID_ORB),bcxmax(MGRID_ORB),bcymin(MGRID_ORB),bcymax(MGRID_ORB),
     &wk(MWORK*MGRID_ORB*MGRID_ORB)

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
        call r8mkbicubw(xorb_grid,ngrid_orbx,yorb_grid,ngrid_orby,orb_num(1,1,1,iorb),MGRID_ORB,
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

      subroutine read_orb_dot
c Written by A.D.Guclu, Feb 2004.
c Reads in 2-dimensional basis fns info for circular quantum dots.

      use all_tools_mod

      implicit real*8(a-h,o-z)


      common /dim/ ndim
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contrl_per/ iperiodic,ibasis
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

c do some debugging. not sure if all these are necessary to check
c but take no chance for now
      if(ndim.ne.2) stop 'ndim must be 2 for quantum dots'
      if(nforce.ne.1) stop 'nforce must be 1 for quantum dots'
      if(nloc.ne.-1) stop 'nloc must be -1 for quantum dots'
      if(numr.ne.0) stop 'numr must be 0 in read_orb_dot'
      if(inum_orb.ne.0) stop 'inum_orb must be 0 for quantum dots'

      if(ibasis.eq.3) then
        call read_orb_dot_fd
      elseif(ibasis.eq.4 .or. ibasis.eq.5) then
        call read_orb_dot_gauss
      else
        stop 'In read_orb_dot: only ibasis=3,4,5 allowed'
      endif
      return
      end
c-----------------------------------------------------------------------

      subroutine read_orb_dot_fd
c Written by A.D.Guclu, Feb 2004.
c Reads in quantum dot orbitals in Fock-Darwin basis set
c with quantum numbers n (quasi-Landau level) and m (angular mom.)
c The definition of "Landau level" is different in quantum Hall litt
c and quantum dot litt. For quantum dots we use the notation n_fd
c (for Fock-Darwin), for projected composite fermions we use n_cf
c just for convenience.
c For Fock-Darwin states, zex is more than an exponential parameter.
c It serves as a multiplicatif optimization factor for the spring constant.


      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /basis3/ n_fd(MBASIS),m_fd(MBASIS),n_cf(MBASIS),ncfmax
!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE)
     &,n2s(MCTYPE),n2p(-1:1,MCTYPE)
     &,n3s(MCTYPE),n3p(-1:1,MCTYPE),n3d(-2:2,MCTYPE)
     &,n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE),n4f(-3:3,MCTYPE)
     &,n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE),n5f(-3:3,MCTYPE)
     &,n5g(-4:4,MCTYPE)
     &,n6d(-2:2,MCTYPE),n6f(-3:3,MCTYPE),n6g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,n7g(-4:4,MCTYPE),n7h(-5:5,MCTYPE),n7i(-6:6,MCTYPE)
     &,n8i(-6:6,MCTYPE),n8j(-7:7,MCTYPE)
     &,n9k(-8:8,MCTYPE)
     &,n10l(-9:9,MCTYPE)
     &,n11m(-10:10,MCTYPE)
     &,n12n(-11:11,MCTYPE)
     &,n13o(-12:12,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /compferm/ emagv,nv,idot

c read the total number of "quasi-Landau" levels:
      read(5,*) nlandau
c next nlandau lines to read represent: n #m m1 m2 m3 ....
c for instance if we want landau levels with in the first LL; m=0,2,3, and
c in the second LL; m=1,2 then the input file should be:
c      0 3 0 2 3
c      1 2 1 2
      icount=0
      ncfmax=0
      do 10 in=1,nlandau
        read(5,*) n,nm,(m_fd(i),i=icount+1,icount+nm)
        do 5 i=icount+1,icount+nm
          n_fd(i)=n
          n_cf(i)=n
          if(m_fd(i).lt.0) n_cf(i)=n-m_fd(i)
          if(n_cf(i).gt.ncfmax) ncfmax=n_cf(i)
    5   enddo
        icount=icount+nm
   10 enddo
      if(icount.ne.nbasis) stop 'nbasis doesnt match the basis set'
      if(ncfmax.gt.MBASIS) stop 'ncfmax.gt.MBASIS. this is a problem in cbasis_fns.f'
c      if(ncfmax.gt.6 .and. idot.eq.3)
      if(ncfmax.gt.6)   ! idot not defined at this point
     &  write(6,'(''WARNING: landau levels 7 and 8 can cause numerical problems in projected cfs'')')

c read orbital coefficients
c      write(6,'(/,(12a10))') (n_fd(j),m_fd(j),j=1,nbasis)
      write(6,'(''orbital coefficients'')')
      do 20 iorb=1,norb
        read(5,*) (coef(j,iorb,1),j=1,nbasis)
        write(6,'(12f10.6)') (coef(j,iorb,1),j=1,nbasis)
   20 enddo

      write(6,'(''screening constants'')')
      read(5,*) (zex(i,1),i=1,nbasis)
      write(6,'(12f10.6)') (zex(i,1),i=1,nbasis)
      do 30 i=1,nbasis
c zex only used for idot=0, and should not be < 0. To be safe set it to 1.
c idot not defined at this point
        if(zex(i,1).le.0.d0) then
          write(6,'(''WARNING: exponent zex set to 1'')')
          zex(i,1)=1
        endif
   30 enddo


      return
      end

c-----------------------------------------------------------------------

      subroutine read_orb_dot_gauss
c Written by A.D.Guclu, Apr 2006.
c Reads in quantum dot orbitals in gaussian basis set
c the witdh of gaussians is given by zex*we

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /orbpar/ oparm(MOTYPE,MBASIS,MWF)
      common /optimo/ iwo(MORB,MOTYPE),nparmo(MOTYPE),nparmot,notype

      do it=1,notype
        read(5,*)  (oparm(it,ib,1),ib=1,nbasis)
        write(*,*) (oparm(it,ib,1),ib=1,nbasis)
      enddo

      do ib=1,nbasis
        if(oparm(3,ib,1).lt.0.d0) then
          write(6,'(''WARNING: exponent oparm(3,ib,1) set to 1'')')
          oparm(3,ib,1)=1
        endif
      enddo

      if(norb.ne.nbasis) stop
     &  'norb must be equal to nbasis in read_orb_dot_gauss'


c read orbital coefficients
c      write(6,'(/,(12a10))') (n_fd(j),m_fd(j),j=1,nbasis)
      write(6,'(''orbital coefficients'')')

c      if(norb.le.100) then     !  we don't need a linear combination of gaussians.
c        do 20 iorb=1,norb     !  turn this feature on if ever needed.
c          read(5,*) (coef(j,iorb,1),j=1,nbasis)
c   20   enddo
c      else
        write(6,'(''Assuming basis set=orbitals for 2D-gaussian orbitals'')')
        do 40 iorb=1,norb
          do 30 j=1,nbasis
            if(iorb.eq.j) then
              coef(j,iorb,1)=1
            else
              coef(j,iorb,1)=0
            endif
   30     enddo
   40   enddo

c      endif

      return
      end


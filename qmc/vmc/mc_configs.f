!=============================================================================
      subroutine mc_configs_read !JT
!-----------------------------------------------------------------------------
!=============================================================================

      use contr3_mod
      implicit none



      if(index(mode,'mpi').ne.0) then
         call mc_configs_read_mpi
      else
         call mc_configs_read_notmpi
      endif

      return !JT
!-----------------------------------------------------------------------
      entry mc_configs_write !JT

      if(index(mode,'mpi').ne.0) then
         call mc_configs_write_mpi
      else
         call mc_configs_write_notmpi
      endif

      return
      end

!-----------------------------------------------------------------------
      subroutine mc_configs_read_notmpi !JT
! If isite<=0 reads initial MC configuration from mc_configs_start
!         >=1 gets initial MC configuration by calling subroutine sites
! Write mc_configs_new at end of run to provide configurations for fit optimization

      use all_tools_mod
      use walkers_mod
      use atom_mod
      use config_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use pairden_mod
      implicit real*8(a-h,o-z)

      dimension nsite(ncent)

! Truncate fbias so that it is never negative, and the quantity
! sampled is never negative
!     fbias=dmin1(two,dmax1(zero,fbias))

      call alloc ('xold', xold, 3, nelec)

      if(irstar.ne.1) then

! if isite<=0 then get initial configuration from mc_configs_start
! if isite>=1 then get initial configuration from sites routine
        if(isite.ge.1) goto 393
        open(9,err=393,file='mc_configs_start')
        rewind 9
        read(9,*,end=393,err=393) ((xold(k,i),k=1,ndim),i=1,nelec)
        write(6,'(/,''initial configuration read from unit 9'')')
        close(9)
        goto 395

  393   l=0
        do 394 i=1,ncent
          nsite(i)=int(znuc(iwctype(i))+0.5d0)
          l=l+nsite(i)
          if(l.gt.nelec) then
            nsite(i)=nsite(i)-(l-nelec)
            l=nelec
          endif
  394     continue
        if(l.lt.nelec) nsite(1)=nsite(1)+(nelec-l)
        call sites(xold,nelec,nsite)
        call object_modified ('xold') ! JT
!JT        write(6,'(/,''initial configuration from sites'')')

  395   continue

! fix the position of electron i=ifixe for pair-density calculation:
        if(ifixe.gt.0) then
          do 398 k=1,ndim
  398       xold(k,ifixe)=xfix(k)
        endif


! If we are moving one electron at a time, then we need to initialize
! xnew, since only the first electron gets initialized in metrop
!       do 400 i=1,nelec
!         do 400 k=1,ndim
! 400       xnew(k,i)=xold(k,i)
      endif

! If nconf_new > 0 then we want to write nconf_new configurations from each processor for a future
! optimization or dmc calculation. So figure out how often we need to write a
! configuration to produce nconf_new configurations. If nconf_new = 0
! then set up so no configurations are written.
      if(nconf_new.ne.0) then
        open(7,form='formatted',file=trim(file_mc_configs_out))
        rewind 7
!       write(7,'(i5)') nconf_new
      endif
      return
!-----------------------------------------------------------------------
      entry mc_configs_write_notmpi !JT

! write out last configuration to unit mc_configs_start
      open(9,status='unknown',file='mc_configs_start')
      write(9,*) ((xold(k,i),k=1,ndim),i=1,nelec)
      close(9)

      return
      end

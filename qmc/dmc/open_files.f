      subroutine open_files
! Written by Cyrus Umrigar
! Open files that are different for serial and parallel runs
! Also do other things that differ for serial and parallel runs.

! Diffusion Monte Carlo and wave function optimization
! Warning: Check that we are fitting to the right linear combination of energy and variance.
! Used both for all-electon and 1-electron move versions.

      use all_tools_mod
      use walkers_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use contr3_mod
      use config_dmc_mod
      use contrldmc_mod
      implicit real*8(a-h,o-z)

      call get_initial_walkers

      if(irstar.ne.1) then
        if(ipr.gt.-2) then
!         open(11,file='walkalize',status='new')
          open(11,file='walkalize')
          rewind 11
!         write(11,*)
!    &    'Move line 2*nstep*nblkeq+1 here and delete this line'
          write(11,'(i3,'' nblkeq to be added to nblock at file end'')')
     &    nblkeq
        endif
      endif

! If nconf_new > 0 then we want to write nconf_new configurations from each processor for a future
! optimization or dmc calculation. So figure out how often we need to write a
! configuration to produce nconf_new configurations. If nconf_new = 0
! then set up so no configurations are written.
      if(nconf_new.eq.0) then
!       ngfmc=2*nstep*nblk
       else
!       ngfmc=max(1,(nstep*nblk)*nconf/nconf_new)
        open(7,form='formatted',file='mc_configs_new')
        rewind 7
      endif

      return
      end

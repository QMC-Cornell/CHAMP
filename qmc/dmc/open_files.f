      subroutine open_files
c Written by Cyrus Umrigar
c Open files that are different for serial and parallel runs
c Also do other things that differ for serial and parallel runs.

c Diffusion Monte Carlo and wave function optimization
c Warning: Check that we are fitting to the right linear combination of energy and variance.
c Used both for all-electon and 1-electron move versions.

      use all_tools_mod
      use walkers_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use contr3_mod
      use config_dmc_mod
      use contrldmc_mod
      implicit real*8(a-h,o-z)


!JT      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_global,nconf_new,isite,idump,irstar
!JT      common /contr3/ mode
!JT      common /contrldmc/ tau,rttau,taueff(MFORCE),tautot,nfprod,idmc,ipq
!JT     &,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
!JT      common /config_dmc/ xoldw(3,MELEC,MWALK,MFORCE),voldw(3,MELEC,MWALK,MFORCE),
!JT     &psidow(MWALK,MFORCE),psijow(MWALK,MFORCE),peow(MWALK,MFORCE),peiow(MWALK,MFORCE),d2ow(MWALK,MFORCE)

      call get_initial_walkers

      if(irstar.ne.1) then
        if(ipr.gt.-2) then
c         open(11,file='walkalize',status='new')
          open(11,file='walkalize')
          rewind 11
c         write(11,*)
c    &    'Move line 2*nstep*nblkeq+1 here and delete this line'
          write(11,'(i3,'' nblkeq to be added to nblock at file end'')')
     &    nblkeq
        endif
      endif

c If nconf_new > 0 then we want to write nconf_new configurations from each processor for a future
c optimization or dmc calculation. So figure out how often we need to write a
c configuration to produce nconf_new configurations. If nconf_new = 0
c then set up so no configurations are written.
      if(nconf_new.eq.0) then
c       ngfmc=2*nstep*nblk
       else
c       ngfmc=max(1,(nstep*nblk)*nconf/nconf_new)
        open(7,form='formatted',file='mc_configs_new')
        rewind 7
      endif

      return
      end

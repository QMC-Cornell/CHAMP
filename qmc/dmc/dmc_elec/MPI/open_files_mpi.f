      subroutine open_files_mpi
c Written by Cyrus Umrigar
c Open files that are different for serial and parallel runs
c Also do other things that differ for serial and parallel runs, such as setting random number seeds

# if defined (MPI)
      use all_tools_mod
      use dmc_mod
      use walkers_mod
      use dets_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use contr3_mod
      use periodic_mod
      use config_dmc_mod
      use contrldmc_mod
      implicit real*8(a-h,o-z)

      character*27 fmt
      character*20 filename

!JT      common /dim/ ndim

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_global,nconf_new,isite,idump,irstar
!JT      common /contr3/ mode
!JT      common /contrldmc/ tau,rttau,taueff(MFORCE),tautot,nfprod,idmc,ipq
!JT     &,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
!JT      common /config_dmc/ xoldw(3,MELEC,MWALK,MFORCE),voldw(3,MELEC,MWALK,MFORCE),
!JT     &psidow(MWALK,MFORCE),psijow(MWALK,MFORCE),peow(MWALK,MFORCE),peiow(MWALK,MFORCE),d2ow(MWALK,MFORCE)
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn

!JT      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
!JT     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
!JT     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
!JT     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
!JT     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
!JT     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
!JT     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
!JT     &,znuc_sum,znuc2_sum,vcell,vcell_sim
!JT     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
!JT     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
!JT     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

c set the random number seed, setrn already called in read_input
      if(irstar.ne.1) then
        do 95 id=1,(ndim*nelec)*idtask
  95      rnd=rannyu(0)
      endif

      if(ipr.gt.-2 .and. irstar.eq.0) then
        if(idtask.le.9) then
          write(filename,'(''walkalize.'',i1)') idtask
         elseif(idtask.le.99) then
          write(filename,'(''walkalize.'',i2)') idtask
         elseif(idtask.le.999) then
          write(filename,'(''walkalize.'',i3)') idtask
         elseif(idtask.le.9999) then
          write(filename,'(''walkalize.'',i4)') idtask
         elseif(idtask.le.99999) then
          write(filename,'(''walkalize.'',i5)') idtask
         elseif(idtask.le.999999) then
          write(filename,'(''walkalize.'',i6)') idtask
         else
          stop 'idtask > 999999'
        endif
        open(11,file=filename)
        rewind 11
c       write(11,*) 'Move line nstep*(2*nblkeq+nblk)+1 here and delete this line'
        write(11,'(i3,'' nblkeq to be added to nblock at file end'')') nblkeq
      endif

      call get_initial_walkers

c initialize sums and averages and reset nconf_global if there is one global population on all processors
      if(irstar.ne.1 .and. (mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3')) then
        call my_second(1,'zeres0')
        call zeres0_dmc
c This is now done in read_input
c       if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') then
c         nconf_global=nconf_global*nproc
c       endif
      endif

c If nconf_new > 0 then we want to dump configurations for a future
c optimization or dmc calculation. So figure out how often we need to write a
c configuration to produce nconf_new configurations. If nconf_new = 0
c then set up so no configurations are written.
      if(nconf_new.eq.0) then
c       ngfmc=2*nstep*nblk
       else
c       ngfmc=max(1,(nstep*nblk)*nconf/nconf_new)
        if(idtask.lt.10) then
          write(filename,'(i1)') idtask
         elseif(idtask.lt.100) then
          write(filename,'(i2)') idtask
         elseif(idtask.lt.1000) then
          write(filename,'(i3)') idtask
         elseif(idtask.lt.10000) then
          write(filename,'(i4)') idtask
         else
          stop 'idtask>=10000'
        endif
        filename='mc_configs_new'//filename(1:index(filename,' ')-1)
        open(7,form='formatted',file=filename)
        rewind 7
c       write(7,'(i5)') nconf_new
      endif

# endif
      return
      end

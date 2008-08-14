      subroutine opt_wf
c Written by Cyrus Umrigar

c Variational Monte Carlo -- both Metropolis and wave function optimization

c Warning: Check that we are fitting to the right linear combination of energy and variance.

c Various types of Metropolis moves can be done, including a few
c versions of directed Metropolis in spherical polar coordinates.
c Also, one or all electrons can be moved at once.

      use all_tools_mod    ! JT
      use montecarlo_mod   ! JT
      use optimization_mod ! JT
      use eloc_mod ! JT
      use mpi_mod ! JT

      implicit real*8(a-h,o-z)

!JT      parameter (one=1.d0,four=4.d0)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_global,nconf_new,isite,idump,irstar
      common /config/ xold(3,MELEC),xnew(3,MELEC),vold(3,MELEC)
     &,vnew(3,MELEC),psi2o(MFORCE),psi2n(MFORCE),eold(MFORCE),enew(MFORCE)
     &,peo,pen,peio,pein,tjfn,tjfo,psido,psijo
     &,rmino(MELEC),rminn(MELEC),rvmino(3,MELEC),rvminn(3,MELEC)
     &,rminon(MELEC),rminno(MELEC),rvminon(3,MELEC),rvminno(3,MELEC)
     &,nearesto(MELEC),nearestn(MELEC),delttn(MELEC)
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /doefp/ nefp
      common /contrl_opt2/ igradhess,iadd_diag_opt
      common /gradhess/ grad(MPARM),grad_var(MPARM),hess(MPARM,MPARM),hess_var(MPARM,MPARM),gerr(MPARM),
     &add_diag(3),energy(3),energy_sigma(3),energy_err(3),force(3),force_err(3),
     &eig_min,eig_max,p_var,tol_energy,nopt_iter,nblk_max

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap  !JT
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt


      dimension ene_var(3)

      call my_second(0,'begin ')

      pi=four*datan(one)

c     open(8,form='formatted',file='tape8')
c     rewind 8

c     mode='vmc         '
!JT      call read_input
c If add_diag(1) <= 0 turn OFF add_diag optimization and use fixed, add_diag
c equal to abs(add_diag) specified in input.  Usually you do NOT want to do this
c because optimization of add_diag adds only about 20% to the time since the
c number of blocks used for the correlated sampling is one-tenth of those used
c for calculating the gradient vector and hessian, overlap and hamiltonian matrices,
c and because optimization of add_diag makes the method totally stable.
      if(add_diag(1).gt.0) then
        iadd_diag_opt=1
       else
        add_diag(1)=abs(add_diag(1))
        iadd_diag_opt=0
      endif

c     if(nefp.gt.0) call readbas

      if(nforce.gt.1) then
c force parameters
        call readforce
c parameters for secondary geometry wave function
        call wf_secondary
       else
        nwftype=1
        iwftype(1)=1
      endif

c Initialize starting MC configuration
      call mc_configs_read

c If we are moving one electron at a time, then we need to initialize
c xnew, since only the first electron gets initialized in metrop
      if(irstar.ne.1) then
        do 400 i=1,nelec
          do 400 k=1,ndim
  400       xnew(k,i)=xold(k,i)
      endif

c Do wavefunction optimization if nopt_iter>0
      if(nopt_iter.gt.0) then

        if(idtask.eq.0) then
          open(2,file='wavefn_new')
         else
          open(2,file='/dev/null')
        endif

c Copy wavefn. for first value of add_diag to other values for use in correlated sampling
        nforce=3
        nwftype=3
        call wf_copy2
        call wf_copy
c initialize asymptotic values related to scalek(2) and scalek(3)
        call set_scale_dist(ipr,2)
        call set_scale_dist(ipr,3)

c Do optimization iterations
c increase_nblk is the number of consecutive optim. iterations for which nblk has not been increased
        increase_nblk=0
        energy_plus_err_best=1.d99
        do 430 iopt_iter=1,nopt_iter

c Do vmc computing gradient and hessian but no correlated sampling
          igradhess=1
          nforce=1
          nwftype=1
          iadd_diag_loop1=0
          iadd_diag_loop2=0
  405     call vmc

          ene_var(1)=(1-p_var)*energy(1)+p_var*energy_sigma(1)**2

          call optim_options

c Check if wavefunction got signficantly worse.
c Although we have overwritten grad and hess with new ones obtained from a worse wave function, we
c have not overwritten grad_sav and hess_sav and the wavefunction.  So, go back to the previous wf,
c increase add_diag(1) and use the previous grad and hess to create a new wf.
          if(iopt_iter.ge.2 .and. (energy_sigma(1).gt.(2.5d0-p_var)*energy_sigma_sav
     &    .or. energy(1)-energy_sav.gt.3*(1+p_var)*(sqrt(energy_err(1)**2+energy_err_sav**2)))) then
            iadd_diag_loop1=iadd_diag_loop1+1
            if(iadd_diag_loop1.gt.6) stop 'sigma went up considerably and iadd_diag_loop1>6'
            write(6,'(/,''Going back to previous wavefn. to generate new grad, hess, ham, ovlp'',/)')
            call wf_restore
            if(iadd_diag_opt.ne.0) then
              add_diag(1)=100*add_diag(1)
             else
              add_diag(1)=1.d-6*100**iadd_diag_loop2
              write(6,'(''Temporarily change add_diag from 0 to 1.d-6*100^iadd_diag_loop2 and back to 0 after calling new_param'')')
            endif
            write(6,'(''Energy_sigma got significantly worse in grad_hess call to vmc so increase add_diag(1) to'',
     &      1pd11.4)') add_diag(1)
            write(6,'(''add_diag_log_min,add_diag_min='',f7.3,1p,d9.2)')
     &      dlog10(add_diag(1)),add_diag(1)
            call new_param(1,1,1,iflag1)
            if(iadd_diag_opt.eq.0) add_diag(1)=0
c just in case mc config is in crazy place, reset mc_configs by calling sites
            isite=1
            call mc_configs_read
            iadd_diag_loop2=iadd_diag_loop2+1
            if(iadd_diag_loop2.gt.6) then
              write(6,'(''iadd_diag_loop2>4'')')
              stop 'iadd_diag_loop2>6'
            endif
            goto 405
          endif
          iadd_diag_loop2=0

c Wavefn. got better or at least not significantly worse so save various quantities.
          energy_sav=energy(1)
          energy_sigma_sav=energy_sigma(1)
          energy_err_sav=energy_err(1)
          ene_var_sav=(1-p_var)*energy(1)+p_var*energy_sigma(1)**2

c If this is the best yet, save it.  Since we are primarily interested in the energy we always use
c that as part of the criterion.  By adding in energy_err we favor those iterations where the energy
c has a smaller error, either because of a reduction in sigma and Tcorr or because nblk is increasing.
c If p_var!=0 then we add that to the criterion too.
          energy_plus_err=energy(1)+3*energy_err(1)+p_var*energy_sigma(1)
          if(energy_plus_err.lt.energy_plus_err_best) then
            energy_plus_err_best=energy_plus_err
            call wf_best_save
          endif

c Save wavefn for iadd_diag=1
          call wf_save

c Setup and save Hamiltonian and overlap matrices for linear method, and, gradient and Hessian for Newton.
c Normalize them and save the normalization for transforming the parameter variations.
c Save grad_sav,hess_sav,ham_sav,ovlp_sav,renorm_ovlp
c In grad_sav and hess_sav, use appropriate linear combination of energy and variance.
c Evaluate the eigenvalues of the Hessian of the objective function (linear comb of energy and variance).
          call ham_ovlp_grad_hess

c Turn correlated sampling on/off
          if(iadd_diag_opt.eq.1) then

c Reset add_diag
            add_diag(1)=1.d-8

c Make sure that the add_diag values are not tiny compared to eig_min
c Done in ham_ovlp_grad_hess now
c           add_diag(1)=max(add_diag(1),1.d-1*eig_min)
c           add_diag(2)=0.1d0*add_diag(1)
c           add_diag(3)=10.d0*add_diag(1)

            iadd_diag_loop3=0
            iadd_diag_loop4=0
c Compute 3 optimized wavefunctions with different values of add_diag
  408       add_diag(2)=0.1d0*add_diag(1)
            add_diag(3)=10*add_diag(1)
c Call with iadd_diag=1 last, because we get the parameters for the other values
c of iadd_diag by imcrementing the old parameters for iadd_diag=1.
            call new_param(2,1,0,iflag2)
            call new_param(3,0,0,iflag3)
            call new_param(1,0,0,iflag1)
            call set_scale_dist(ipr,1)
            call set_scale_dist(ipr,2)
            call set_scale_dist(ipr,3)
            if(iflag1.eq.1 .or. iflag2.eq.1 .or. iflag3.eq.1) then
              call wf_restore
c             call set_scale_dist(ipr,1)
              add_diag(1)=10**(iflag1+iflag2+iflag3)*add_diag(1)
              write(6,'(''change in params. too large or a(2) or b(2) < -scalek, add_diag(1) increased to'',1pd12.4)') add_diag(1)
              iadd_diag_loop4=iadd_diag_loop4+1
              if(iadd_diag_loop4.lt.9) then
                goto 408
               else
                stop 'iadd_diag_loop4 too large'
              endif
            endif

c Do vmc not computing gradient and hessian but do correlated sampling. Use smaller nblk
            igradhess=0
            nforce=3
            nwftype=3
            nblk_sav=nblk
            nblk=max(10,nblk/10)
            call vmc
            nblk=nblk_sav

c This is the objective function being optimized
            do 411 iadd_diag=1,3
  411         ene_var(iadd_diag)=(1-p_var)*energy(iadd_diag)+p_var*energy_sigma(iadd_diag)**2

c If ene_sigma(1) or energy_sigma(3) is much worse than before increase add_diag and loop back
c For variance minimization we allow it to be at most 1.5 times worse
c For energy   minimization we allow it to be at most 2.5 times worse
c and
c If energy(1) or energy(3) is much worse than before increase add_diag and loop back
c For variance minimization we allow it to be at most 6 std dev. worse
c For energy   minimization we allow it to be at most 3 std dev. worse

            if(energy_sigma(1).gt.(2.5d0-p_var)*energy_sigma_sav .or. energy_sigma(3).gt.(2.5d0-p_var)*energy_sigma_sav
     &      .or. energy(1)-energy_sav.gt.3*(1+p_var)*(sqrt(energy_err(1)**2+energy_err_sav**2))
     &      .or. energy(3)-energy_sav.gt.3*(1+p_var)*(sqrt(energy_err(3)**2+energy_err_sav**2))) then
              call wf_restore
c             call set_scale_dist(ipr,1)
              iadd_diag_loop3=iadd_diag_loop3+1
              add_diag(1)=100*add_diag(1)
              write(6,'(''Objective function went up too much, add_diag(1) increased to'',1pd12.4)') add_diag(1)
              if(iadd_diag_loop3.lt.9 .and. add_diag(1).lt.1.d10) then
c just in case mc config is in crazy place, reset mc_configs by calling sites
                isite=1
                call mc_configs_read
                goto 408
               else
                write(6,'(''iadd_diag_loop3 or add_diag(1) too large'')')
                stop 'iadd_diag_loop3 or add_diag(1) too large'
              endif
            endif

c Find optimal a_diag
            add_diag1_sav=add_diag(1)
c           call quad_min(energy,add_diag,force,force_err,3,eig_min,eig_max)
c           call quad_min(energy_sav,energy_err_sav,energy,energy_err,force,force_err,ene_var,add_diag,3,eig_min,eig_max,p_var)
c           call quad_min(energy_sav,energy_err_sav,ene_var,3)
            call quad_min(ene_var,3)


c Restore wavefn for iadd_diag=1 before updating with optimal add_diag
            call wf_restore

          endif
c end of 1st if(iadd_diag_opt.eq.1) then

c Calculate optimized wavefunction for the optimal a_diag
c If add_diag(1).ge.0.1d0*add_diag1_sav then use ipr_new=2 to put _new subscript when writing
c parameters regardless of the value of iflag1
          if(add_diag(1).ge.0.1d0*add_diag1_sav .or. iadd_diag_opt.eq.1) then
            call new_param(1,0,2,iflag1)
           else
            call new_param(1,0,1,iflag1)
          endif
          call set_scale_dist(ipr,1)

          if(iadd_diag_opt.eq.1) then

c It is still possible that move is too large if add_diag_log_min was set to
c a smaller value than add_diag_log(1)-1 in quad_min.f.  In quad_min.f we allow values as
c small as add_diag_log(1)-2.d0), i.e. setting add_diag(1) to 0.01*add_diag(1).
c So check iflag1 if this is the case.
            if(iflag1.ne.0 .and. add_diag(1).lt.0.1d0*add_diag1_sav) then
              write(6,'(''Increasing add_diag(1) because iflag.ne.0'')')
              call wf_restore
c             add_diag(1)=10*add_diag(1)
              add_diag(1)=0.1d0*add_diag1_sav
              call new_param(1,0,1,iflag1)
              call set_scale_dist(ipr,1)
              if(iflag1.ne.0) stop 'iflag.ne.0 should not occur here'
            endif

c Check if optimized wavefunction energies for 3 different values of a_diag
c are less than the tolerance
            energy_min=9.d90
            energy_max=-9.d90
            do 420 i=1,3
              energy_min=min(energy_min,energy(i))
  420         energy_max=max(energy_max,energy(i))
            e_diff=energy_max-energy_min
            write(6,'(''iopt_iter,e_diff='',i4,d12.4)') iopt_iter,e_diff

c Increase nblk if near convergence to value needed to get desired statistical error
            increase_nblk=increase_nblk+1
            add_diag_ratio=(add_diag(1)+2*max(-eig_min,0.d0))/add_diag(1)
            if(add_diag_ratio*e_diff.lt.tol_energy .and.
     &      force_err(2).lt.tol_energy .and. force_err(3).lt.tol_energy) then
c             nblk_tmp=nblk*min(10.d0,max(1.d0,(energy_err_sav/tol_energy)**2))
              nblk_tmp=nblk*max(1.d0,(energy_err_sav/tol_energy)**2)
              nblk_tmp=int_round(nblk_tmp)
              nblk_tmp=min(nblk_tmp,nblk_max)
              if(nblk_tmp.gt.nblk) then
                increase_nblk=0
                nblk=nblk_tmp
                write(6,'(''nblk reset to'',i8,9d12.4)') nblk,energy_err(1),tol_energy
              endif
              write(6,'(''The difference in the energies for different add_diag is converged to'',d12.4)') tol_energy
              goto 440
            endif
c At minimum increase nblk by a factor of 2 every other iteration
            if(increase_nblk.eq.2 .and. nblk.lt.nblk_max) then
              increase_nblk=0
              nblk=min(2*nblk,nblk_max)
              write(6,'(''nblk reset to'',i8,9d12.4)') nblk
            endif

          endif
c end of 2nd if(iadd_diag_opt.eq.1) then

        write(6,*)
  430   write(2,*)
  440   nforce=1
      endif

c Do vmc with either input params. if nopt_iter=0 or with optim. param if nopt_iter>0
c If nopt_iter is negative then we want to calculate the gradient and Hessian for later
c analysis outside the program.
      if(nopt_iter.lt.0) then
        igradhess=1
       else
        igradhess=0
      endif
      call vmc

c if nopt_iter=-1 we are doing one iteration for predicting new parameters
c without any checking of whether they are good or not.
      if(nopt_iter.lt.0) then
        call optim_options
        call ham_ovlp_grad_hess
        call new_param(1,1,1,iflag1)
c       call set_scale_dist(ipr,1)
      endif

c If this is the best yet, save it.  Since we are primarily interested in the energy we use
c that as part of the criterion.  By adding in energy_err we favor those iterations where the energy
c has a smaller error, either because of a reduction in sigma and Tcorr or because nblk is increasing.
c If p_var!=0 then we add that to the criterion too.
      if(nopt_iter.ne.0) then
        energy_plus_err=energy(1)+3*energy_err(1)+p_var*energy_sigma(1)
        write(6,'(''energy_plus_err,energy_plus_err_best'',2f10.5)') energy_plus_err,energy_plus_err_best
        if(energy_plus_err.lt.energy_plus_err_best) then
          energy_plus_err_best=energy_plus_err
          call wf_best_save
        endif
        call wf_best_write
      endif

      return
      end
c-----------------------------------------------------------------------
      function int_round(i)
c Find the nearest larger "round" integer
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)

      ipow=10**int((dlog10(dfloat(i))))
      rest=dfloat(i)/ipow
      if(rest.le.2.d0) then
        j=2
       elseif(rest.le.4.d0) then
        j=4
       elseif(rest.le.5.d0) then
        j=5
       elseif(rest.le.8.d0) then
        j=8
       else
        j=10
      endif

      int_round=j*ipow
c     write(6,*) i,rest,j,ipow,int_round

      return
      end

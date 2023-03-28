      subroutine vmc
! Written by Cyrus Umrigar

! Program to do variational monte carlo calculations on
! atoms and molecules.
! Various types of Metropolis moves can be done, including a few
! versions of directed Metropolis in spherical polar coordinates.
! Also, one or all electrons can be moved at once.
! Currently this program contains
! 1s, 2s, 2p, 3s, 3p, 3d, 4s,  and 4p  Slater or gaussian basis functions.
! and sa, pa, da asymptotic functions
! MELEC must be >= number of electrons

! The dimension of the Slater matrices in determinant is taken from
! above assuming equal number of up and down spins.
! That is the Slater matrices are dimensioned (MELEC/2)**2.
! MELEC would have to be correspondingly larger if spin
! polarized calculations were attempted.

      use all_tools_mod
      use main_menu_mod
      use intracule_mod
      use montecarlo_mod
      use eloc_mod
      use average_mod
      use control_mod
      use print_mod
      use deriv_exp_mod
      use walkers_mod
      use atom_mod
      use config_mod
      use dorb_mod
      use coefs_mod
      use dets_mod
      use basis1_mod
      use contrl_mod
      use const_mod
      use const2_mod
      use dim_mod
      use numbas_mod
      use contr2_mod
      use forcepar_mod
      use wfsec_mod
      use pseudo_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use bparm_mod
      use contr3_mod
      use qua_mod
      use forcest_mod
      use pars_mod
      use jaspar1_mod
      use jaspar2_mod
      use distance_mod
      implicit real*8(a-h,o-z)
      integer fflag
      character*80 fmt

      common /fflags/ fflag

      common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4

! common block variables:

!   /const/
!        nelec  = number of electrons
!        pi     = 3.14159...
!        hb     = hbar**2/(2m)
!        delta  = side of box in which metropolis steps are made
!        deltai = 1/delta
!        fbias  = force bias parameter
!   /contrl/
!        nstep  = number of metropolis steps/block
!        nblk   = number of blocks od nstep steps after the
!                equilibrium steps
!        nblkeq = number of equilibrium blocks
!        nconf_global  = target number of mc configurations (dmc only)
!        nconf_new = number of mc configurations generated for optim and dmc
!        idump  =  1 dump out stuff for a restart
!        irstar =  1 pick up stuff for a restart
!   /config/
!        xold   = current position of the electrons
!        xnew   = new position after a trial move
!        vold   = grad(psi)/psi at current position
!        vnew   = same after trial move
!        psi2o  = psi**2 at current position
!        psi2n  = same after trial move
!        eold   = local energy at current position
!        enew   = same after trial move
!        peo    = local potential at current position
!        pen    = same after trial move
!        peio   = local e-e interaction potential at current position
!        pein   = same after trial move
!        tjfo   = Jackson Feenberg kinetic energy at current position
!        tjfn   = same after trial move
!        psido  = determinantal part of wave function
!        psijo  = log(Jastrow)
!        pot_ee_old(ielec) = interaction potential felt by ith electron
!                              at current position
!        pot_ee_new(ielec) = same after trial move
!   /coefs/
!        coef   = read in coefficients of the basis functions
!                 to get the molecular orbitals used in determinant
!        nbasis = number of basis functions read in
!   /basis/
!        ncent  = number of centers
!        zex    = screening constants for each basis function
!        cent   = positions of each center
!        pecent = potential energy of the centers
!        znuc   = charge of the nuclei (centers)
!        n1s    = number of 1s functions at each center
!        n2s    = number of 2s functions at each center
!        n2p    = number of 2p functions of each type at each center
!        n3s    = number of 3s functions at each center
!        n3p    = number of 3p functions of each type at each center
!        n3dzr  = number of z**2-r**2 d functions at each center
!        n3dx2  = number of x**2-y**2 d functions at each center
!        n3dxy  = number of xy d functions at each center
!        n3dxz  = number of xz d functions at each center
!        n3dyz  = number of yz d functions at each center
!        n4s    = number of 4s functions at each center
!        n4p    = number of 4p functions of each type at each center
!   /dets/
!        cdet   = coefficients of the determinants
!        ndet   = number of determinants of molecular orbitals
!                 used
!        nup    = number of up spin electrons
!        ndn    = number of down spin electrons
!   /jaspar/
!        Jastrow function is dexp(cjas1*rij/(1+cjas2*rij)) if ijas=1

!   Other variables main program
!        title  = title of run
!        date   = date of run
!        eunit  = energy units
!        sitsca = scaling factor to set up initial configuration of
!                 electrons on sites
!        nsite  = number of electrons to put on each site initially
!        isite  = flag if 1 then take initial configuration from
!                 sites routine


!      initial printing
       write(6,*)
       write(6,'(a)') '*********** START VMC CALCULATION  ***********'
       write(6,*)

       call object_average_request ('eloc_av')
       call object_variance_request ('eloc_av_var')
       call object_average_request ('eloc_sq_av')
       call object_error_request ('error_sigma')

      call print_list_of_averages_and_errors

      nparma_read=2+max(0,norda-1)
      nparmb_read=2+max(0,nordb-1)
      nparmc_read=nterms4(nordc)

! Temporary print out of Jastrow params.
!     do 156 iwf=1,nforce
!     write(6,'(''iwf='',i3)') iwf
!     do 152 ict=1,nctype
! 152   write(6,'(''a='',9f10.6)') (a4(i,ict,iwf),i=1,nparma_read)
!     do 154 isp=nspin1,nspin2b
! 154   write(6,'(''b='',9f10.6)') (b(i,isp,iwf),i=1,nparmb_read)
!     do 156 ict=1,nctype
! 156   write(6,'(''c='',9f10.6)') (c(i,ict,iwf),i=1,nparmc_read)

! If we are moving one electron at a time, then we need to initialize
! xnew, since only the first electron gets initialized in metrop
      call alloc ('xnew', xnew, 3, nelec)
      if(irstar.ne.1) then
        do 400 i=1,nelec
          do 400 k=1,ndim
  400       xnew(k,i)=xold(k,i)
      endif

! We need to set pot_ee_old(:)
!  also, zero out rshift(k,i,ic) to start with
      rshift(:,:,:) = 0.
      call distances(xold, pe, pei)
      pot_ee_old = pot_ee    ! this is an array assignment

! If nconf_new > 0 then we want to write nconf_new configurations from each processor for a future
! optimization or dmc calculation. So figure out how often we need to write a
! configuration to produce nconf_new configurations. If nconf_new = 0
! then set up so no configurations are written.
      if(nconf_new.eq.0) then
        ngfmc=2*nstep*nblk
       else
        ngfmc=max(1,(nstep*nblk)/nconf_new)
      endif

! zero out estimators and averages
      if(irstar.ne.1) then
!JT        call my_second(1,'zerest')
         call zerest
      endif

! check if restart flag is on. If so then read input from
! dumped data to restart

      if(irstar.eq.1) then
        open(10,file='restart_vmc',err=401,form='unformatted')
        goto 402
  401   stop 'restart_vmc empty: not able to restart'
  402   rewind 10
        call startr
        close(unit=10)
      endif

! Start equilibration steps
      call print_cpu_time_in_seconds ('Beginning of equilibration')

! if there are equilibrium steps to take, do them here
! skip equilibrium steps if restart run
! imetro = 1 original force-bias
!        = 5 spherical-polar with exponential/linear T
!        = 6 spherical-polar with slater T
!        = 7 spherical-polar with slater T
      if(nblkeq.ge.1.and.irstar.ne.1) then
        l=0

        l_equilibration = .true.
        do 420 i=1,nblkeq
          do 410 j=1,nstep
            current_walker_weight = 1
!           call object_modified ('current_walker_weight')
            l=l+1
            if(nloc.gt.0) call rotqua
            if(imetro.eq.1) then
              call metrop(l)
             elseif(imetro.eq.5) then
              call metrop_polar(l)
             elseif(imetro.eq.6) then
              call metrop_slat(l)
            endif
  410     continue
  420   call acuest
        l_equilibration = .false.

!       Equilibration steps done. Zero out estimators again.
        call print_cpu_time_in_seconds ('End       of equilibration')
        call zerest
      endif

! now do averaging steps
      l=0
      do 440 i=1,1000000  !JT
        block_iterations_nb = block_iterations_nb + 1  !JT
        do 430 j=1,nstep
          step_iterations_nb = step_iterations_nb + 1   !JT
          current_walker_weight = 1
!         call object_modified ('current_walker_weight')
          l=l+1
          if(nloc.gt.0) call rotqua
          if(imetro.eq.1) then
            call metrop(l)
           elseif(imetro.eq.5) then
            call metrop_polar(l)
           elseif(imetro.eq.6) then
             call metrop_slat(l)
          endif

!JT       call object_modified_by_index (xold_index)  !JT

!         accumulate data for averages and statistical errors
!         call compute_block_averages       !JT old
          call compute_averages_step        !JT new

          call compute_averages_walk_step   !JT

!         write out configuration for optimization/dmc/gfmc here
          if (nconf_new /= 0) then
           if(mod(l,ngfmc).eq.1 .or. ngfmc.eq.1) then
!           write(fmt,'(a1,i6,a21)')'(',ndim*nelec,'f14.8,i3,d12.4,f12.5)'
            write(fmt,'(''('',i6,''f13.8,i3,es12.4,f12.5,2es12.4)'')') ndim*nelec
            write(7,fmt) ((xold(k,jj),k=1,ndim),jj=1,nelec), &
     &      int(sign(1.d0,psido)),log(dabs(psido))+psijo,eold(1),psido,psijo
           endif
          endif

!         write out configurations in Scemama's format
          if (l_write_walkers) then
           if(mod(l,write_walkers_step) == 0) then
            do jj = 1, nelec
              if (l_write_walkers_modified_format) then
               write(file_walkers_out_unit,'(3(F11.7,X))') (xold(k,jj),k=1,ndim)
              else
               write(file_walkers_out_unit,'(3(F7.3,X))') (xold(k,jj),k=1,ndim)
              endif
            enddo
             write(file_walkers_out_unit,*) dexp(psijo)*psido
             if (l_write_walkers_modified_format) then
              write(file_walkers_out_unit,*) eold(1)
             endif
           endif
          endif

  430   continue


!      call compute_global_averages   !JT old
      call compute_averages_block   !JT new
      call compute_averages_walk_block   !JT

      call compute_covariances  !JT
      call compute_variances  !JT
      call compute_errors  !JT

!     write at each block
      call objects_print_at_each_block
      call routines_write_block !JT

!     compute averages and statitical errors
      call acuest


!     save walkers for DMC run
      if (l_generate_walkers_from_vmc) then
        call save_vmc_walkers_for_dmc
      endif

!     exit loop if nblk and threshold on statistical error reached
!JT      if (block_iterations_nb .ge. nblk .and. eloc_av_err .gt. 0 .and. eloc_av_err .le. error_threshold) then    !JT
!     eloc_av_err .ne. eloc_av_err is true if eloc_av_err is NaN
      if (block_iterations_nb .ge. nblk .and. ((eloc_av_err .le. error_threshold) .or. (eloc_av_err .ne. eloc_av_err))) then    !JT
        exit                              !JT
      endif                               !JT
  440 continue

!      if( nforce == 1) then
!         write(6, *) "fp: biased sigma", sigma
!         write(6, *) "fp: biased energy", eloc_av
!         write(6, *) "fp: biased energy error", sqrt(eloc_av_var)
!         write(6, *) "fp: biased sigma error", error_sigma
!      end if
!     reinitilization at the end of MC iterations
      block_nb = block_iterations_nb !JT save final block number
      step_iterations_nb  = 0   !JT
      block_iterations_nb = 0   !JT
      call reinit_averages_and_errors !JT
      call reinit_routines_write_block !JT

      call print_cpu_time_in_seconds ('End       of accumulation ')

! print out final results
      call finwrt

! Presently grad_hess_jas_fin needs to be called after finwrt
!     if(ijasderiv.eq.1 .and. idtask.eq.0) then
!       passes=dfloat(iblk)*dfloat(nstep)*dfloat(nproc)
!       efin=ecum(1)/passes
!       call grad_hess_fin(passes,efin)
!     endif

! Write out last MC configuration on mc_configs_start
      call mc_configs_write

! if dump flag is on then dump out data for a restart
      if(idump.eq.1) then
        open(10,form='unformatted',file='restart_vmc')
        rewind 10
        call dumper
        close(unit=10)
      endif

      close(unit=9)
!     close(unit=5)
!     close(unit=6)
      if(nconf_new.ne.0) close(unit=7)

      return
      end

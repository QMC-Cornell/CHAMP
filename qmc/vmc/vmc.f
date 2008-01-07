      subroutine vmc
c Written by Cyrus Umrigar

c Program to do variational monte carlo calculations on
c atoms and molecules.
c Various types of Metropolis moves can be done, including a few
c versions of directed Metropolis in spherical polar coordinates.
c Also, one or all electrons can be moved at once.
c Currently this program contains
c 1s, 2s, 2p, 3s, 3p, 3d, 4s,  and 4p  Slater or gaussian basis functions.
c and sa, pa, da asymptotic functions
c MELEC must be >= number of electrons
c MBASIS  must be >= number of basis functions
c MDET must be >= number of determinants used
c MCENT  must be >= number of centers used

c The dimension of the Slater matrices in determinant is taken from
c above assuming equal number of up and down spins.
c That is the Slater matrices are dimensioned (MELEC/2)**2.
c MELEC would have to be correspondingly larger if spin
c polarized calculations were attempted.

      use all_tools_mod        !JT
      use main_menu_mod        !JT
      use intracule_mod        !JT
      use montecarlo_mod       !JT
      use eloc_mod             !JT
      use average_mod          !JT
      use control_mod          !JT
      use print_mod            !JT
      use deriv_exp_mod        !JT
      use walkers_mod          !JT

      implicit real*8(a-h,o-z)
      integer fflag
      character*16 mode
      character*25 fmt

c     include '../fit/fit.h'

      common /dim/ ndim
      common /fflags/ fflag
      common /pars/ a00,a20,a21,eps_fock,c0000,c1110,c2000,
     &   xm1,xm2,xm12,xms,xma,Z
      common /rlobxy/ rlobx(nsplin), rloby(nsplin), rloby2(nsplin)
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contr3/ mode

      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fcum(MFORCE),fcm2(MFORCE)
      common /forcewt/ wsum(MFORCE),wcum(MFORCE)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /const2/ deltar,deltat
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /config/ xold(3,MELEC),xnew(3,MELEC),vold(3,MELEC)
     &,vnew(3,MELEC),psi2o(MFORCE),psi2n(MFORCE),eold(MFORCE),enew(MFORCE)
     &,peo,pen,peio,pein,tjfn,tjfo,psido,psijo
     &,rmino(MELEC),rminn(MELEC),rvmino(3,MELEC),rvminn(3,MELEC)
     &,rminon(MELEC),rminno(MELEC),rvminon(3,MELEC),rvminno(3,MELEC)
     &,nearesto(MELEC),nearestn(MELEC),delttn(MELEC)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE)
     &,n5f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar1/ cjas1(MWF),cjas2(MWF)
      common /jaspar2/ a1(MPARMJ,3,MWF),a2(MPARMJ,3,MWF)
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /bparm/ nspin2b,nocuspb
      common /estsum/ esum1,esum(MFORCE),pesum,peisum,tpbsum,tjfsum,r2sum,accsum
      common /estcum/ ecum1,ecum(MFORCE),pecum,peicum,tpbcum,tjfcum,r2cum,acccum,iblk
      common /est2cm/ ecm21,ecm2(MFORCE),pecm2,peicm2,tpbcm2,tjfcm2,r2cm2
      common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4
      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /doefp/ nefp

c common block variables:

c   /const/
c        nelec  = number of electrons
c        pi     = 3.14159...
c        hb     = hbar**2/(2m)
c        delta  = side of box in which metropolis steps are made
c        deltai = 1/delta
c        fbias  = force bias parameter
c   /contrl/
c        nstep  = number of metropolis steps/block
c        nblk   = number of blocks od nstep steps after the
c                equilibrium steps
c        nblkeq = number of equilibrium blocks
c        nconf  = target number of mc configurations (dmc only)
c        nconf_new = number of mc configurations generated for optim and dmc
c        idump  =  1 dump out stuff for a restart
c        irstar =  1 pick up stuff for a restart
c   /config/
c        xold   = current position of the electrons
c        xnew   = new position after a trial move
c        vold   = grad(psi)/psi at current position
c        vnew   = same after trial move
c        psi2o  = psi**2 at current position
c        psi2n  = same after trial move
c        eold   = local energy at current position
c        enew   = same after trial move
c        peo    = local potential at current position
c        pen    = same after trial move
c        peio   = local e-e interaction potential at current position
c        pein   = same after trial move
c        tjfo   = Jackson Feenberg kinetic energy at current position
c        tjfn   = same after trial move
c        psido  = determinantal part of wave function
c        psijo  = log(Jastrow)
c   /coefs/
c        coef   = read in coefficients of the basis functions
c                 to get the molecular orbitals used in determinant
c        nbasis = number of basis functions read in
c   /basis/
c        ncent  = number of centers
c        zex    = screening constants for each basis function
c        cent   = positions of each center
c        pecent = potential energy of the centers
c        znuc   = charge of the nuclei (centers)
c        n1s    = number of 1s functions at each center
c        n2s    = number of 2s functions at each center
c        n2p    = number of 2p functions of each type at each center
c        n3s    = number of 3s functions at each center
c        n3p    = number of 3p functions of each type at each center
c        n3dzr  = number of z**2-r**2 d functions at each center
c        n3dx2  = number of x**2-y**2 d functions at each center
c        n3dxy  = number of xy d functions at each center
c        n3dxz  = number of xz d functions at each center
c        n3dyz  = number of yz d functions at each center
c        n4s    = number of 4s functions at each center
c        n4p    = number of 4p functions of each type at each center
c   /dets/
c        cdet   = coefficients of the determinants
c        ndet   = number of determinants of molecular orbitals
c                 used
c        nup    = number of up spin electrons
c        ndn    = number of down spin electrons
c   /jaspar/
c        Jastrow function is dexp(cjas1*rij/(1+cjas2*rij)) if ijas=1

c   Other variables main program
c        title  = title of run
c        date   = date of run
c        eunit  = energy units
c        sitsca = scaling factor to set up initial configuration of
c                 electrons on sites
c        nsite  = number of electrons to put on each site initially
c        isite  = flag if 1 then take initial configuration from
c                 sites routine


!      initial printing
       write(6,*)
       write(6,'(a)') '*********** START VMC CALCULATION  ***********'
       write(6,*)

!     sigma
      call object_associate ('error_sigma', error_sigma) !JT
      call object_average_request ('eloc_sq_av') !JT
      call object_error_request ('error_sigma') !JT
!     variance on averaged local energy
      call object_associate ('eloc_av', eloc_av) !JT
      call object_associate ('eloc_av_var', eloc_av_var) !JT
      call object_variance_request ('eloc_av_var')

      call print_list_of_averages_and_errors

      nparma_read=2+max(0,norda-1)
      nparmb_read=2+max(0,nordb-1)
      nparmc_read=nterms4(nordc)

c Temporary print out of Jastrow params.
c     do 156 iwf=1,nforce
c     write(6,'(''iwf='',i3)') iwf
c     do 152 ict=1,nctype
c 152   write(6,'(''a='',9f10.6)') (a4(i,ict,iwf),i=1,nparma_read)
c     do 154 isp=nspin1,nspin2b
c 154   write(6,'(''b='',9f10.6)') (b(i,isp,iwf),i=1,nparmb_read)
c     do 156 ict=1,nctype
c 156   write(6,'(''c='',9f10.6)') (c(i,ict,iwf),i=1,nparmc_read)

c If we are moving one electron at a time, then we need to initialize
c xnew, since only the first electron gets initialized in metrop
      if(irstar.ne.1) then
        do 400 i=1,nelec
          do 400 k=1,ndim
  400       xnew(k,i)=xold(k,i)
      endif

c If nconf_new > 0 then we want to dump configurations for a future
c optimization or dmc calculation. So figure out how often we need to write a
c configuration to produce nconf_new configurations. If nconf_new = 0
c then set up so no configurations are written.
      if(nconf_new.eq.0) then
        ngfmc=2*nstep*nblk
       else
        ngfmc=max(1,(nstep*nblk)/nconf_new)
      endif

c zero out estimators and averages
      if(irstar.ne.1) then
!JT        call my_second(1,'zerest')
        call zerest
      endif

c check if restart flag is on. If so then read input from
c dumped data to restart

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

c if there are equilibrium steps to take, do them here
c skip equilibrium steps if restart run
c imetro = 1 original force-bias
c        = 5 spherical-polar with exponential/linear T
c        = 6 spherical-polar with slater T
c        = 7 spherical-polar with slater T
      if(nblkeq.ge.1.and.irstar.ne.1) then
        l=0

        do 420 i=1,nblkeq
          do 410 j=1,nstep
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

c       Equilibration steps done. Zero out estimators again.
        call print_cpu_time_in_seconds ('End       of equilibration')
        call zerest
      endif

c now do averaging steps
      l=0
      do 440 i=1,1000000  !JT
        block_iterations_nb = block_iterations_nb + 1  !JT
        do 430 j=1,nstep
         step_iterations_nb = step_iterations_nb + 1   !JT
        l=l+1
        if(nloc.gt.0) call rotqua
        if(imetro.eq.1) then
          call metrop(l)
         elseif(imetro.eq.5) then
          call metrop_polar(l)
         elseif(imetro.eq.6) then
           call metrop_slat(l)
        endif

        call object_modified_by_index (xold_index)  !JT

!       accumulate data for averages and statitical errors
        call compute_averages        !JT
        call compute_averages_walk_step   !JT

c       write out configuration for optimization/dmc/gfmc here
        if (nconf_new /= 0) then
         if(mod(l,ngfmc).eq.1 .or. ngfmc.eq.1) then
          if(ndim*nelec.lt.100) then
           write(fmt,'(a1,i2,a21)')'(',ndim*nelec,'f14.8,i3,d12.4,f12.5)'
          else
           write(fmt,'(a1,i3,a21)')'(',ndim*nelec,'f14.8,i3,d12.4,f12.5)'
          endif
          write(7,fmt) ((xold(k,jj),k=1,ndim),jj=1,nelec),
     &    int(sign(1.d0,psido)),log(dabs(psido))+psijo,eold(1)
         endif
        endif

!       write out configurations in Scemama's format
        if (l_write_walkers) then
         if(mod(l,write_walkers_step) == 0) then
          do jj = 1, nelec
           write(file_walkers_out_unit,'(3f)') (xold(k,jj),k=1,ndim)
          enddo
           write(file_walkers_out_unit,*) dexp(psijo)*psido
         endif
        endif

  430   continue

!     compute averages and statitical errors
      call acuest
      call compute_averages_walk_block   !JT
      call compute_covariances  !JT
      call compute_variances  !JT
      call compute_errors  !JT

!     write at each block
      call objects_print_at_each_block
      call writing_routines !JT

!     exit loop if nblk and threshold on statistical error reached
      if (block_iterations_nb .ge. nblk .and. eloc_av_err .gt. 0 .and. eloc_av_err .le. error_threshold) then    !JT
        exit                              !JT
      endif                               !JT

  440 continue

!     reinitilization at the end of MC iterations
      block_nb = block_iterations_nb !JT save final block number
      step_iterations_nb  = 0   !JT
      block_iterations_nb = 0   !JT
      call reinit_averages_and_errors !JT
      call reinit_writing_routines !JT

      call print_cpu_time_in_seconds ('End       of accumulation ')

c print out final results
      call finwrt

c Presently grad_hess_jas_fin needs to be called after finwrt
c     if(ijasderiv.eq.1 .and. idtask.eq.0) then
c       passes=dfloat(iblk)*dfloat(nstep)*dfloat(nproc)
c       efin=ecum(1)/passes
c       call grad_hess_fin(passes,efin)
c     endif

c Write out last MC configuration on mc_configs_start
      call mc_configs_write

c if dump flag is on then dump out data for a restart
      if(idump.eq.1) then
        open(10,form='unformatted',file='restart_vmc')
        rewind 10
        call dumper
        close(unit=10)
      endif

      close(unit=9)
c     close(unit=5)
c     close(unit=6)
      if(nconf_new.ne.0) close(unit=7)

      return
      end

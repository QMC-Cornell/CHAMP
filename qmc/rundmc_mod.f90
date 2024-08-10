module rundmc_mod
    use constants_mod
    use walker_type_mod
    use psi_utils_mod
    use stats_type_mod
    use green_functions_mod
!    use mpi
    use mpi_mod, only: idtask
    use dim_mod, only: ndim
    use const_mod, only: nelec, ipr
    use contrldmc_mod, only: nfprod
    use contrl_mod, only: nconf_global, nconf
    implicit none

contains

    subroutine rundmc(nblock, nstep, tau, etrial, eest, eestfrag, undo_pop_ctrl, &
            walkers, nwalk, stats, l_eq)
        use iterat_mod, only: ipass
        use branch_mod, only: ff, MWALK
        use contrldmc_mod, only: nfprod, icross, tmoves, itau_integ, itau_eff, idmc, rtrttau, iacc_rej
        use split_join_mod, only: split_join
        use age_mod, only: ioldest, ioldestmx
        use stats_mod, only: nodecr, try_int, dr2ac, dr2un, acc, acc_int
        use mpi_mod, only: nproc
        use contr3_mod, only: mode
        use pseudo_mod, only: nloc
        use stats_index_mod
        use debug_shared_mod, only: iw, is
        use estcum_mod, only: iblk
        use control_mod, only: l_reset_etrial, limit_rewt_dmc, l_print_rewt
        use fragments_mod
        use basic_tools_mod
        use redistribute_mod
        implicit none
        integer, intent(in)                     :: nblock, nstep
        real(dp), intent(in)                    :: tau, etrial
        real(dp), intent(inout)                 :: undo_pop_ctrl, eest, eestfrag(:)
        type(walker_t), target, intent(inout)   :: walkers(:)
        integer, intent(inout)                  :: nwalk
        type(stats_t), intent(inout)            :: stats
        integer                 :: ie, ifrag, IERROR, npass_print, iacc(nelec)
!        integer                 :: iw, is, ie
        logical                 :: l_accept, l_eq
        type(walker_t), pointer :: walker
        real(dp)                :: xo(ndim, nelec), xn(ndim, nelec)
        real(dp)                :: pop_ctrl, pexp, ngen_pop_ctrl
        real(dp)                :: wtarg, wtot, wo, wf, wg, esigma, eint, smax
        real(dp)                :: rannyu, po, pn, psirat, pacc(nelec), ptot
        real(dp)                :: taunow, taueff, tauint
        real(dp)                :: ddiffusesq(nelec)
        real(dp)                :: rn(nelec) 
        real(dp)                :: bavg(stats%nest), bwgt(stats%nest)
        real(dp)                ::  avg(stats%nest),  err(stats%nest), err1(stats%nest)
        real(dp)                :: ddavgfrag(nfrag), ddtotfrag(nfrag)
        real(dp)                :: taunowfrag(nfrag+1), tauintfrag(nfrag+1), tauefffrag(nfrag+1)
        real(dp)                :: ees, eintn, einto, emin, enetau, etri, tint
        type(psi_t), target     :: psin, psio
        real(dp) :: eint_nofrag, enetau_nofrag, dwt_nofrag

        !Initialization of a few variables
        call psio%init
        call psin%init
        wtarg=nconf_global
        pexp=exp(-1d0/nfprod)
        ngen_pop_ctrl=ceiling(1d0/tau)
        try_int=0
        nodecr=0
        dr2ac=0
        dr2un=0
        acc=0
        acc_int=0
        iblk=0

        if (idtask.EQ.0) &
        write(6,'(t5,''egnow'', &
                & t15, ''egave'',  t21,  ''(egerr)'', &
                & t32, ''peave'',  t38,  ''(peerr)'', &
                & t49, ''tpbave'', t55,  ''(tpberr)'',&
                & t66, ''tjfave'', t72,  ''(tjferr)'',&
                & t83, ''fgave'',  t89,  ''(fgerr)'', &
                & t101,''npass'',  t111, ''wgsum'',t121,''ioldest'')')

!        undo_pop_ctrl=dexp((eest-etrial)*taueff/(1-pexp))
!        wtot=wtarg
        wtot=0
        do iw=1,nwalk
            wtot=wtot+walkers(iw)%weight
        enddo
        if (l_mode_dmc_mov1_mpi2) then
            call MPI_Allreduce(MPI_IN_PLACE,wtot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        endif
        taueff=tau
        if (l_fragments) tauefffrag=tau

        !Begin main loop over all Monte Carlo blocks and steps
        do is=1,nblock*nstep
            pop_ctrl=dexp((eest-etrial)*taueff)*(wtarg/wtot)**(1/ngen_pop_ctrl)
            if (ipr.LE.-8) then
                write(6,'(''pop_ctrl,1/eigv,eest,etrial,taueff,wtarg,wtot='',8f16.8)') &
                pop_ctrl,dexp((eest-etrial)*taueff),eest,etrial,taueff,wtarg,wtot
            endif
            if (l_fragments) then
                pop_ctrl=(wtarg/wtot)**(1/ngen_pop_ctrl)
                pop_ctrl=pop_ctrl*dexp(sum((eestfrag-etrialfrag)*tauefffrag))
            endif
!            undo_pop_ctrl=(undo_pop_ctrl)**pexp/pop_ctrl
            ipass=ipass+1
            ff(mod(ipass,nfprod))=1/pop_ctrl
            undo_pop_ctrl=product(ff)

            if (idmc.LE.0) then
                pop_ctrl=1
                undo_pop_ctrl=1
            endif

            !Rotate grid points for pseudopotential systems
            if (nloc.GT.0) call rotqua

            !Begin loop over all walkers in a generation
            wtot=0
            ioldest=-1
            do iw=1,nwalk
                walker=>walkers(iw)
                wo=walker%weight
                xo=walker%x
                !psio initialized to wave function at position walker%x
                psio=walker%psi
                walker%age=walker%age+1
                l_accept=.FALSE.

                !do_tmoves executes all Nelec T-moves and updates both xo and
                !psio. l_accept is true if any of the T-moves are accepted (if
                !l_accept is true at the end of the entire Monte Carlo step, we
                !calculate the wave function at the new position, otherwise we
                !skip that calculation to save time)
                if (tmoves) l_accept = do_tmoves(iw,xo,psio)

                !Begin loop over electrons
                ddiffusesq=0
                ptot=1
                iacc=0
                do ie=1,nelec
                    !The Assaraf, Moroni, Filippi method cannot efficiently
                    !calculate all wave function derivatives for all electrons
                    !after each one-electron move. Thus, we calculate the
                    !derivatives for only the "ie"th electron at the beginning of
                    !the "ie"th step of the electron loop using the subroutine
                    !"calc_1e_derivatives".
                    call calc_1e_derivatives(ie,xo,psio) 
                    !sample_1e_move samples from the stochastic part of the
                    !Green's function to propose a new electron configuration xn.
                    !xn is 3-Nelec dimensional and differs from xo only in the
                    !"ie"th electron position.
                    call sample_1e_move(tau,ie,psio,xo,ddiffusesq(ie),xn)
                    !move_1e uses Sherman-Morrison to calculate psin, the wave
                    !function at the position xn, from psio.
                    call move_1e(ie,xo,psio,xn,psin)
                    !sample_1e_prob calculates the proposal probability for the
                    !forward and backward moves.
                    po=sample_1e_prob(tau,ie,psio,xo,xn) !forward move
                    pn=sample_1e_prob(tau,ie,psin,xn,xo) !backward move
                    ptot=ptot*po
                    psirat=(psin%det/psio%det)**2*dexp(2*(psin%jas-psio%jas))
                    pacc(ie)=min(1d0, psirat*pn/po)

                    if (iacc_rej.EQ.0) pacc(ie)=1d0

                    !The following lines calculate statistics that are printed at
                    !the end of the run.
                    try_int=try_int+1
                    if (psin%det*psio%det.LE.0) then
                        if (icross.LE.0) pacc(ie)=0
                        nodecr=nodecr+1
                    endif
                    dr2ac=dr2ac+pacc(ie)*sum((xo(:,ie)-xn(:,ie))**2)
                    dr2un=dr2un+         sum((xo(:,ie)-xn(:,ie))**2)
                    acc=acc+pacc(ie)

                    if (ipr.LE.-9) then
                        write(6, '(''ie,pacc,psirat,po,pn= '',1i3,4f16.8)') &
                        ie,pacc(ie),psirat,po,pn
                        write(6, '(''xo,xn= '',6f16.8)') xo(:,ie),xn(:,ie)
                        write(6, '(''psido,psidn,psijo,psijn= '',4f16.8)') &
                        psio%det,psin%det,psio%jas,psin%jas
                    endif 

                    !The Metropolis-Hastings accept-reject step.
                    if (rannyu(0).LT.pacc(ie)) then
                        xo=xn
                        call copy_psi_1e(psio,psin,ie)
!                        psio=psin
                        acc_int=acc_int+1
                        iacc(ie)=1
                        l_accept=.TRUE.
                    endif
                enddo !do ie=1, nelec
                
                if (l_accept) then 
                    walker%age=0
                    call psi_at(xo,psio)
                endif

! Effective tau for branching
                if (itau_eff.GE.1) then
                    if (itau_eff.EQ.2) then
                        taunow=tau*sum(iacc*ddiffusesq)/sum(ddiffusesq)
                    else
                        taunow=tau*sum(pacc*ddiffusesq)/sum(ddiffusesq)
                    endif
                    !if (itau_eff.EQ.1) taunow=tau*sum(pacc*ddiffusesq)/sum(ddiffusesq)
                    !if (itau_eff.EQ.3) taunow=taueff
                    if (itau_eff.GT.3) call die('rundmc_mod', 'itau_eff cannot be greater than 3.')
                    if (l_fragments) then
                       ddavgfrag=0
                       ddtotfrag=0
                       do ie=1,nelec
                           ifrag=iwfragelec(ie)
                           if (itau_eff.EQ.2) then
                               ddavgfrag(ifrag)=ddavgfrag(ifrag)+iacc(ie)*ddiffusesq(ie)
                           else
                               ddavgfrag(ifrag)=ddavgfrag(ifrag)+pacc(ie)*ddiffusesq(ie)
                           endif
                           ddtotfrag(ifrag)=ddtotfrag(ifrag)+ddiffusesq(ie)
                       enddo
                       taunowfrag(1:nfrag)=tau*ddavgfrag/ddtotfrag
                       taunowfrag(nfrag+1)=taunow
                    endif
                else
                    taunow=tau
                    if (l_fragments) taunowfrag=tau
                endif
!                wo=wo*max(0d0, psio%det/walker%psi%det)*sho_green_function(etrial,taunow,walker%x,xo)/ptot !WARNING: no accrej

                if ((.NOT.l_eq).AND.(is.GT.max(10,nint(10.d0/tau)))) then
                    esigma=dsqrt(dble(wtarg))*stats%sigma(elocalg_)
                    !if (mode.eq.'dmc_mov1_mpi2'.OR.mode.eq.'dmc_mov1_mpi3') then
                    !    esigma=esigma*dsqrt(dble(nproc))
                    !endif
                else
                    esigma=0.2d0*dsqrt(dble(nelec))
                endif

                if (itau_integ.GE.1) then
                    if (itau_eff.EQ.3) then
                        tauint=taueff
                        if (l_fragments) tauintfrag=tauefffrag
                    else
                        tauint=taunow
                        if (l_fragments) tauintfrag=taunowfrag
                    endif
                else
                    tauint=tau
                    if (l_fragments) tauintfrag=tau
                endif

                !Reweighting step
                if (idmc.GT.0) then 
                    if (l_fragments) then
                        enetau=0
                        do ifrag=1,nfrag+1
                            etri=etrialfrag(ifrag)
                            ees =eestfrag(ifrag)
                            tint=merge(tauintfrag(ifrag),taunow,ifrag.LE.nfrag)
                            if ((.NOT.l_eq).AND.(is.GT.max(10,nint(10.d0/tau)))) then
                                esigma=dsqrt(dble(wtarg))*stats%sigma(offset_enefrag_+ifrag)
                                if (mode.eq.'dmc_mov1_mpi2'.OR.mode.eq.'dmc_mov1_mpi3') then
                                    esigma=esigma*dsqrt(dble(nproc))
                                endif
                            else
                                if (ifrag.LE.nfrag) then
                                  esigma=0.2d0*dsqrt(dble(nelecfrag(ifrag)))
                                else
                                  esigma=0.2d0*dsqrt(dble(nfrag))
                                endif
                            endif
                            eintn=ene_int2(etri,ees,esigma,tau,tint,      psio,ifrag)
                            einto=ene_int2(etri,ees,esigma,tau,tint,walker%psi,ifrag)
                            !emin =ees-etri-limit_rewt_dmc*esigma/dsqrt(tau)
                            emin =ees-etri-limit_rewt_dmc*esigma/rtrttau
                            eint =limit_ene_int((eintn+einto)/2,emin)
                            if (itau_eff.EQ.3) then
                                enetau=enetau+eint*tauefffrag(ifrag)
                            else
                                enetau=enetau+eint*taunowfrag(ifrag)
                            endif
                            if (ipr.LE.-8) then
                                write(6,'(''eint,taunow,etrial,eest,energy_sigma= '',6f16.8)') &
                                eint,taunowfrag(ifrag),etri,ees,esigma
                            endif
                        enddo
                        wo=wo*dexp(-enetau)
                        wo=wo*pop_ctrl
                    else
                        eint=ene_int_rewt(etrial,eest,esigma,tau,tauint,&
                            walker%psi,psio)
!                        eintn=ene_int_rewt(etrial,eest,esigma,tau,tauint,psio)
!                        eint=(eintn+walker%eint)/2
                        if (itau_eff.EQ.3) then
                            enetau=eint*taueff
                        else
                            enetau=eint*taunow
                        endif
                        !wo=wo*dexp(-eint*taunow)
                        wo=wo*dexp(-enetau)
                        wo=wo*pop_ctrl
                        if (ipr.LE.-8) then
                            write(6,'(''eint,taunow,etrial,eest,energy_sigma= '',6f16.8)') &
                            eint,taunow,etrial,eest,esigma
                        endif
                    endif
                endif

                if (ipr.LE.-8) then
                    write(6,'(''iw,wo,wg,enetau,pop_ctrl= '',1i8,6f16.8)') iw,wo,wo*undo_pop_ctrl,enetau,pop_ctrl
                    !write(0,'(''idtask,iw,wo,wg= '',2i8,2f16.8)') idtask,iw,wo,wo*undo_pop_ctrl
                endif

                walker%weight=wo
                walker%x=xo
                walker%psi=psio !TODO: this makes an extra copy. avoid this?
    
                rn=norm2(xo,DIM=1)
                wf=wo/pop_ctrl
                wg=wo*undo_pop_ctrl
                !Collect statistics for this step. In the line
                !call stats % accuest1(elocalg_, wg, psio%eloc)
                !elocalg_ is an INDEX, defined manually in commons/stats_index_mod.f90
                !wg is the WEIGHT used in the weighted average
                !psio%eloc: is the VALUE of the observable being averaged
                call stats % accuest1(elocalg_, wg, psio%eloc)
                call stats % accuest1(elocalf_, wf, psio%eloc)
                call stats % accuest1(elocal0_, wo, psio%eloc)
                call stats % accuest1(ekinpb_,  wg, psio%ekinpb)
                call stats % accuest1(ekinjf_,  wg, psio%ekinjf)
                call stats % accuest1(epot_,    wg, psio%epot)
                call stats % accuest1(epot_ee_, wg, psio%epot_ee)
                call stats % accuest1(taueff_,  wg, taunow)
                call stats % accuest1(R1_,      wg, sum(rn   )/nelec)
                call stats % accuest1(R2_,      wg, sum(rn**2)/nelec)
                call stats % accuest1(R3_,      wg, sum(rn**3)/nelec)
                call stats % accuest1(R4_,      wg, sum(rn**4)/nelec)
                call stats % accuest1(RI_,      wg, sum(1/rn )/nelec)
                call stats % accuest1(wtsq_,   1d0, wg**2)
                if (l_fragments) then
                    do ifrag=1,nfrag+1
                        call stats % accuest1(offset_enefrag_+ifrag, wg, psio%enefrag(ifrag))
                        call stats % accuest1(offset_taufrag_+ifrag, wg,   taunowfrag(ifrag))
                    enddo
                endif

                wtot=wtot+wo
                ioldest=max(ioldest, walker%age)
            enddo !iw=1, nwalk
            call stats % end_step !Finish collecting statistics for this step.
            if (l_mode_dmc_mov1_mpi2) then
                call MPI_Allreduce(MPI_IN_PLACE,ioldest,1,         MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
                call MPI_Allreduce(MPI_IN_PLACE,   wtot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
                call redistribute(walkers)
            endif
            !Calculate statistics for the total weight of this step.
            call stats % accuest(ws0_, 1d0, wtot)
            call stats % accuest(wsf_, 1d0, wtot/pop_ctrl)
            call stats % accuest(wsg_, 1d0, wtot*undo_pop_ctrl)
            call split_join(walkers, nwalk)
    
            if (mod(is,nstep).EQ.0) then !end block
                iblk=iblk+1
                call stats % end_block(bavg, bwgt, avg, err, err1) !Finish collecting statistics for this block.
!                if (idtask.EQ.0) then
                    if (l_mode_dmc_mov1_mpi2) then
                      npass_print=is
                    else
                      npass_print=is*nproc
                    endif
                    write(6,'(f10.5,4(f10.5,''('',i5,'')''),17x,3i10)') &
                        bavg(elocalg_),&
                        avg(elocalg_), nint(100000*err(elocalg_)),&
                        avg(epot_),    nint(100000*err(epot_)),&
                        avg(ekinpb_),  nint(100000*err(ekinpb_)),&
                        avg(ekinjf_),  nint(100000*err(ekinjf_)),&
                        npass_print, nint(bwgt(elocalg_)/nproc), ioldest
                    if (l_fragments) then
                        do ifrag=1,nfrag+1
                            write(6,'(''ifrag,ene,tau,ene_tcorr='',1I3,2(f10.5,''('',i5,'')''),1f8.2)') &
                                ifrag,avg(offset_enefrag_+ifrag), nint(100000*err(offset_enefrag_+ifrag)), &
                                      avg(offset_taufrag_+ifrag), nint(100000*err(offset_taufrag_+ifrag)), &
                                      (err(offset_enefrag_+ifrag)/max(err1(offset_enefrag_+ifrag),tiny(0d0)))**2
                        enddo
                    endif
!                endif
                if (.not.l_reset_etrial .and. iblk>1 .and. iblk<5) then
                  if(abs(avg(elocalg_)-etrial).gt.3*err(elocalg_)) then
                    write(6,'(''Warning: abs(eest-etrial).gt.3*egerr, &
                               &eest,etrial,egerr='',3es12.4,'' &
                               &Fix: use better etrial'')')  &
                        avg(elocalg_),etrial,err(elocalg_)
                  endif
                elseif(.not.l_reset_etrial .and. iblk>1 .and. iblk==5) then
                  if(abs(avg(elocalg_)-etrial).gt.5*err(elocalg_)) then
                    write(6,'(''Warning: abs(eest-etrial).gt.5*egerr, &
                               &eest,etrial,egerr='',3es12.4,'' &
                               &Fix: use better etrial'')') &
                        avg(elocalg_),etrial,err(elocalg_)
                    call die ('rundmc', &
                             'Warning: abs(eest-etrial).gt.5*egerr. &
                             &Fix: use better etrial')
                  endif
                endif
            endif

            ioldestmx=max(ioldestmx, ioldest)
            taueff=stats % avg(taueff_)
            eest  =stats % avg(elocalg_)
            if (l_fragments) then
                do ifrag=1,nfrag+1
                    eestfrag(ifrag)   = stats % avg(offset_enefrag_+ifrag)
                    tauefffrag(ifrag) = stats % avg(offset_taufrag_+ifrag)
                enddo
            endif
        enddo !is=1, nblock*nstep
    end subroutine rundmc
end module rundmc_mod

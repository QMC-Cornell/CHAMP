module psi_utils_mod
    use constants_mod, only: dp
    use const_mod, only: nelec
    use dim_mod, only: ndim
    use atom_mod, only: ncent
    use contrldmc_mod, only: l_psi_approx, l_tau_diffusion, l_modified_adrift
    use psi_type_mod
    implicit none

contains
!------------------------------------------------------------------------------!

    subroutine psi_at(x, psi)
        use contrl_opt_mod, only: nparm
        use const_mod, only: nelec
        use pseudo_mod, only: nloc 
        use psi_dummy_mod, only: psi_dummy
        use fragments_mod, only: iwfragelec, enefrag, l_fragments
        implicit none
        real(dp), intent(in)          :: x(:,:)
        type(psi_t), intent(inout)    :: psi
        real(dp)    :: d2lnpsi, pei
        real(dp)    :: divv(size(x,2)), denergy(nparm)
        integer :: ie, k

        if (nloc.GT.0) call set_newvar_pointers(1,psi_dummy) !nonloc uses orbe to calculate nonlocal contribution to pe
        call set_oldvar_pointers(psi)
        call hpsi(x,psi%det,psi%jas,psi%grad,divv,d2lnpsi,psi%epot,psi%epot_ee,psi%eloc,denergy,1)
!        call copy_globalo_to_psi(psi)
        call nullify_pointers
        if (l_tau_diffusion.OR.l_psi_approx.OR.l_modified_adrift) then
            do ie=1,nelec
                call calc_approx_hess(ie,x,psi)
            enddo
        endif
        psi%ekinpb=psi%eloc-psi%epot
        psi%ekinjf=psi%ekinpb/2+sum(psi%grad**2)/4
        if (l_fragments) then
            psi%iwfragelec=iwfragelec !needed to calculate einto in rundmc
            psi%enefrag=enefrag !ditto  
        endif
    end subroutine psi_at
!------------------------------------------------------------------------------!

    subroutine move_1e(ie,xo,psio,xn,psin)
        use contrl_opt_mod, only: nparm
        use dete_mod, only: tupn, tdnn
        use dets_mod, only: nup  
        use gamma_mod
        implicit none
        integer, intent(in)           :: ie
        type(psi_t), intent(inout), target    :: psio
        type(psi_t), intent(inout), target    :: psin
        real(dp), intent(in)          :: xo(:,:), xn(:,:)
        integer :: i, j

        if (l_tau_diffusion.OR.l_psi_approx.OR.l_modified_adrift) then
!Without the following somewhat expensive copy, one must be more careful about 
!how the pointers are handled within the later call to calc_approx_hess.
!Because calc_approx_hess recalculates some arrays that depend on the
!positions of both the 'up' spin and 'down' spin electrons, e.g. the 'y'
!matrix, some pointers must point to arrays in psin, but others must point
!to arrays in psio. More specifically, if an up spin electron is moved, then
!the pointer 'invdn' should point to psio%invdn, since after the call to
!hpsie the contents of psin%invdn are totally meaningless. On the other hand,
!'invup' should point to psin%invup, and 'invupn' should point to
!psi_dummy%invup.
!The following is a simple but expensive way to get around this problem.
!Since the 2nd derivatives are calculated numerically and thus the program is
!already rather slow, this copy is proportionately less impactful to overall
!performance than if tau_diffusion or psi_approx were not used.
            psin=psio
        else
            psin%r_en    = psio%r_en
            psin%r_ee    = psio%r_ee
            psin%rvec_en = psio%rvec_en
            psin%rvec_ee = psio%rvec_ee
        endif

        call set_oldvar_pointers(psio)
        call set_newvar_pointers(ie,psin) !WARNING: must use set_newvar_pointers after set_oldvar_pointers
        call hpsie(ie,xn,psin%det,psin%jas,psin%grad,psin%lapl)
        if (ie.LE.nup) then           
            psin%tup(:,occup)=tupn
        else         
            psin%tdn(:,occdn)=tdnn
        endif  
        call nullify_pointers
        if (l_tau_diffusion.OR.l_psi_approx.OR.l_modified_adrift) call calc_approx_hess(ie,xn,psin)
    end subroutine move_1e
!------------------------------------------------------------------------------!

    subroutine set_oldvar_pointers(psi)
        use jaso_mod, only: lapjo, lapjijo, fsumo, fso, fjo, fijo, d2o, d2ijo
        use orb_mod, only: orb, dorb, ddorb
        use deriv_fast_mod, only: aiup, deta_up, tup, detup, invup,&
                                  aidn, deta_dn, tdn, detdn, invdn,&
                                  chi, yup, ydn
        use distance_mod, only: rvec_en, rvec_ee, r_en, r_ee
        use qua_mod, only: quadr, quadx 
        implicit none
        type(psi_t), intent(inout), target :: psi

        !jas vars
        lapjo   => psi%lapj
        lapjijo => psi%lapjij
        fsumo   => psi%fsum
        fso     => psi%fs  
        fjo     => psi%fj  
        fijo    => psi%fij 
        d2o     => psi%d2  
        d2ijo   => psi%d2ij

        !det vars
        orb     => psi%orb 
        dorb    => psi%dorb
        ddorb   => psi%ddorb

        aiup    => psi%aiup   
        deta_up => psi%deta_up
        tup     => psi%tup    
        detup   => psi%detup  
        invup   => psi%invup  

        aidn    => psi%aidn   
        deta_dn => psi%deta_dn
        tdn     => psi%tdn    
        detdn   => psi%detdn  
        invdn   => psi%invdn  

        chi     => psi%chi
        yup     => psi%yup
        ydn     => psi%ydn

        !other vars
        r_en    => psi%r_en
        r_ee    => psi%r_ee
        rvec_en => psi%rvec_en
        rvec_ee => psi%rvec_ee
        quadr   => psi%quadr
        quadx   => psi%quadx
    end subroutine set_oldvar_pointers
!------------------------------------------------------------------------------!

    subroutine set_newvar_pointers(ie,psi)
        use jasn_mod, only: lapjn, lapjijn, fsumn, fsn, fjn, fijn, d2n, d2ijn
        use orbe_mod, only: orbe, dorbe, ddorbe
        use dete_mod, only: aiupn, deta_upn, tupn, detupn, invupn,&
                            aidnn, deta_dnn, tdnn, detdnn, invdnn,&
                            chin, yupn, ydnn
        use deriv_fast_mod, only: aiup, deta_up, tup, detup, invup,&
                                  aidn, deta_dn, tdn, detdn, invdn
        use dets_mod, only: nup, ndn
        use distance_mod, only: rvec_en, rvec_ee, r_en, r_ee
        use gamma_mod
        implicit none
        integer, intent(in) :: ie
        type(psi_t), intent(inout), target :: psi

        lapjn   => psi%lapj
        lapjijn => psi%lapjij
        fsumn   => psi%fsum
        fsn     => psi%fs  
        fjn     => psi%fj  
        fijn    => psi%fij 
        d2n     => psi%d2  
        d2ijn   => psi%d2ij

        orbe    => psi%orb(ie,:)
        dorbe   => psi%dorb(:,ie,:)
        ddorbe  => psi%ddorb(ie,:)

        aiupn    => psi%aiup
        deta_upn => psi%deta_up
!        tupn     => psi%tup(:,occup)
!        call set_2d_pointer(tupn, psi%tup(:,occup))
!        tupn     => psi%tup
        if (.NOT.associated(tupn)) allocate(tupn(nup,noccup))
        tupn     = psi%tup(:,occup)
        detupn   => psi%detup
        invupn   => psi%invup

        aidnn    => psi%aidn   
        deta_dnn => psi%deta_dn
!        tdnn     => psi%tdn(:,occdn)
!        call set_2d_pointer(tdnn, psi%tdn(:,occdn))
!        tdnn     => psi%tdn
        if (.NOT.associated(tdnn)) allocate(tdnn(ndn,noccdn))
        tdnn     = psi%tdn(:,occdn)
        detdnn   => psi%detdn  
        invdnn   => psi%invdn  

        chin => psi%chi
        yupn => psi%yup
        ydnn => psi%ydn

        !other
        r_en    => psi%r_en
        r_ee    => psi%r_ee
        rvec_en => psi%rvec_en
        rvec_ee => psi%rvec_ee
    end subroutine set_newvar_pointers
!------------------------------------------------------------------------------!

    subroutine set_2d_pointer(ptr, mat)
        implicit none
        real(dp), pointer :: ptr(:,:)
        real(dp), target :: mat(:,:)

        ptr => mat
    end subroutine

!------------------------------------------------------------------------------!
    subroutine nullify_pointers
        use jaso_mod, only: lapjo, lapjijo, fsumo, fso, fjo, fijo, d2o, d2ijo
        use jasn_mod, only: lapjn, lapjijn, fsumn, fsn, fjn, fijn, d2n, d2ijn
        use orb_mod, only: orb, dorb, ddorb
        use orbe_mod, only: orbe, dorbe, ddorbe
        use deriv_fast_mod, only: aiup, deta_up, tup, detup, invup,&
                                  aidn, deta_dn, tdn, detdn, invdn,&
                                  chi, yup, ydn
        use dete_mod, only: aiupn, deta_upn, tupn, detupn, invupn,&
                            aidnn, deta_dnn, tdnn, detdnn, invdnn,&
                            chin, yupn, ydnn
        use distance_mod, only: rvec_en, rvec_ee, r_en, r_ee
        use qua_mod, only: quadr, quadx 
        implicit none

        !jas vars
        nullify(lapjo  )
        nullify(lapjijo)
        nullify(fsumo  )
        nullify(fso    )
        nullify(fjo    )
        nullify(fijo   )
        nullify(d2o    )
        nullify(d2ijo  )

        nullify(lapjn  )
        nullify(lapjijn)
        nullify(fsumn  )
        nullify(fsn    )
        nullify(fjn    )
        nullify(fijn   )
        nullify(d2n    )
        nullify(d2ijn  )

        !det vars
        nullify(orb    )
        nullify(dorb   )
        nullify(ddorb  )

        nullify(orbe   )
        nullify(dorbe  )
        nullify(ddorbe )

        nullify(aiup   )
        nullify(deta_up)
        nullify(tup    )
        nullify(detup  )
        nullify(invup  )

        nullify(aidn   )
        nullify(deta_dn)
        nullify(tdn    )
        nullify(detdn  )
        nullify(invdn  )

        nullify(chi    )
        nullify(yup    )
        nullify(ydn    )

        nullify(aiupn  )
        nullify(deta_upn)
!        nullify(tupn   )
        nullify(detupn )
        nullify(invupn )

        nullify(aidnn  )
        nullify(deta_dnn)
!        nullify(tdnn   )
        nullify(detdnn )
        nullify(invdnn )

        nullify(chin   )
        nullify(yupn   )
        nullify(ydnn   )

        !other vars
        nullify(r_en   )
        nullify(r_ee   )
        nullify(rvec_en)
        nullify(rvec_ee)
        nullify(quadr  )
        nullify(quadx  )
    end subroutine nullify_pointers
!------------------------------------------------------------------------------!

    subroutine copy_psi_1e(psio,psin,ie)
        use dets_mod, only: nup
        use dorb_mod, only: ndetup, ndetdn
        use gamma_mod
        use deriv_fast_mod, only: orddn, ordup
        use dete_mod, only: tdnn, tupn
        implicit none
        type(psi_t), intent(inout) :: psio
        type(psi_t), intent(in) :: psin
        integer, intent(in) :: ie
        integer :: i, j, idet, ord

!        psio=psin

        psio%det = psin%det
        psio%jas = psin%jas

        psio%d2  = psin%d2
        psio%fsum= psin%fsum
        do i=1,nelec
          psio%lapj(i)=psin%lapj(i)
          psio%fj(1,i)=psin%fj(1,i)
          psio%fj(2,i)=psin%fj(2,i)
          psio%fj(3,i)=psin%fj(3,i)
        enddo

        do j=1,ie
          psio%fs(ie,j)=psin%fs(ie,j)
        enddo

        do j=ie+1,nelec
          psio%fs(j,ie)=psin%fs(j,ie)
        enddo

        do j=1,nelec
          psio%lapjij(ie,j)=psin%lapjij(ie,j)
          psio%lapjij(j,ie)=psin%lapjij(j,ie)
          psio%fij(1,ie,j) =psin%fij(1,ie,j)
          psio%fij(2,ie,j) =psin%fij(2,ie,j)
          psio%fij(3,ie,j) =psin%fij(3,ie,j)
          psio%fij(1,j,ie) =psin%fij(1,j,ie)
          psio%fij(2,j,ie) =psin%fij(2,j,ie)
          psio%fij(3,j,ie) =psin%fij(3,j,ie)
        enddo

        psio%orb(ie,:)    = psin%orb(ie,:)
        psio%dorb(:,ie,:) = psin%dorb(:,ie,:)
        psio%ddorb(ie,:)  = psin%ddorb(ie,:)

        if (ie.le.nup) then
          psio%aiup         = psin%aiup
          psio%deta_up      = psin%deta_up
!          psio%tup(:,occup) = psin%tup(:,occup)
          psio%tup(:,occup) = tupn
          psio%detup        = psin%detup
          do idet=1,ndetup
            ord = ordup(idet)
            do i=1,ord
              do j=1,ord
                psio%invup(j+ord*(i-1),idet) = psin%invup(j+ord*(i-1),idet)
              enddo
            enddo
          enddo
        else
          psio%aidn         = psin%aidn
          psio%deta_dn      = psin%deta_dn
!          psio%tdn(:,occdn) = psin%tdn(:,occdn)
          psio%tdn(:,occdn) = tdnn
          psio%detdn        = psin%detdn
          do idet=1,ndetdn
            ord = orddn(idet)
            do i=1,ord
              do j=1,ord
                psio%invdn(j+ord*(i-1),idet) = psin%invdn(j+ord*(i-1),idet)
              enddo
            enddo
          enddo
        endif
        psio%chi = psin%chi
        psio%yup = psin%yup
        psio%ydn = psin%ydn

        !other
        !TODO: don't need to copy all of these, right?
        psio%r_en    = psin%r_en
        psio%r_ee    = psin%r_ee
        psio%rvec_en = psin%rvec_en
        psio%rvec_ee = psin%rvec_ee
    end subroutine copy_psi_1e
!------------------------------------------------------------------------------!

    logical function do_tmoves(iw,xo,psio)
        use qua_mod, only: iaccept_tmove
        use config_dmc_mod, only: xoldw, voldw, psidow, psijow
        use psi_dummy_mod, only: psi_dummy
        implicit none
        integer, intent(in)         :: iw
        real(dp), intent(inout)     :: xo(:,:)
        type(psi_t), intent(inout)  :: psio

        xoldw(:,:,iw,1)=xo
        voldw(:,:,iw,1)=psio%grad
        psidow(iw,1)=psio%det
        psijow(iw,1)=psio%jas
        call set_newvar_pointers(1,psi_dummy)
        call set_oldvar_pointers(psio)
!        call copy_psi_to_globalo(psio)
        call getvps(psio)
        call tmove
!        call copy_globalo_to_psi(psio)
        call nullify_pointers
        xo=xoldw(:,:,iw,1)
        psio%grad=voldw(:,:,iw,1)
        psio%det=psidow(iw,1)
        psio%jas=psijow(iw,1)
        do_tmoves=(iaccept_tmove.NE.0)
    end function do_tmoves
!------------------------------------------------------------------------------!

    subroutine calc_approx_hess(ie,x,psi)
        use psi_dummy_mod, only: psi_dummy
        implicit none
        integer, intent(in)         :: ie
        real(dp), intent(in)        :: x(:,:)                                
        type(psi_t), intent(inout)  :: psi
        real(dp)    :: d2psi
        integer     :: k

        call set_oldvar_pointers(psi)
        call set_newvar_pointers(ie,psi_dummy)
!        call copy_psi_to_globalo(psi)
        call d2psi_numericale(1d-3,ie,x,psi%det,psi%jas,psi%grad(:,ie),d2psi)
        call nullify_pointers
        psi%hess(:,:,ie)=0d0
        do k=1,size(x,1)
            psi%hess(k,k,ie)=d2psi
        enddo
    end subroutine calc_approx_hess
!------------------------------------------------------------------------------!

    real(dp) function numerical_lapl(ie,x,psi) result(lapl)
        use psi_dummy_mod, only: psi_dummy
        implicit none
        integer, intent(in)         :: ie
        real(dp), intent(in)        :: x(:,:)
        type(psi_t), intent(inout)  :: psi
        real(dp)                    :: laplx, laply, laplz

        call set_oldvar_pointers(psi)
        call set_newvar_pointers(ie,psi_dummy)
        call d2psi_numericale(1d-3,ie,x,psi%det,psi%jas,[1d0,0d0,0d0],laplx)
        call d2psi_numericale(1d-3,ie,x,psi%det,psi%jas,[0d0,1d0,0d0],laply)
        call d2psi_numericale(1d-3,ie,x,psi%det,psi%jas,[0d0,0d0,1d0],laplz)
        call nullify_pointers
        lapl=laplx+laply+laplz
    end function numerical_lapl
!------------------------------------------------------------------------------!

    subroutine calc_1e_derivatives(ie,x,psi)
        use dets_mod, only: nup
        use dete_mod, only: eval_grad, eval_lapl
        implicit none
        integer, intent(in)         :: ie
        real(dp), intent(in)        :: x(:,:)
        type(psi_t), intent(inout)  :: psi

        if (ie.le.nup) then 
            call eval_grad(ie, psi%dorb(:,ie,:),psi%aiup,psi%tup,psi%yup,&
                           psi%grad(:,ie))
            if (l_tau_diffusion.OR.l_psi_approx) then
                call eval_lapl(ie,psi%ddorb(  ie,:),psi%aiup,psi%tup,psi%yup,&
                               psi%lapl)
            endif
        else
            call eval_grad(ie, psi%dorb(:,ie,:),psi%aidn,psi%tdn,psi%ydn,&
                           psi%grad(:,ie))
            if (l_tau_diffusion.OR.l_psi_approx) then
                call eval_lapl(ie,psi%ddorb(  ie,:),psi%aidn,psi%tdn,psi%ydn,&
                               psi%lapl)
            endif
        endif
        psi%lapl=psi%lapl+2*sum(psi%fj(:,ie)*psi%grad(:,ie))+psi%lapj(ie)
        psi%grad(:,ie)=psi%grad(:,ie)+psi%fj(:,ie)
        if (l_tau_diffusion.OR.l_psi_approx.OR.l_modified_adrift) call calc_approx_hess(ie,x,psi)
    end subroutine calc_1e_derivatives
!------------------------------------------------------------------------------!


    subroutine getvps(psi)
        use pseudo_mod, only: nloc
        implicit none
        type(psi_t), intent(in) :: psi
        integer :: ie

        !calculate psp
            do ie=1,nelec
              if(nloc.eq.1) then
                call getvps_fahy(psi%r_en,ie)
               elseif(nloc.ge.2 .and. nloc.le.5) then
                call getvps_champ(psi%r_en,ie)
               elseif(nloc.eq.6) then
                call getvps_gauss(psi%r_en,ie)
               else
                stop 'nloc < 1 or 6 < nloc'
              endif
            enddo
    end subroutine getvps
!------------------------------------------------------------------------------!

    subroutine d2psi_numericale(h, ie, x, psidn2, psijn2, dir, d2psi)
        use dim_mod, only: ndim
        use atom_mod, only: ncent
        use const_mod, only: nelec
        use distance_mod, only: rvec_ee, rvec_en, r_ee, r_en
        implicit none
    
        real(dp), intent(in)    :: h, x(:,:), dir(3), psidn2, psijn2
        integer, intent(in)     :: ie
        real(dp), intent(out)   :: d2psi
        real(dp) ::   lapl
        real(dp) ::      d(size(x,1), size(x,2))
        real(dp) ::   vnew(size(x,1), size(x,2))
        real(dp) ::   xnow(size(x,1), size(x,2))
    
!        real(dp) :: divn, lapn
        real(dp) :: psidn1, psidn3
        real(dp) :: psijn1, psijn3

        real(dp) :: rvec_ee_save(ndim,nelec*(nelec-1)/2)
        real(dp) :: rvec_en_save(ndim,nelec,ncent)
        real(dp) :: r_ee_save(nelec*(nelec-1)/2)
        real(dp) :: r_en_save(nelec,ncent)

        real(dp) :: a, b, c

        rvec_ee_save=rvec_ee
        rvec_en_save=rvec_en
        r_ee_save=r_ee
        r_en_save=r_en
    
        d = 0
        d(:,ie) = dir/norm2(dir)

        xnow=x+h*d
        call hpsie(ie,xnow,psidn3,psijn3,vnew,lapl)
        xnow=x-h*d
        call hpsie(ie,xnow,psidn1,psijn1,vnew,lapl)
!        xnow=x
!        call hpsie(ie,xnow,psidn2,psijn2,vnew,lapl)
    
        a=psidn1*dexp(psijn1-psijn2)/psidn2
        b=psidn3*dexp(psijn3-psijn2)/psidn2
        c=a+b
        d2psi=c-2
!        d2psi = psidn1*dexp(psijn1-psijn2)/psidn2 - 2 + psidn3*dexp(psijn3-psijn2)/psidn2
        d2psi = d2psi/(h**2)

        rvec_ee=rvec_ee_save
        rvec_en=rvec_en_save
        r_ee=r_ee_save
        r_en=r_en_save
    
        if (isnan(d2psi)) then
          write(6,'(''warning: d2psi is nan: d2psi, psidn, psijn'',99e20.12)') d2psi, psidn2, psijn2
          write(0,'(''warning: d2psi is nan: d2psi, psidn, psijn'',99e20.12)') d2psi, psidn2, psijn2
        endif
    end subroutine
end module psi_utils_mod

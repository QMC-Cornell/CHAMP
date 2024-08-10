  subroutine determinante(iel,x,rvec_en,r_en,grad,lap,psi)
    use dete_mod
    use constants_mod
    use control_mod
    use dorb_mod
    use slatn_mod
    use orbe_mod
    use dets_mod
    use slater_mod
    use const_mod
    use contr2_mod
    use wfsec_mod
    use contrl_per_mod
    use contrldmc_mod, only: l_tau_diffusion, l_psi_approx
    implicit none

    integer iel, i
    real(dp) x(3,*),rvec_en(3,nelec,*),r_en(nelec,*)
    real(dp) grad(3,nelec),lap,div,psi

!get orbital values and derivatives for electron iel
    if(iperiodic.eq.0 .or. iperiodic.eq.1) then

      if(inum_orb.eq.0) then
        call orbitals_loc_ana_grade(iel,rvec_en,r_en,orbe,dorbe,ddorbe)
      else
        call orbitals_loc_num_grade(iel,x,orbe,dorbe,ddorbe)
      endif

    else

      if(inum_orb.eq.0) then
        call orbitals_pw_grade(x(1,iel),orbe,dorbe,ddorbe)
      else
        call orbitals_period_num_grade(x(1,iel),orbe,dorbe,ddorbe)
      endif

    endif

!update psid
    call dete_update(iel)
    if (iel.le.nup) then
      psi = chin*deta_upn*deta_dn
      call eval_grad(iel,  dorbe, aiupn, tupn, yupn, grad(:,iel))
      if (l_tau_diffusion.OR.l_psi_approx) then
        call eval_lapl(iel, ddorbe, aiupn, tupn, yupn, lap)
      endif
    else
      psi = chin*deta_up*deta_dnn
      call eval_grad(iel,  dorbe, aidnn, tdnn, ydnn, grad(:,iel))
      if (l_tau_diffusion.OR.l_psi_approx) then
        call eval_lapl(iel, ddorbe, aidnn, tdnn, ydnn, lap)
      endif
    endif
  end subroutine

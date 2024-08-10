module green_functions_mod
    use constants_mod
    use psi_utils_mod
    use pseudo_mod, only: nloc
    implicit none

    real(dp), parameter :: pi       = 3.141592653589793d0
    real(dp), parameter :: eps      = 1d-10
    real(dp), parameter :: adrift0  = 1d-1
contains
!------------------------------------------------------------------------------!

    subroutine sample_1e_move(tau,ie,psio,xo,ddiffusesq,xn)
        use contrldmc_mod, only: l_tau_diffusion, l_psi_approx
        implicit none
        real(dp), intent(in)     :: tau, xo(:,:)
        integer, intent(in)     :: ie
        type(psi_t), intent(in) :: psio
        real(dp), intent(out)    :: ddiffusesq, xn(:,:)

        xn=xo
        if (nloc.EQ.0) then !all-electron calculation
            call sample_1e_move_drift_diffusion_with_cusp(&
               tau,psio%grad(:,ie),psio%hess(:,:,ie),xo(:,ie),&
               psio%rvec_en(:,ie,:),psio%r_en(ie,:),ddiffusesq,xn(:,ie))
        else if (l_tau_diffusion) then
            call sample_1e_move_taudif(&
                tau,psio%grad(:,ie),psio%hess(:,:,ie),psio%lapl,xo(:,ie),&
                ddiffusesq,xn(:,ie))
        else if (l_psi_approx) then
            call sample_1e_move_psia(&
                tau,psio%grad(:,ie),psio%hess(:,:,ie),psio%lapl,xo(:,ie),&
                ddiffusesq,xn(:,ie))
        else
            call sample_1e_move_drift_diffusion(&
                tau,psio%grad(:,ie),psio%hess(:,:,ie),xo(:,ie),&
                ddiffusesq,xn(:,ie))
        endif
    end subroutine sample_1e_move
!------------------------------------------------------------------------------!

    real(dp) function sample_1e_prob(tau,ie,psio,xo,xn)
        use contrldmc_mod, only: l_tau_diffusion, l_psi_approx
        implicit none
        real(dp), intent(in)     :: tau, xo(:,:), xn(:,:)
        integer, intent(in)     :: ie
        type(psi_t), intent(in) :: psio

        if (nloc.EQ.0) then
            sample_1e_prob = sample_1e_prob_drift_diffusion_with_cusp(&
                tau,psio%grad(:,ie),psio%hess(:,:,ie),xo(:,ie),&
                psio%rvec_en(:,ie,:),psio%r_en(ie,:),xn(:,ie))
        else if (l_tau_diffusion) then
            sample_1e_prob = sample_1e_prob_taudif(&
                tau,psio%grad(:,ie),psio%hess(:,:,ie),psio%lapl,&
                xo(:,ie),xn(:,ie))
        else if (l_psi_approx) then
            sample_1e_prob = sample_1e_prob_psia(&
                tau,psio%grad(:,ie),psio%hess(:,:,ie),psio%lapl,&
                xo(:,ie),xn(:,ie))
        else
            sample_1e_prob = sample_1e_prob_drift_diffusion(&
                tau,psio%grad(:,ie),psio%hess(:,:,ie),xo(:,ie),&
                xn(:,ie))
        endif
    end function sample_1e_prob
!------------------------------------------------------------------------------!

    subroutine sample_1e_move_drift_diffusion(tau,grado,hesso,xo,ddiffusesq,xn)
        implicit none
        real(dp), intent(in)     :: tau, grado(:), hesso(:,:), xo(:)
        real(dp), intent(out)    :: ddiffusesq, xn(:)
        real(dp)                 :: vav(size(xo)), xdiffuse(size(xo)), gauss
        integer                 :: k

        vav=averaged_1e_velocity(tau,grado,hesso)
        xdiffuse=dsqrt(tau)*[ (gauss(), k=1,size(xo)) ]
        xn=xo+vav*tau+xdiffuse
        ddiffusesq=sum(xdiffuse**2)
    end subroutine sample_1e_move_drift_diffusion
!------------------------------------------------------------------------------!

    subroutine sample_1e_move_drift_diffusion_with_cusp(tau,grado,hesso,xo,xen,ren,ddiffusesq,xn)
        use atom_mod, only: ncent, cent, znuc, iwctype
        implicit none
        real(dp), intent(in)     :: tau, grado(:), hesso(:,:), xo(:), xen(:,:), ren(:)
        real(dp), intent(out)    :: ddiffusesq, xn(:)
        real(dp)                 :: vav(size(xo)), xdiffuse(size(xo)), gauss
        real(dp)                :: z(size(xo)), x(size(xo)), zi, zf, pgaus, xf, xdrift(size(xo))
        real(dp)                :: rannyu
        real(dp)                :: r, costht, sintht, phi, zeta
        integer                 :: k, iwnuc
!Find index of nucleus closest to electron being moved
        iwnuc=minloc(ren, DIM=1)
!Place z-axis along direction from nearest nucleus to electron and
!x-axis along direction of angular component of velocity.
        z=xen(:,iwnuc)
        z=z/norm2(z)
        x=grado-sum(grado*z)*z
        if (norm2(x).LT.eps) then
            x(1)=eps*(1-z(1)**2)
            x(2)=eps*(-z(1)*z(2))
            x(3)=eps*(-z(1)*z(3))
        endif
        x=x/norm2(x)
!Calculate averaged velocity  
        vav=averaged_1e_velocity(tau,grado,hesso,xen)
!Prob. of sampling exponential rather than gaussian is
!half*derfc(zf/dsqrt(two*tau)) = half*(two-derfc(-zf/dsqrt(two*tau)))
!We use both expressions because under AIX the former is faster if zf>0
!and the latter is faster if zf<0.
!Note that if adrift is always set to 1 then it may be better to use
!vavvt rather than tau since the max drift distance is dsqrt(2*tau/adrift),
!so both the position of the gaussian and its width are prop to dsqrt(tau)
!if tau is used in derfc, and so qgaus does not tend to 1 for large tau.
!However, we are using a variable adrift that can be very small and then
!using tau is the better choice.
        zi=ren(iwnuc)
        zf=zi+sum(vav*z)*tau
        if (zf.GE.0) then
            pgaus=1-derfc(zf/dsqrt(2*tau))/2
        else
            pgaus=derfc(-zf/dsqrt(2*tau))/2
        endif
!The electron cannot drift past the nucleus
        zf=max(0d0,zf)
!The drift bends toward the nucleus
        xf=sum(vav*x)*tau*min(1d0,zf/((zf+zi)/2))
        xdrift=cent(:,iwnuc)+xf*x+zf*z
!Sample from drift-diffusion Green's function with probability pgaus.
!With probability 1-pgaus, sample from exponential centered at nearest nucleus.
        zeta=dsqrt(1/tau+znuc(iwctype(iwnuc))**2)
        if (rannyu(0).lt.pgaus) then
            xdiffuse=dsqrt(tau)*[ (gauss(), k=1,size(xo)) ]
            xn=xdrift+xdiffuse
        else
            r=(-half/zeta)*dlog(rannyu(0)*rannyu(0)*rannyu(0))
            costht=two*(rannyu(0)-half)
            sintht=sqrt(one-costht*costht)
            phi=two*pi*rannyu(0)
            xdiffuse = [r*sintht*cos(phi), r*sintht*sin(phi), r*costht]
            xn=cent(:,iwnuc)+xdiffuse
        endif
        ddiffusesq=sum((xn-xdrift)**2)
    end subroutine sample_1e_move_drift_diffusion_with_cusp
!------------------------------------------------------------------------------!

    subroutine sample_1e_move_taudif(tau,grado,hesso,laplo,xo,ddiffusesq,xn)
        use contrldmc_mod, only: l_method_of_images
        implicit none
        real(dp), intent(in)     :: tau, grado(:), hesso(:,:), laplo, xo(:)
        real(dp), intent(out)    :: ddiffusesq, xn(:)
        real(dp)    :: vnorm, vhat(size(xo)), vav(size(xo))
        real(dp)    :: d2psi, divpll, divprp
        real(dp)    :: taupll, tauprp
        real(dp)    :: dnode, xdif(size(xo)), xpll(size(xo)), ximg(size(xo))
        real(dp)    :: gauss, rannyu, pacc
        integer     :: k, cnt

        vnorm=norm2(grado)
        vhat=grado/vnorm
        d2psi=sum(vhat*matmul(hesso,vhat))
        divpll=d2psi-sum(grado**2)
        divprp=(laplo-d2psi)/(size(xo)-1)
        taupll=tau_diffusion(tau,divpll)
        tauprp=tau_diffusion(tau,divprp)
        if (l_method_of_images) then
            if (d2psi.LE.0) then
                dnode=abs((-vnorm+dsqrt(vnorm**2-2*d2psi))/d2psi)
            else
                dnode=1d3
            endif
        endif
        cnt=0
        do while (.TRUE.)
            xdif=[(gauss(),k=1,size(xo))]
            xpll=sum(xdif*vhat)*vhat
            xdif=dsqrt(taupll)*xpll+dsqrt(tauprp)*(xdif-xpll)
            vav=averaged_1e_velocity(tau,grado,hesso)
            xn=xo+vav*tau+xdif
            if (.NOT.l_method_of_images) exit

            ximg=xo-2*dnode*vhat-vav*tau
            pacc=1-dexp((sum(xdif*vhat)**2-sum((xn-ximg)*vhat)**2)/(2*taupll))
            if (rannyu().LT.pacc) exit

            cnt=cnt+1
            if (cnt.GE.1d3) then
                write(6,'(''warning: early exit from rejection method loop'')')
                write(0,'(''warning: early exit from rejection method loop'')')
            endif
        end do
        ddiffusesq=sum(xdif**2)
    end subroutine sample_1e_move_taudif
!------------------------------------------------------------------------------!

    subroutine sample_1e_move_psia(tau,grado,hesso,laplo,xo,ddiffusesq,xn)
        use contrldmc_mod, only: l_method_of_images
        implicit none
        real(dp), intent(in)     :: tau, grado(:), hesso(:,:), laplo, xo(:)
        real(dp), intent(out)    :: ddiffusesq, xn(:)
        real(dp)    :: rej_method_factor
        real(dp)    :: vhat(size(xo)), d2psio
        real(dp)    :: a, c, cprp
        real(dp)    :: taupll, tauprp, normpll, normprp, pdif, pacc
        real(dp)    :: gauss, rannyu
        real(dp)    :: xpeak(size(xo)), xdif(size(xo)), xpll(size(xo)), dx(size(xo))
        real(dp)    :: greensfn
        integer     :: cnt, k

        rej_method_factor=2d0
        call fit_psia(grado,hesso,a,c)
        xpeak=a+(-a-sign(dsqrt(a**2+4*tau*(c*tau+1)),a))/(2*(c*tau+1))
        vhat=grado/norm2(grado)
        d2psio=sum(vhat*matmul(hesso,vhat))
        cprp=(d2psio-laplo)/(ndim-1)
        taupll=1/(c+1/tau)
        tauprp=1/(cprp+1/tau)
        if ((tauprp.LE.0).OR.(tauprp.GE.3*tau)) then
          cprp=0d0
          tauprp=tau
        endif
        normpll=dsqrt(2*pi*taupll)
        normprp=dsqrt(2*pi*tauprp)**(size(xo)-1)
        do while (.TRUE.)
          xdif=[(gauss(),k=1,size(xo))]
          xpll=sum(xdif*vhat)*vhat
          pdif=dexp(-sum(xdif**2)/2)/(normpll*normprp)
          xdif=dsqrt(taupll)*xpll+dsqrt(tauprp)*(xdif-xpll)
          dx=xpeak*vhat+xdif
          greensfn=psia_greensfn(laplo,d2psio,grado,a,c,dx)
          xn=xo+xpeak*vhat+xdif
          pacc=greensfn/(rej_method_factor*pdif)

          if (pacc.GT.1) then
            rej_method_factor=1.1d0*greensfn/pdif
            write(6,'(''warning: rej_method_factor too small'')')
            write(6,'(''rej_method_factor reset to '', f13.5)') &
                rej_method_factor
            write(0,'(''warning: rej_method_factor too small'')')
            write(0,'(''rej_method_factor reset to '', f13.5)') &
                rej_method_factor
          endif

          if (rannyu().LT.pacc) exit

          cnt=cnt+1
          if (cnt.GE.1e3) then
            write(6,'(''warning: early exit from rejection method loop'')')
            write(0,'(''warning: early exit from rejection method loop'')')
            exit
          endif
        enddo
        ddiffusesq=sum(xdif**2)
    end subroutine sample_1e_move_psia

    real(dp) function sample_1e_prob_drift_diffusion(tau,grado,hesso,xo,xn)
        implicit none
        real(dp), intent(in)     :: tau, grado(:), hesso(:,:), xo(:), xn(:)
        real(dp)                 :: vav(size(xo))

        vav=averaged_1e_velocity(tau,grado,hesso)
        sample_1e_prob_drift_diffusion=&
            dexp(-sum((xn-xo-vav*tau)**2)/(2*tau))/dsqrt(2*pi*tau)**size(xo)
    end function sample_1e_prob_drift_diffusion
!------------------------------------------------------------------------------!

    real(dp) function sample_1e_prob_drift_diffusion_with_cusp(tau,grado,hesso,xo,xen,ren,xn)
        use atom_mod, only: cent, znuc, iwctype
        implicit none
        real(dp), intent(in)     :: tau, grado(:), hesso(:,:), xo(:), xen(:,:), ren(:), xn(:)
        real(dp)                 :: z(size(xo)), x(size(xo)), zi, zf, xf
        real(dp)                 :: vav(size(xo)), xdrift(size(xo))
        real(dp)                 :: pgaus, drift_diffusion_prob, zeta
        integer                  :: iwnuc
!Find index of nucleus closest to electron being moved
        iwnuc=minloc(ren, DIM=1)
!Place z-axis along direction from nearest nucleus to electron and
!x-axis along direction of angular component of velocity.
        z=xen(:,iwnuc)
        z=z/norm2(z)
        x=grado-sum(grado*z)*z
        if (norm2(x).LT.eps) then
            x(1)=eps*(1-z(1)**2)
            x(2)=eps*(-z(1)*z(2))
            x(3)=eps*(-z(1)*z(3))
        endif
        x=x/norm2(x)
!Calculate averaged velocity  
        vav=averaged_1e_velocity(tau,grado,hesso,xen)
!Prob. of sampling exponential rather than gaussian is
!half*derfc(rtry/dsqrt(two*tau)) = half*(two-derfc(-rtry/dsqrt(two*tau)))
!We use both expressions because under AIX the former is faster if rtry>0
!and the latter is faster if rtry<0.
!Note that if adrift is always set to 1 then it may be better to use
!vavvt rather than tau since the max drift distance is dsqrt(2*tau/adrift),
!so both the position of the gaussian and its width are prop to dsqrt(tau)
!if tau is used in derfc, and so qgaus does not tend to 1 for large tau.
!However, we are using a variable adrift that can be very small and then
!using tau is the better choice.
        zi=ren(iwnuc)
        zf=zi+sum(vav*z)*tau
        if (zf.GE.0) then
            pgaus=1-derfc(zf/dsqrt(2*tau))/2
        else
            pgaus=derfc(-zf/dsqrt(2*tau))/2
        endif
!The electron cannot drift past the nucleus
        zf=max(0d0,zf)
!The drift bends toward the nucleus
        xf=sum(vav*x)*tau*min(1d0,zf/((zf+zi)/2))
        xdrift=cent(:,iwnuc)+xf*x+zf*z
!Sample from drift-diffusion Green's function with probability pgaus.
!With probability 1-pgaus, sample from exponential centered at nearest nucleus.
        drift_diffusion_prob=dexp(-sum((xn-xdrift)**2)/(2*tau))/dsqrt(2*pi*tau)**size(xo)
        zeta=dsqrt(1/tau+znuc(iwctype(iwnuc))**2)
        !rttau=dsqrt(tau)
        !zeta=znuc(iwctype(iwnuc))+(0.34/rttau)*(1+(zi/rttau)/(1+zi/rttau))
        sample_1e_prob_drift_diffusion_with_cusp=&
            pgaus*drift_diffusion_prob + (1-pgaus)*zeta**3*dexp(-2*zeta*norm2(xn-cent(:,iwnuc)))/pi
!        sample_1e_prob_drift_diffusion=&
!            dexp(-sum((xn-xo-vav*tau)**2)/(2*tau))/dsqrt(2*pi*tau)**size(xo)
    end function sample_1e_prob_drift_diffusion_with_cusp
!------------------------------------------------------------------------------!

    real(dp) function sample_1e_prob_taudif(tau,grado,hesso,laplo,xo,xn)
        use contrldmc_mod, only: l_method_of_images
        use debug_shared_mod, only: iw, is
        implicit none
        real(dp), intent(in) :: tau, grado(:), hesso(:,:), laplo, xo(:), xn(:)
        real(dp) :: vnorm, vhat(size(xo)), vav(size(xo))
        real(dp) :: d2psi, divpll, divprp, dnode
        real(dp) :: taupll, tauprp
        real(dp) :: xdif(size(xo)), xprp(size(xo)), ximg(size(xo))
        real(dp) :: prob, pdiffo, pdiffi, normpll, normprp, normimg

        vnorm=norm2(grado)
        vhat=grado/vnorm
        d2psi=sum(vhat*matmul(hesso,vhat))
        if (l_method_of_images) then
            if (d2psi.LE.0) then
                dnode=abs((-vnorm+dsqrt(vnorm**2-2*d2psi))/d2psi)
            else
                dnode=1d3
            endif
        endif
        divpll=d2psi-sum(grado**2)
        divprp=(laplo-d2psi)/(size(xo)-1)
        taupll=tau_diffusion(tau,divpll)
        tauprp=tau_diffusion(tau,divprp)
        vav=averaged_1e_velocity(tau,grado,hesso)
        xdif=xn-xo-vav*tau
        xprp=xdif-sum(xdif*vhat)*vhat
        ximg=xo-2*dnode*vhat-vav*tau
        normpll=dsqrt(2*pi*taupll)
        normprp=dsqrt(2*pi*tauprp)**(size(xo)-1)
        prob=dexp(-sum(xprp**2)/(2*tauprp))/normprp
        if (l_method_of_images) then
            pdiffo=dexp(-sum(    xdif *vhat)**2/(2*taupll))/normpll
            pdiffi=dexp(-sum((xn-ximg)*vhat)**2/(2*taupll))/normpll
            normimg=erf((dnode+norm2(vav*tau))/dsqrt(2*taupll))
            prob=prob*max(0d0,pdiffo-pdiffi)/normimg
        else
            pdiffo=dexp(-sum(xdif*vhat)**2/(2*taupll))/normpll
            prob=prob*pdiffo
        endif
        sample_1e_prob_taudif=prob
    end function sample_1e_prob_taudif
!------------------------------------------------------------------------------!

    real(dp) function sample_1e_prob_psia(tau,grado,hesso,laplo,xo,xn)
        implicit none
        real(dp), intent(in) :: tau, grado(:), hesso(:,:), laplo, xo(:), xn(:)
        real(dp) :: a, c, d2pll

        call fit_psia(grado,hesso,a,c)
        d2pll=sum(grado*matmul(hesso,grado))/sum(grado**2)
        sample_1e_prob_psia=psia_greensfn(laplo,d2pll,grado,a,c,xn-xo)
    end function sample_1e_prob_psia
!------------------------------------------------------------------------------!

    real(dp) function ene_int_rewt(etrial,eest,esigma,tau,tauint,psio,psin)
        use control_mod, only: limit_rewt_dmc, c_rewt, ene_int, l_print_rewt
        use iterat_mod, only: ipass
        use contrldmc_mod, only: rtrttau, pow_rewt
        implicit none
        real(dp), intent(in)    :: etrial, eest, esigma, tau, tauint
        type(psi_t), intent(in) :: psio, psin
        real(dp)                :: eint, emin, ecut, edifo, edifn, ewto, ewtn
        real(dp)                :: ecuto, ecutn
        real(dp)                :: facto, factn, frato, fratn, frewto, frewtn
        real(dp)                :: fo, fn, so, sn
        real(dp)                :: vavo(size(psio%grad,1),size(psio%grad,2)) 
        real(dp)                :: vavn(size(psin%grad,1),size(psin%grad,2))
        integer                 :: ipr_sav=0
        real(dp) :: vav2sumn, v2sumn, vav2sumo, v2sumo

        vavo=averaged_velocity(tau,psio%grad,psio%hess,psio%rvec_en)
        vavn=averaged_velocity(tau,psin%grad,psin%hess,psin%rvec_en)
        edifo=psio%eloc-eest
        edifn=psin%eloc-eest
!        ecuto=max(edifo,-limit_rewt_dmc*esigma)
!        ecutn=max(edifn,-limit_rewt_dmc*esigma)

        if     (index(ene_int,'ene_int_v').NE.0) then
!            facto=max(dsqrt(sum(psio%grad**2)/nelec)*tauint*c_rewt,1e-9)
!            factn=max(dsqrt(sum(psin%grad**2)/nelec)*tauint*c_rewt,1e-9)
            facto=max(dsqrt(sum(psio%grad**2))/(nelec**pow_rewt)*tauint*c_rewt,1e-9)
            factn=max(dsqrt(sum(psin%grad**2))/(nelec**pow_rewt)*tauint*c_rewt,1e-9)
        elseif (index(ene_int,'ene_int_e').NE.0) then
!            ecuto=max(edifo,-limit_rewt_dmc*esigma)
!            ecutn=max(edifn,-limit_rewt_dmc*esigma)
!            facto=max(abs(ecuto)*tauint*c_rewt/esigma,1e-9)
!            factn=max(abs(ecutn)*tauint*c_rewt/esigma,1e-9)
            facto=max(abs(edifo)*tauint*c_rewt/esigma,1e-9)
            factn=max(abs(edifn)*tauint*c_rewt/esigma,1e-9)
        endif

        if     (ene_int=='ene_int_v9') then
            ewto=eest+edifo*(dsqrt(pi)/2)*erf(facto)/facto
            ewtn=eest+edifn*(dsqrt(pi)/2)*erf(factn)/factn
            frewto=(dsqrt(pi)/2)*erf(facto)/facto
            frewtn=(dsqrt(pi)/2)*erf(factn)/factn
        elseif (ene_int=='ene_int_e9') then
            ewto=eest+edifo*(dsqrt(pi)/2)*erf(facto)/facto
            ewtn=eest+edifn*(dsqrt(pi)/2)*erf(factn)/factn
        elseif (ene_int=='alfe') then
            ecut=0.2d0*dsqrt(nelec/tau)
            ewto=eest+min(ecut,max(-ecut,edifo))
            ewtn=eest+min(ecut,max(-ecut,edifn))
        elseif (ene_int=='no_ene_int') then
            ewto=psio%eloc
            ewtn=psin%eloc
        elseif (ene_int=='unr93') then
            frato=dsqrt(sum(vavo**2)/sum(psio%grad**2))
            fratn=dsqrt(sum(vavn**2)/sum(psin%grad**2))
            frewto=frato
            frewtn=fratn
            ewto=eest+edifo*frato
            ewtn=eest+edifn*fratn
        else
            write(6,'(''ene_int not set correctly in input'')')
            stop 'ene_int not set correctly in input'
        endif

        if (l_print_rewt) then
            write(6,'(''ifrag, -edif, frewt, vloc'',1I3,3d12.4),'':''') 0, -edifn, frewtn, dsqrt(sum(psin%grad**2))
            write(6,'(''ifrag, -edif, frewt, vloc'',1I3,3d12.4),'':''') 0, -edifo, frewto, dsqrt(sum(psio%grad**2))
        endif

        eint=(ewto+ewtn)/2-etrial
!        fo=fcut_UNR93(tauint,psio%grad,psio%hess,psio%rvec_en)
!        fn=fcut_UNR93(tauint,psin%grad,psin%hess,psin%rvec_en)
!        so=etrial-eest+(eest-psio%eloc)*fo
!        sn=etrial-eest+(eest-psin%eloc)*fn
!        eint=-(so+sn)/2

        emin=eest-etrial-limit_rewt_dmc*esigma/rtrttau
        if (eint.LT.emin) then
            ipr_sav=ipr_sav+1
            if(ipr_sav.LE.3) then
                write(6,'(''Warning: eint<eest-etrial-limit_rewt_dmc*esigma/rtrttau: esigma,eint,ewto,ewtn,fratio(iw,ifr),fration='',6d12.4)') &
     &          esigma,eint,ewto,ewtn,frato,fratn
                if(ipr_sav.EQ.1) write(6,'(''This should add a totally negligible positive bias to the energy'')')
            elseif(ipr_sav.EQ.4) then
                write(6,'(''Warning: Additional warning msgs. of eint<eest-etrial-limit_rewt_dmc*esigma suppressed'')')
            endif
            eint=emin
        endif

        ene_int_rewt=eint
    end function ene_int_rewt
!------------------------------------------------------------------------------!

    real(dp) function limit_ene_int(eint,emin)
        implicit none
        real(dp), intent(in) :: eint, emin
        integer :: ipr_sav=0

        if (eint.LT.emin) then
            ipr_sav=ipr_sav+1
            if(ipr_sav.LE.3) then
                write(6,'(''Warning: eint<emin=eest-etrial-limit_rewt_dmc*esigma): eint,emin='',6d12.4)') &
     &          eint,emin
                if(ipr_sav.EQ.1) write(6,'(''This should add a totally negligible positive bias to the energy'')')
            elseif(ipr_sav.EQ.4) then
                write(6,'(''Warning: Additional warning msgs. of eint<eest-etrial-limit_rewt_dmc*esigma suppressed'')')
            endif
            limit_ene_int=emin
        else
            limit_ene_int=eint
        endif
    end function limit_ene_int
!------------------------------------------------------------------------------!

    real(dp) function ene_int2(etrial,eest,esigma,tau,tauint,psi,ifrag)
        use fragments_mod, only: l_fragments, iwfragelec, nelecfrag, nfrag
        use control_mod, only: limit_rewt_dmc, c_rewt_frag, ene_int, l_print_rewt
        use const_mod, only: nelec
        use contrldmc_mod, only: pow_rewt
        implicit none
        real(dp), intent(in) :: etrial, eest, esigma, tau, tauint
        type(psi_t), intent(in) :: psi
        integer, intent(in) :: ifrag
        real(dp) :: edif, fact, ewt, ecut, frat, frewt
        real(dp), allocatable :: grad(:,:), hess(:,:,:), vav(:,:)
        integer :: ie, nel, i

        if (ifrag.EQ.nfrag+1) then
            edif=psi%enefrag(ifrag)-eest
            fact=abs(edif)*tauint*c_rewt_frag(ifrag)/esigma
            frewt=(dsqrt(pi)/2)*erf(fact)/fact
!            ewt=eest+edif*(dsqrt(pi)/2)*erf(fact)/fact
            ewt=eest+edif*frewt
            ene_int2=ewt-etrial
            if (l_print_rewt) then
                write(6,'(''ifrag, -edif, frewt, abs(edif)/esigma'',1I3,3d12.4),'':''') ifrag, -edif, frewt, abs(edif)/esigma
            endif
            return
        endif
        
        if (l_fragments) then
!            nelecfrag=0
!            do ie=1,nelec
!                if (iwfragelec(ie).NE.ifrag) cycle
!                nelecfrag=nelecfrag+1
!            enddo
            nel=nelecfrag(ifrag)
        else
            nel=nelec
        endif
        allocate(grad(size(psi%grad,1),nel))
        allocate(hess(size(psi%hess,1),size(psi%hess,2),nel))

        edif=psi%enefrag(ifrag)-eest
        i=0
        do ie=1,nelec
            if (l_fragments.AND.(psi%iwfragelec(ie).NE.ifrag)) cycle
            i=i+1
            grad(:,i)  =psi%grad(:,ie)
            hess(:,:,i)=psi%hess(:,:,ie)
        enddo

        if     (ene_int=='ene_int_v9') then
            fact=max(dsqrt(sum(grad**2))/(nel**pow_rewt)*tauint*c_rewt_frag(ifrag),1e-9)
            ewt=eest+edif*(dsqrt(pi)/2)*erf(fact)/fact
            frewt=(dsqrt(pi)/2)*erf(fact)/fact
        elseif (ene_int=='alfe') then
            ecut=0.2d0*dsqrt(nel/tau)
            ewt=eest+min(ecut,max(-ecut,edif))
        elseif (ene_int=='no_ene_int') then
            ewt=psi%enefrag(ifrag)
        elseif (ene_int=='unr93') then
            allocate(vav(size(grad,1),size(grad,2)))
            vav=averaged_velocity(tau,grad,hess,psi%rvec_en)
            frat=dsqrt(sum(vav**2)/sum(grad**2))
            ewt=eest+edif*frat
            frewt=frat
        else
            write(6,'(''ene_int not set correctly in input'')')
            stop 'ene_int not set correctly in input'
        endif

        if (l_print_rewt) then
            write(6,'(''ifrag, -edif, frewt, vloc'',1I3,3d12.4),'':''') ifrag, -edif, frewt, dsqrt(sum(grad**2))
        endif
        ene_int2=ewt-etrial
    end function ene_int2
!------------------------------------------------------------------------------!

    real(dp) function fcut_UNR93(taunow,grad,hess,xen)
        implicit none
        real(dp), intent(in) :: taunow, grad(:,:), hess(:,:,:), xen(:,:,:)

        fcut_UNR93=dsqrt(sum(averaged_velocity(taunow,grad,hess,xen)**2)/sum(grad**2))
    end function fcut_UNR93
!------------------------------------------------------------------------------!

function averaged_velocity(tau,grad,hess,xen)
        implicit none
        real(dp), intent(in)    :: tau, grad(:,:), hess(:,:,:)
        real(dp), intent(in), optional :: xen(:,:,:)
        real(dp)                :: averaged_velocity(size(grad,1),size(grad,2))
        integer                 :: ie

        do ie=1,size(grad,2)
            if (present(xen)) then
                averaged_velocity(:,ie)=&
                averaged_1e_velocity(tau,grad(:,ie),hess(:,:,ie),xen(:,ie,:))
            else
                averaged_velocity(:,ie)=&
                averaged_1e_velocity(tau,grad(:,ie),hess(:,:,ie))
            endif
        enddo
    end function averaged_velocity
!------------------------------------------------------------------------------!

function averaged_1e_velocity(tau,grad,hess,xen)
        use atom_mod, only: znuc, iwctype
        use control_mod, only: adrift_default, n_drift_param
        use contrldmc_mod, only: l_modified_adrift
        implicit none
        real(dp), intent(in)    :: tau, grad(:), hess(:,:)
        real(dp), intent(in), optional :: xen(:,:)
        real(dp)    :: vnorm, vhat(size(grad)), adrift, vav
        real(dp)    :: averaged_1e_velocity(size(grad))
        real(dp)    :: xmin(size(grad)), rmin, hafzr2
        integer     :: iwnuc, ic
        real(dp)    :: dnuc,dmin

        vnorm=norm2(grad)
        vhat=grad/vnorm
        if (l_modified_adrift) then
            adrift=(vnorm**2-sum(vhat*matmul(hess,vhat)))/((n_drift_param-1)*vnorm**2)
            adrift=max(0.05d0, adrift)
        else if (present(xen).AND.(nloc.EQ.0)) then
            iwnuc=minloc(norm2(xen, DIM=1), DIM=1)
            xmin=xen(:,iwnuc)
            rmin=norm2(xmin)
            hafzr2=(znuc(iwctype(iwnuc))*rmin/2)**2
            adrift=(1+eps+sum(xmin*grad)/(vnorm*rmin))/2&
                    +adrift0*hafzr2/(1+hafzr2)
        else
            adrift=adrift_default
        endif
        if (adrift*vnorm**2*tau.LE.0d0) then
            vav=vnorm
        else
            !vav=(-1+dsqrt(1+2*adrift*vnorm**2*tau))/(adrift*vnorm*tau)
            vav=(-1+(1+n_drift_param*adrift*vnorm**2*tau)**(1d0/n_drift_param))/(adrift*vnorm*tau)
        endif
        averaged_1e_velocity=vav*vhat
    end function averaged_1e_velocity
!------------------------------------------------------------------------------!

    function sho_green_function(etrial,taunow,xo,xn) result(gf)
        implicit none
        real(dp), intent(in)    :: etrial, taunow, xo(:,:), xn(:,:)
        integer                 :: ie, k
        real(dp)                :: gf, a, b
        
        gf=dexp(etrial*taunow)
        do ie=1,size(xo,2)
            do k=1,size(xo,1)
                a=xo(k,ie)
                b=xn(k,ie)
                gf=gf*dexp(-(cosh(taunow)*(a**2+b**2)-2*a*b)/(2*sinh(taunow)))/dsqrt(2*pi*sinh(taunow))
            enddo
        enddo
    end function sho_green_function
!------------------------------------------------------------------------------!

    real(dp) function tau_diffusion(tau,divv)
        implicit none
        real(dp), intent(in) :: tau, divv
        real(dp) :: tmp1, tmp2

        if (divv.ge.0) then
          tau_diffusion = tau
        else
          tmp1=(tau-1/divv)*dsqrt(-divv)*erf(1/sqrt(-2*tau*divv))
          tmp2=dsqrt(2*tau/pi)*dexp(1/(2*tau*divv))
          tau_diffusion=3*tau-1/divv-(tmp1+tmp2)**2
        endif
    end function tau_diffusion
!------------------------------------------------------------------------------!

    subroutine fit_psia(grad, hess, a, c)
        use constants_mod, only: dp
        use contrldmc_mod, only: tau
        implicit none
        real(dp), intent(in) :: grad(:), hess(:,:)
        real(dp), intent(out) :: a, c
        real(dp) :: vnorm, v2, vhat(size(grad))
        real(dp) :: d2psi, rad
        real(dp) :: anear, afar
        real(dp) :: cnear, cfar
        real(dp) :: cmin, cmax
    
        vnorm=norm2(grad)
        vhat=grad/vnorm
        d2psi=sum(vhat*matmul(hess,vhat))
        v2=vnorm**2
        rad=v2+8*(v2-d2psi)
        vnorm = dsqrt(v2)
        if (rad.LT.0) then
          a=-1d0/vnorm
          c=0d0   
          return
        endif
        cmax=10d0
        cmin=-0.5d0/tau
        anear=-(-vnorm+dsqrt(rad))/(2*(v2-d2psi))
        anear=(1+1.0d0*sqrt(tau)*min(anear**4,1.d0))*anear
        cnear=(1+vnorm*anear)/(anear**2)
        afar=-(-vnorm-dsqrt(rad))/(2*(v2-d2psi))
        afar=(1+1.0d0*sqrt(tau)*min(afar**4,1.d0))*afar
        cfar = (1 + vnorm*afar)/(afar**2)
        if ((0.LT.cnear).AND.(cnear.LT.cmax)) then
          a=anear
          c=cnear
        else if ((0.LT.cfar).AND.(cfar.LT.cmax)) then
          a=afar
          c=cfar
        else if ((cnear.LT.0).AND.(cfar.LT.0)) then
          if(abs(cnear).le.abs(cfar)) then
            a=anear
            c=max(cmin,cnear)
          else
            a=afar
            c=max(cmin,cfar)
          endif
        else
          a=anear
          c=max(cmin, min(cmax, cnear))
        endif
    end subroutine fit_psia
!------------------------------------------------------------------------------!

    function gamma_incomplete(a,x)
        use constants_mod, only: dp
        implicit none
        real(dp), intent(in) :: a, x
        real(dp) :: gamma_incomplete
        real(dp), external :: gammai
        integer :: iflag
        
        gamma_incomplete = -gammai(a,x,x**a*exp(-x),iflag)
        if (iflag.eq.1) gamma_incomplete = gamma(a)+gamma_incomplete
    end function gamma_incomplete
!------------------------------------------------------------------------------!

    subroutine psia_greensfn_norm(t, a, c, norm)
        use constants_mod !dp
        implicit none
    
        real(dp), intent(in) :: t, a, c
        real(dp), intent(out) :: norm
        real(dp) :: k1,k2,n1,n2,n3,gi
        
        k1 = c*t+1
        k2 = a**2/(2*t*k1)
        if (k1.LE.0d0) then
          write(5,'(''warning: psia is not integrable: tau, c'')') t, c
          write(6,'(''warning: psia is not integrable: tau, c'')') t, c
        endif
    !    n1 = dsqrt(2*pi*t)*exp((c*t/(c*t+1))*c*a**2/2)
        n1 = exp((c*t/(c*t+1))*c*a**2/2)
        if (k2.LT.1d2) then
          gi = gamma_incomplete(5d-1, k2)/dsqrt(pi)
        else !if k2 >= 1d2, gamma_incomplete(5d-1, k2) <= 1d-43
          gi = 0d0 
        endif
        n2 = k1*(erf(dsqrt(k2))-1+gi)+1
        n3 = k1**1.5
        norm = n1*n2/n3
    end subroutine psia_greensfn_norm
!------------------------------------------------------------------------------!

    subroutine psia_greensfn_norm_old(t, a, c, norm)
        use constants_mod !dp
        implicit none
    
        real(dp), intent(in) :: t, a, c
        real(dp), intent(out) :: norm
        real(dp) :: k1,k2,n1,n2,n3,gi
        
        k1 = c*t+1
        k2 = a**2/(2*t*k1)
    !    n1 = dsqrt(2*pi*t*a**2)*exp(-(k1-1)*k2)
        n1 = abs(a)*exp(-(k1-1)*k2)
        if (n1.EQ.0d0) then !gamma_incomplete takes a long time calculating zero with large enough args (k2 ~= 1E+11)
            norm = 0d0
            return
        endif
        n2 = k1*(erf(dsqrt(k2)) + gamma_incomplete(5d-1, k2)/dsqrt(pi)-1)+1
        n3 = k1**1.5
        norm = n1*n2/n3
    end subroutine psia_greensfn_norm_old
!------------------------------------------------------------------------------!

    real(dp) function psia_greensfn(lap,d2pll,grad,a,c,dx)
        use constants_mod, only: dp
        use contrldmc_mod, only: tau, l_psi_approx_norm_old
        implicit none
        real(dp), intent(in) :: lap, d2pll, a, c, grad(:), dx(:)
        real(dp) :: xpll, xprp2, xpeak
        real(dp) :: psia_pll, psia_prp
        real(dp) :: gauss_xo, gauss_im, vhat(3), avg, var
        real(dp) :: cprp, tau_prp, norm_prp, norm_pll
    
        vhat = grad/norm2(grad)
        xpll = dot_product(dx, vhat)
        xprp2 = sum(dx**2) - xpll**2

        if (l_psi_approx_norm_old) then
          psia_pll = abs(xpll-a)*exp(-c*(xpll-a)**2/2)
          if (a.lt.0) then !reject moves across estimated node
              if (xpll.lt.a) psia_pll = 0d0
          else
              if (xpll.gt.a) psia_pll = 0d0
          endif
        else
          psia_pll = max(0d0, (1-xpll/a)*exp(-c*(xpll-a)**2/2 + c*a**2/2))
        endif
        cprp = (d2pll-lap)/(size(grad)-1) 
!        if (.not.psi_approx_prp) cprp = 0d0
        tau_prp = 1d0/(cprp + 1d0/tau)
        if ((tau_prp.LE.0).OR.(tau_prp.GE.3*tau)) then
          cprp = 0d0
          tau_prp = tau
        endif
        psia_prp = exp(-cprp*xprp2/2)
        norm_prp = 1d0/(cprp*tau+1d0)
    
        gauss_xo = exp(-sum(dx**2)/(2*tau))/(2*pi*tau)**1.5d0
        gauss_im = exp(-((xpll-2*a)**2 + xprp2)/(2*tau))/(2*pi*tau)**1.5d0
        
        !TODO: psia_greensfn_norm can be slow. Should speed up gammai
        if (l_psi_approx_norm_old) then
          call psia_greensfn_norm_old(tau, a, c, norm_pll)
        else
          call psia_greensfn_norm(tau, a, c, norm_pll)
        endif
        if (norm_pll.GT.0d0) then
          psia_greensfn=&
              psia_pll*psia_prp*(gauss_xo - gauss_im)/(norm_pll*norm_prp)
        else !in the large 'a' limit, the sampled distribution is a gaussian
          avg = c*a*tau/(c*tau + 1)
          var = tau/(c*tau + 1)
          psia_greensfn=&
              psia_prp*dexp(-(xpll-avg)**2/(2*var)-xprp2/(2*tau))/(dsqrt(2*pi*var)*2*pi*tau_prp)
        endif
    end function psia_greensfn
end module green_functions_mod

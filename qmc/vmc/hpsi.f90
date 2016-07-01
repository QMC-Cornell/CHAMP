!     subroutine hpsi(coord,psid,psij,energy,d2psi,pe,velocity,div_v,ifr)
      subroutine hpsi(coord,psid,psij,velocity,div_v,d2psi,pe,pei,energy,denergy,ifr)
! Written by Cyrus Umrigar
! adapted to complex wavefunctions by A.D.Guclu, Feb2004.
! Calculates determinantal and Jastrow parts of psi, velocity, divergence of V,
! Laplacian of psi, potential energy, total energy and derivatives of the energy
! wrt the wavefunction parameters.

! Note that for a nonlocal psp. the potential energy at a given MC config depends on the wavefn. params. So:
! For a non-optimization run
! --------------------------
!      distances
!      |
!      jastrow
!      |             getvps_xxx
! hpsi /            /        nonlocd
!      \nonloc_pot---nonloc-/
!                           \nonlocj

! For an optimization run
! -----------------------
!      distances
!      |
!      deriv_jastrow
!      |                  getvps_xxx
! hpsi /                 /               nonlocd
!      \deriv_nonloc_pot---deriv_nonloc-/
!                                       \deriv_nonlocj
! Note there is no deriv_nonlocd because the additional object needed for CSF optimization is calculated in nonlocd.

      use all_tools_mod
      use control_mod
      use fitdet_mod
      use eloc_mod
      use psi_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use contrl_opt2_mod
      use wfsec_mod
      use vj_mod
      use pseudo_mod
      use contrl_per_mod
      use deriv_csf_mod
      use derivjas_mod
      use contr3_mod
      use distance_mod
      use vd_mod
      use optimo_mod
      use contrl_opt_mod
      implicit real*8(a-h,o-z)

! complex local:
      complex*16 cvd(3,nelec),cvk(3,nelec),cvelocity(3,nelec)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot
      dimension coord(3,*),velocity(3,nelec)
      dimension div_vj(nelec),div_vk(nelec),div_vd(nelec),div_v(nelec),dpe(nparm),denergy(nparm)

      iwf=iwftype(ifr)

! Distances are needed for Jastrow, determinant(e-n only), and for potential energy.
! pe_ee is always computed in distances. (pot_ee is called from distances if iperiodic != 0)
! pe_en is computed in distances if nloc=0.
! pe_en is computed in nonloc_pot if nloc !=0 and iperiodic=0.
! pe_en(loc) is computed in distances and pe_en(nonloc) in nonloc_pot if nloc !=0 and iperiodic!=0.
! (pot_en is called from distances if iperiodic != 0)
      call distances(coord,pe,pei)

! get contribution from jastrow
      if(ianalyt_lap.eq.1) then
!       if((igradhess.eq.0 .and. index(mode,'fit').eq.0) .or. ifr.gt.1) then
        if((igradhess.eq.0 .or. ifr.gt.1) .and. .not. l_opt_jas) then
           call jastrow(coord,vj,d2j,div_vj,psij)
        else
          call deriv_jastrow(coord,vj,d2j,div_vj,psij)
       endif
      else
!       if((igradhess.eq.0 .and. index(mode,'fit').eq.0) .or. ifr.gt.1) then
        if((igradhess.eq.0 .or. ifr.gt.1) .and. .not. l_opt_jas) then
          call jastrow_num(coord,vj,d2j,div_vj,psij)
         else
          stop 'cannot do derivatives wrt wavefn. parameters when doing numerical laplacian'
        endif
      endif

!      write(6,'(''After calling jastrow or jastrow_num, d2j = '', f12.8)') d2j
!      write(6,'(''div_vj(i)='',90f10.6)') (div_vj(i),i=1,nelec)
!      write(6,'(''x(1,i)='',90f6.2)') (coord(1,i),i=1,nelec)
!      write(6,'(''x(2,i)='',90f6.2)') (coord(2,i),i=1,nelec)
!      write(6,'(''v(1,i)='',90f10.6)') (vj(1,i),i=1,nelec)
!      write(6,'(''v(2,i)='',90f10.6)') (vj(2,i),i=1,nelec)
      call object_modified_by_index (vj_index)  !JT
      sum_lap_lnj = d2j            !JT
      call object_modified_by_index (sum_lap_lnj_index) !JT

! get contribution from determinants
!   first check if we need to call determinant()/cdeterminant()
      idodet=1
      if(index(mode,'fit').ne.0 .and. nloc.le.0) then
        if(isaved.eq.1) idodet=0
      endif

      if(ibasis.eq.3) then
! note:
! d2psi = real(lap(psi)) - dabs(cv)**2 = jackson feenberg kin.en.
! velocity = real(cv)

        if(idodet.eq.1) then    ! this part should go in a subroutine..
          if(idot.eq.0 .or. idot.eq.1) then
            call cdeterminant(coord,rvec_en,r_en,cvd,d2d,div_vd,psid)
          elseif(idot.eq.3) then
            call cdeterminant_cf(coord,rvec_en,r_en,cvd,d2d,div_vd,psid)
          else                 ! no determinant for laughlin wfs
            d2d=0.d0
            psid=1.d0
            do 3 i=1,nelec
              div_vd(i)=0.d0
              do 3 k=1,ndim
    3           cvd(k,i)=dcmplx(0.d0,0.d0)
          endif
! save only if we are doing fit"
          if(index(mode,'fit').ne.0) then
            d2d_sav(iconfg)=d2d
            psid_sav(iconfg)=psid
            do 5 i=1,nelec
              div_vd_sav(i,iconfg)=div_vd(i)
              do 5 k=1,ndim
    5           cvd_sav(k,i,iconfg)=cvd(k,i)
          endif

          if(idot.gt.0) then
            call cjastrow(cvk,d2k,div_vk,psik)
            if(index(mode,'fit').ne.0) then
              psik_sav(iconfg)=psik
              d2k_sav(iconfg)=d2k
              do 6 i=1,nelec
                div_vk_sav(i,iconfg)=div_vk(i)
                do 6 k=1,ndim
    6             cvk_sav(k,i,iconfg)=cvk(k,i)
            endif
          endif

         else

          d2d=d2d_sav(iconfg)
          psid=psid_sav(iconfg)
          do 7 i=1,nelec
            div_vd(i)=div_vd_sav(i,iconfg)
            do 7 k=1,ndim
    7          cvd(k,i)=cvd_sav(k,i,iconfg)
          if(idot.gt.0) then
            psik=psik_sav(iconfg)
            d2k=d2k_sav(iconfg)
            do 8 i=1,nelec
              div_vk(i)=div_vk_sav(i,iconfg)
              do 8 k=1,ndim
    8           cvk(k,i)=cvk_sav(k,i,iconfg)
          endif
        endif
! combine to get derivatives of full wavefunction
! Gradient and Laplacian of ln(JD) are the sum of the gradient and Laplacian of ln(J) and ln(D).
        do 10 i=1,nelec
          div_v(i)=div_vj(i)+div_vd(i)
          do 10 k=1,ndim
            if(ipr.ge.2) write(6,'(''ifr,cvd,vj'',i2,9d12.5)') ifr,cvd(k,i),vj(k,i)
            cvelocity(k,i)=dcmplx(vj(k,i),0.d0)+cvd(k,i)
   10   continue
        d2psi=d2d+d2j

! calculate extra pieces due to composite fermion vortices
        if(idot.gt.0) then
          psid=psid*psik
          d2psi=d2psi+d2k
          do 13 i=1,nelec
            div_v(i)=div_v(i)+div_vk(i)
            do 13 k=1,ndim
              cvelocity(k,i)=cvelocity(k,i)+cvk(k,i)
! complex jastrow factor complicates calculation of jackson feenberg energy
              temp=2*dimag(cvk(k,i))*(dimag(cvk(k,i))+2*dimag(cvd(k,i)))
              d2psi=d2psi-temp
              div_v(i)=div_v(i)-temp
   13     continue
        endif

! calculate velocity
        velocity2=0.d0
        do 14 i=1,nelec
          do 14 k=1,ndim
            velocity(k,i)=dreal(cvelocity(k,i))
            velocity2=velocity2+dreal(cvelocity(k,i))**2+dimag(cvelocity(k,i))**2
   14   continue

        if(ipr.ge.2) write(6,'(''coord'',9d12.5)') ((coord(k,ielec),k=1,ndim),ielec=1,nelec)
        if(ipr.ge.2) write(6,'(''rvec_en'',9d12.5)') ((rvec_en(k,ie,1),k=1,ndim),ie=1,nelec)
        if(ipr.ge.2) write(6,'(''ifr,psid,psij'',i2,9d12.5)') ifr,psid,psij
        if(ipr.ge.2) write(6,'(''ifr,d2d,d2j'',i2,9d12.5)') ifr,d2d,d2j
        if(ipr.ge.2) write(6,'(''pe,velocity2,d2psi,emaglz,emagsz'',9d12.5)') &
     &                     pe,pei,velocity2,d2psi,emaglz,emagsz

! calculate local energy. for convenience "magnetic energy" is included in pe.
! it will be separated while printing out the results...
        pe=pe+emag
!       pe=pe+emagsz+emagv+emaglz
        energy=pe-hb*(velocity2+d2psi)

        if(igradhess.ne.0) then

! calculate parameter-derivatives of csf_coefs for optimization
          do 20 iparm=1,nparmcsf+nparmot
            stop 'optimization of csf parameters not possible with complex wfs yet'
!           if(nloc.le.0) then
!             denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)
!           else
!             denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)+dpe(iparm)
!           endif
!           do 20 i=1,nelec
!             do 20 k=1,ndim
!  20           denergy(iparm)=denergy(iparm)-hb*(2*(ddeti_det(k,i,iparm)-deti_det(iparm)*vd(k,i))*vj(k,i))
   20       continue

! calculate parameter-derivatives of Jastrow factor for optimization
! if Psi=J*D, J=exp(f), V=grad(Psi)/Psi, E_L = H*Psi/Psi and _i denotes deriv wrt. parameter c_i, then
! E_L,i = -0.5*(\nabla^2*f_i + 2 grad(f_i).V) + (\hat{V}\psi/psi)_i
          do 25 iparm=1,nparmj+nparms
            if(nloc.le.0) then
              denergy(nparmcsf+nparmot+iparm)=-hb*d2g(iparm)
            else
              denergy(nparmcsf+nparmot+iparm)=-hb*d2g(iparm)+dpe(nparmcsf+nparmot+iparm)
            endif
            do 25 i=1,nelec
              do 25 k=1,ndim
   25           denergy(nparmcsf+nparmot+iparm)=denergy(nparmcsf+nparmot+iparm)-2*hb*(g(k,i,iparm)*velocity(k,i))

        endif

      else
        if(idodet.eq.1) then
          call determinant(coord,rvec_en,r_en,vd,d2d,div_vd,psid)
          call object_modified_by_index (vd_index) !WAS

          if(index(mode,'fit').ne.0) then
            d2d_sav(iconfg)=d2d
            psid_sav(iconfg)=psid
            do 30 i=1,nelec
              div_vd_sav(i,iconfg)=div_vd(i)
              do 30 k=1,ndim
   30           vd_sav(k,i,iconfg)=vd(k,i)
          endif
         else
          d2d=d2d_sav(iconfg)
          psid=psid_sav(iconfg)
          do 35 i=1,nelec
            div_vd(i)=div_vd_sav(i,iconfg)
            do 35 k=1,ndim
   35         vd(k,i)=vd_sav(k,i,iconfg)
        endif
        do 40 i=1,nelec
          div_v(i)=div_vj(i)+div_vd(i)
          do 40 k=1,ndim
            if(ipr.ge.2) write(6,'(''ifr,vd,vj'',i2,9d12.5)') ifr,vd(k,i),vj(k,i)
   40       velocity(k,i)=vj(k,i)+vd(k,i)
        d2psi=d2d+d2j

        if(ipr.ge.2) write(6,'(''coord'',9d12.5)') ((coord(k,ielec),k=1,ndim),ielec=1,nelec)
        if(ipr.ge.2) write(6,'(''ifr,psid,psij'',i2,9d12.5)') ifr,psid,psij
        if(ipr.ge.2) write(6,'(''ifr,d2d,d2j'',i2,9d12.5)') ifr,d2d,d2j

! get pseudo-potential contribution
! nonloc_pot must be called after determinant because psid is needed in nonloc_pot
! This call computes the entire pot_en if iperiodic=0 but only the nonlocal contrib. otherwise.
        if(ipr.ge.2) write(6,'(''hpsi: pe before nonloc_pot'',9f12.5)') pe
        if(nloc.gt.0) then
!         if((igradhess.eq.0 .and. index(mode,'fit').eq.0) .or. ifr.gt.1) then
          if((igradhess.eq.0 .or. ifr.gt.1) .and. .not. l_opt) then
            call nonloc_pot(coord,rshift,rvec_en,r_en,detu,detd,slmui,slmdi,vpsp,psid,pe)
          else
            call deriv_nonloc_pot(coord,rshift,rvec_en,r_en,detu,detd,deti_det,slmui,slmdi,vpsp,psid,pe,dpe,ifr)
          endif
        endif

        if(ipr.ge.3) write(6,'(''hpsi: pe after nonloc_pot'',9f12.5)') pe

      eloc_pot = pe
      call object_modified_by_index (eloc_pot_index) !JT
      call object_modified_by_index (deti_det_index) !JT

! calculate local energy. for convenience "magnetic energy" is included in pe.
! it will be separated while printing out the results...
        pe=pe+emag

! Laplacian(JD)/JD = Laplacian(ln(JD)) + V^2
        energy=pe-hb*d2psi
        do 45 i=1,nelec
          do 45 k=1,ndim
   45       energy=energy-hb*velocity(k,i)**2
        if(ipr.ge.3) write(6,'(''hpsi: ifr,energy='',i2,9f16.10)') ifr,energy,psid,psij,pe
        if(ipr.ge.3) write(6,'(''hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy='',i2,6f16.10,/)') &
     &  ifr,psid,exp(psij),d2psi,energy-pe,pe,energy

        if(igradhess.ne.0 .or. l_opt) then
           ! calculate parameter-derivatives of csf_coefs for optimization
           ! MJO98 6/16 - denergy below is calculated without E_local * deti_det
           ! because expanding H psi_j/psi and psi_j/psi *H psi/psi
           ! allows the potential terms (V) to cancel
          do 50 iparm=1,nparmcsf+1+nparmot
!             print*,'iparm',iparm
            if(nloc.le.0) then
              denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)
!             write(6,*) 'd2deti_det,deti_det,denergy,iparm=' ,d2deti_det(iparm),deti_det(iparm),denergy(iparm),iparm
            else
              if(nparmot.gt.0) stop 'orbital optimization not possible with nonlocal pseudopotentials yet'
              denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)+dpe(iparm)
!             write(6,'(''dpe(iparm)='',9d12.5)') dpe(iparm)
            endif
            do 50 i=1,nelec
              do 50 k=1,ndim
!                write(6,*) '++ ddeti_det,k,i,iparm=',ddeti_det(k,i,iparm),k,i,iparm
   50           denergy(iparm)=denergy(iparm)-hb*(2*(ddeti_det(k,i,iparm)-deti_det(iparm)*vd(k,i))*vj(k,i))

          if(ipr.ge.4) write(6,'(''denergy='',9f10.6)') (denergy(iparm),iparm=1,nparmcsf+nparmot)
!          write(6,'(''denergy='',20f10.6)') (denergy(iparm),iparm=1,nparmcsf+1)
          denergy0 = denergy(nparmcsf+1)

          call object_modified_by_index (denergy0_index)
! calculate parameter-derivatives of Jastrow factor for optimization
! if Psi=J*D, J=exp(f), V=grad(Psi)/Psi, E_L = H*Psi/Psi and _i denotes deriv wrt. parameter c_i, then
! E_L,i = -0.5*(\nabla^2*f_i + 2 grad(f_i).V) + (\hat{V}\psi/psi)_i
          do 55 iparm=1,nparmj+nparms
            if(nloc.le.0) then
              denergy(nparmcsf+nparmot+iparm)=-hb*d2g(iparm)
            else
              denergy(nparmcsf+nparmot+iparm)=-hb*d2g(iparm)+dpe(nparmcsf+nparmot+iparm)
            endif
            do 55 i=1,nelec
              do 55 k=1,ndim
   55           denergy(nparmcsf+nparmot+iparm)=denergy(nparmcsf+nparmot+iparm)-2*hb*(g(k,i,iparm)*velocity(k,i))

          call object_modified_by_index (denergy_index)

        endif

      endif

      return
      end

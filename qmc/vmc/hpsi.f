c     subroutine hpsi(coord,psid,psij,energy,d2psi,pe,velocity,div_v,ifr)
      subroutine hpsi(coord,psid,psij,velocity,div_v,d2psi,pe,pei,energy,denergy,ifr)
c Written by Cyrus Umrigar
c adapted to complex wavefunctions by A.D.Guclu, Feb2004.
c Calculates determinantal and Jastrow parts of psi, velocity, divergence of V,
c Laplacian of psi, potential energy, total energy and derivatives of the energy
c wrt the wavefunction parameters.

      use all_tools_mod
      use control_mod
      use fitdet_mod
      use eloc_mod
      use psi_mod
      use slater_mod
      use optim_mod
!     use periodic_jastrow_mod  !WAS

      implicit real*8(a-h,o-z)

      character*16 mode

c complex local:
      complex*16 cvd(3,MELEC),cvk(3,MELEC),cvelocity(3,MELEC)
c common:
c     complex*16 cvd_sav,cvk_sav

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl_per/ iperiodic,ibasis
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contr3/ mode
      common /contrl_opt2/ igradhess,iadd_diag_opt
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc

!JT      common /slater/ slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD)
!JT     &,fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
!JT     &,fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
!JT     &,detu(MDETUD),detd(MDETUD)
!JT     &,ddeti_deti(3,MELEC,MDETUD),d2edeti_deti(MELEC,MDETUD),deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD),d2deti_det(MPARMD),d2det_det
!JT     &,detij_det(MPARMD,MPARMD)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase
      common /optimo/ iwo(MORB,MOTYPE),nparmo(MOTYPE),nparmot,notype
      common /derivjas/ gvalue(MPARMJ),g(3,MELEC,MPARMJ),d2g(MPARMJ)
     &,go(MELEC,MELEC,MPARMJ)

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot
c Now this common is passed with use fitdet_mod
c     common /fitdet/ cvd_sav(3,MELEC,MDATA),vd_sav(3,MELEC,MDATA),psid_sav(MDATA)
c    &               ,d2d_sav(MDATA),div_vd_sav(MELEC,MDATA),cvk_sav(3,MELEC,MDATA),psik_sav(MDATA)
c    &               ,div_vk_sav(MELEC,MDATA),d2k_sav(MDATA),iconfg,isaved

      common /vj/ vj(3,MELEC)  !JT
c      dimension vd(3,MELEC),
      common /vd/ vd(3,MELEC)   !WAS
      dimension coord(3,*),velocity(3,MELEC)
      dimension div_vj(MELEC),div_vk(MELEC),div_vd(MELEC)
     &,div_v(MELEC),dpe(MPARM),denergy(MPARM)

      iwf=iwftype(ifr)

c Distances are needed for Jastrow, determinant(e-n only), and for potential energy.
c pe_ee is always computed in distances. (pot_ee is called from distances if iperiodic != 0)
c pe_en is computed in distances if nloc=0.
c pe_en is computed in nonloc_pot if nloc !=0 and iperiodic=0.
c pe_en(loc) is computed in distances and pe_en(nonloc) in nonloc_pot if nloc !=0 and iperiodic!=0.
c (pot_en is called from distances if iperiodic != 0)
      call distances(coord,pe,pei)

c get contribution from jastrow
      if(ianalyt_lap.eq.1) then
c       if((igradhess.eq.0 .and. index(mode,'fit').eq.0) .or. ifr.gt.1) then
        if((igradhess.eq.0 .or. ifr.gt.1) .and. .not. l_opt_jas) then
           call jastrow(coord,vj,d2j,div_vj,psij)
        else
          call deriv_jastrow(coord,vj,d2j,div_vj,psij)
       endif
      else
c       if((igradhess.eq.0 .and. index(mode,'fit').eq.0) .or. ifr.gt.1) then
        if((igradhess.eq.0 .or. ifr.gt.1) .and. .not. l_opt_jas) then
          call jastrow_num(coord,vj,d2j,div_vj,psij)
         else
          stop 'cannot do derivatives wrt wavefn. parameters when doing numerical laplacian'
        endif
      endif

      call object_modified_by_index (vj_index)  !JT
      sum_lap_lnj = d2j            !JT
      call object_modified_by_index (sum_lap_lnj_index) !JT

c get contribution from determinants
c   first check if we need to call determinant()/cdeterminant()
      idodet=1
      if(index(mode,'fit').ne.0 .and. nloc.le.0) then
        if(isaved.eq.1) idodet=0
      endif

      if(ibasis.eq.3) then
c note:
c d2psi = real(lap(psi)) - dabs(cv)**2 = jackson feenberg kin.en.
c velocity = real(cv)

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
c save only if we are doing fit"
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
c combine to get derivatives of full wavefunction
c Gradient and Laplacian of ln(JD) are the sum of the gradient and Laplacian of ln(J) and ln(D).
        do 10 i=1,nelec
          div_v(i)=div_vj(i)+div_vd(i)
          do 10 k=1,ndim
            if(ipr.ge.2) write(6,'(''ifr,cvd,vj'',i2,9d12.5)') ifr,cvd(k,i),vj(k,i)
            cvelocity(k,i)=dcmplx(vj(k,i),0.d0)+cvd(k,i)
   10   continue
        d2psi=d2d+d2j

c calculate extra pieces due to composite fermion vortices
        if(idot.gt.0) then
          psid=psid*psik
          d2psi=d2psi+d2k
          do 13 i=1,nelec
            div_v(i)=div_v(i)+div_vk(i)
            do 13 k=1,ndim
              cvelocity(k,i)=cvelocity(k,i)+cvk(k,i)
c complex jastrow factor complicates calculation of jackson feenberg energy
              temp=2*dimag(cvk(k,i))*(dimag(cvk(k,i))+2*dimag(cvd(k,i)))
              d2psi=d2psi-temp
              div_v(i)=div_v(i)-temp
   13     continue
        endif

c calculate velocity
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
        if(ipr.ge.2) write(6,'(''pe,velocity2,d2psi,emaglz,emagsz'',9d12.5)')
     &                     pe,pei,velocity2,d2psi,emaglz,emagsz

c calculate local energy. for convenience "magnetic energy" is included in pe.
c it will be separated while printing out the results...
        pe=pe+emag
c       pe=pe+emagsz+emagv+emaglz
        energy=pe-hb*(velocity2+d2psi)

        if(igradhess.ne.0) then

c calculate parameter-derivatives of csf_coefs for optimization
          do 20 iparm=1,nparmcsf+nparmot
            stop 'optimization of csf parameters not possible with complex wfs yet'
c           if(nloc.le.0) then
c             denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)
c           else
c             denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)+dpe(iparm)
c           endif
c           do 20 i=1,nelec
c             do 20 k=1,ndim
c  20           denergy(iparm)=denergy(iparm)-hb*(2*(ddeti_det(k,i,iparm)-deti_det(iparm)*vd(k,i))*vj(k,i))
   20       continue

c calculate parameter-derivatives of Jastrow factor for optimization
c if Psi=J*D, J=exp(f), V=grad(Psi)/Psi, E_L = H*Psi/Psi and _i denotes deriv wrt. parameter c_i, then
c E_L,i = -0.5*(\nabla^2*f_i + 2 grad(f_i).V) + (\hat{V}\psi/psi)_i
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

c get pseudo-potential contribution
c nonloc_pot must be called after determinant because psid is needed in nonloc_pot
c This call computes the entire pot_en if iperiodic=0 but only the nonlocal contrib. otherwise.
        if(ipr.ge.3) write(6,'(''hpsi: pe before nonloc_pot'',9f12.5)') pe
        if(nloc.gt.0) then
c         if((igradhess.eq.0 .and. index(mode,'fit').eq.0) .or. ifr.gt.1) then
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

c calculate local energy. for convenience "magnetic energy" is included in pe.
c it will be separated while printing out the results...
        pe=pe+emag

c Laplacian(JD)/JD = Laplacian(ln(JD)) + V^2
        energy=pe-hb*d2psi
        do 45 i=1,nelec
          do 45 k=1,ndim
   45       energy=energy-hb*velocity(k,i)**2
        if(ipr.ge.3) write(6,'(''ifr,energy='',i2,9f16.10)') ifr,energy,psid,psij,pe

      if(igradhess.ne.0 .or. l_opt) then

c calculate parameter-derivatives of csf_coefs for optimization
          do 50 iparm=1,nparmcsf+nparmot
            if(nloc.le.0) then
              denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)
c             write(*,*) 'd2deti_det,deti_det,denergy,iparm=' ,d2deti_det(iparm),deti_det(iparm),denergy(iparm),iparm
            else
              if(nparmot.gt.0) stop 'orbital optimization not possible with nonlocal pseudopotentials yet'
              denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)+dpe(iparm)
c             write(6,'(''dpe(iparm)='',9d12.5)') dpe(iparm)
            endif
            do 50 i=1,nelec
              do 50 k=1,ndim
c                write(*,*) '++ ddeti_det,k,i,iparm=',ddeti_det(k,i,iparm),k,i,iparm
   50           denergy(iparm)=denergy(iparm)-hb*(2*(ddeti_det(k,i,iparm)-deti_det(iparm)*vd(k,i))*vj(k,i))

          if(ipr.ge.4) write(6,'(''denergy='',9f10.6)') (denergy(iparm),iparm=1,nparmcsf+nparmot)

c calculate parameter-derivatives of Jastrow factor for optimization
c if Psi=J*D, J=exp(f), V=grad(Psi)/Psi, E_L = H*Psi/Psi and _i denotes deriv wrt. parameter c_i, then
c E_L,i = -0.5*(\nabla^2*f_i + 2 grad(f_i).V) + (\hat{V}\psi/psi)_i
          do 55 iparm=1,nparmj+nparms
            if(nloc.le.0) then
              denergy(nparmcsf+nparmot+iparm)=-hb*d2g(iparm)
            else
              denergy(nparmcsf+nparmot+iparm)=-hb*d2g(iparm)+dpe(nparmcsf+nparmot+iparm)
            endif
            do 55 i=1,nelec
              do 55 k=1,ndim
   55           denergy(nparmcsf+nparmot+iparm)=denergy(nparmcsf+nparmot+iparm)-2*hb*(g(k,i,iparm)*velocity(k,i))

        endif

      endif

      return
      end

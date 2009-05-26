      subroutine deriv_hpsi(coord,psid,psij,energy,denergy,ifr)
c No longer used.
c Written by Claudia Filippi, by modifying hpsi, written by Cyrus Umrigar
c Calculates energy, velocity and Laplacian, and derivative of the energy
c wrt the jastrow parameters.

      use control_mod
      use basic_tools_mod
      use fitdet_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use wfsec_mod
      use pseudo_mod
      use contrl_per_mod
      use derivjas_mod
      use distance_mod
      implicit real*8(a-h,o-z)

!      include '../vmc/vmc.h'
!      include '../vmc/pseudo.h'
!      include '../vmc/force.h'
!      include 'fit.h'

c complex
      complex*16 cvd(3,MELEC),cvk(3,MELEC),cvelocity(3,MELEC)
c     complex*16 cvd_sav,cvk_sav





      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot
c      common /fitdet/ cvd_sav(3,MELEC,MDATA),vd_sav(3,MELEC,MDATA),psid_sav(MDATA)
c     &               ,d2d_sav(MDATA),div_vd_sav(MELEC,MDATA),cvk_sav(3,MELEC,MDATA),psik_sav(MDATA)
c     &               ,div_vk_sav(MELEC,MDATA),d2k_sav(MDATA),iconfg,isaved

      dimension coord(3,MELEC),velocity(3,MELEC)
      dimension vj(3,MELEC),vd(3,MELEC),div_vj(MELEC),div_vd(MELEC),div_vk(MELEC)
     &,dpe(MPARM),denergy(MPARM)

      stop 'deriv_hpsi in fit should no longer be called (check, copied from Cyrus version)'

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
        call deriv_jastrow(coord,vj,d2j,div_vj,psij)
       else
        call jastrow_num(coord,vj,d2j,div_vj,psij)
      endif

c get contribution from determinants
c   first check if we need to call determinant() or cdeterminant()
      if(isaved.eq.1 .and. nloc.le.0) then
        idodet=0
       else
        idodet=1
      endif

      if(ibasis.eq.3) then

        if(idodet.eq.1) then               ! take care of each case 1by1 would make it more clear...
          if(idot.eq.0 .or. idot.eq.1) then
            call cdeterminant(coord,rvec_en,r_en,cvd,d2d,div_vd,psid)
          elseif(idot.eq.3) then
            call cdeterminant_cf(coord,rvec_en,r_en,cvd,d2d,div_vd,psid)
          else
            d2d=0.d0
            psid=1.d0
            do 3 i=1,nelec
              div_vd(i)=0.d0
              do 3 k=1,ndim
    3           cvd(k,i)=(0.d0,0.d0)
          endif
          if(idot.gt.0) then
            call cjastrow(cvk,d2k,div_vk,psik)   ! this is not a variational jastrow factor
          endif

         else

          d2d=d2d_sav(iconfg)
          psid=psid_sav(iconfg)
          do 7 i=1,nelec
            div_vd(i)=div_vd_sav(i,iconfg)
            do 7 k=1,ndim
    7         cvd(k,i)=cvd_sav(k,i,iconfg)
          if(idot.gt.0) then
            psik=psik_sav(iconfg)
            d2k=d2k_sav(iconfg)
            do 8 i=1,nelec
              div_vk(i)=div_vk_sav(i,iconfg)
              do 8 k=1,ndim
    8           cvk(k,i)=cvk_sav(k,i,iconfg)
          endif
        endif

c calculate jastrow and determinantal components
c note : div_v is not used in this routine and should be ignored also in /fitdet/
        do 10 i=1,nelec
          do 10 k=1,ndim
            if(ipr.ge.2) write(6,'(''cvd,vj'',9d12.5)') cvd(k,i),vj(k,i)
            cvelocity(k,i)=vj(k,i)+cvd(k,i)
   10   continue
        d2psi=d2d+d2j

c calculate extra pieces due to composite fermion vortices.
        if(idot.gt.0) then
          psid=psid*psik
          d2psi=d2psi+d2k
          do 12 i=1,nelec
            do 12 k=1,ndim
              cvelocity(k,i)=cvelocity(k,i)+cvk(k,i)
c    complex jastrow factor complicates calculation of jackson feenberg energy:
              temp=2*dimag(cvk(k,i))*(dimag(cvk(k,i))+2*dimag(cvd(k,i)))
              d2psi=d2psi-temp
   12     continue
        endif

c calculate velocity..
        velocity2=0.d0
        do 13 i=1,nelec
          do 13 k=1,ndim
            velocity(k,i)=dreal(cvelocity(k,i))
            velocity2=velocity2+dreal(cvelocity(k,i))**2+dimag(cvelocity(k,i))**2
   13   continue


        if(ipr.ge.2) write(6,'(''coord'',9d12.5)') ((coord(k,ielec),k=1,ndim),ielec=1,nelec)
        if(ipr.ge.2) write(6,'(''psid,psij'',9d12.5)') psid,psij
        if(ipr.ge.2) write(6,'(''d2d,d2j'',9d12.5)') d2d,d2j

c calculate local energy. for convenience "magnetic energy" is included in pe.
        pe=pe+emag
        energy=pe-hb*(velocity2+d2psi)

c calculate parameter-derivatives of csf_coefs for optimization
        do 35 iparm=1,nparmcsf
          stop 'optimization of csf parameters not possible with complex wfs yet'
c        denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)
c          if(nloc.eq.0) then
c            denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)
c          else
c            denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)+dpe(iparm)
c          endif
c          do 35 i=1,nelec
c            do 35 k=1,ndim
c   35         denergy(iparm)=denergy(iparm)-hb*(2*(ddeti_det(k,i,iparm)-deti_det(iparm)*vd(k,i))*vj(k,i))
   35   continue

c calculate parameter-derivatives of Jastrow factor for optimization
c if Psi=J*D, J=exp(f), V=grad(Psi)/Psi, E_L = H*Psi/Psi and _i denotes deriv wrt. parameter c_i, then
c E_L,i = -0.5*(\nabla^2*f_i + 2 grad(f_i).V) + (\hat{V}\psi/psi)_i
        do 40 iparm=1,nparmj+nparms
          if(nloc.eq.0) then
            denergy(nparmcsf+iparm)=-hb*d2g(iparm)
          else
            denergy(nparmcsf+iparm)=-hb*d2g(iparm)+dpe(nparmcsf+iparm)
          endif
          do 40 i=1,nelec
            do 40 k=1,ndim
   40         denergy(nparmcsf+iparm)=denergy(nparmcsf+iparm)-2*hb*(g(k,i,iparm)*velocity(k,i))
c    &    +g(2,i,iparm)*velocity(2,i)+g(3,i,iparm)*velocity(3,i))

      else

        if(idodet.eq.1) then
          call determinant(coord,rvec_en,r_en,vd,d2d,div_vd,psid)
         else
          d2d=d2d_sav(iconfg)
          psid=psid_sav(iconfg)
          do 45 i=1,nelec
            div_vd(i)=div_vd_sav(i,iconfg)
            do 45 k=1,ndim
   45         vd(k,i)=vd_sav(k,i,iconfg)
        endif

c combine to get derivatives of full wavefunction
c Gradient and Laplacian of ln(JD) are the sum of the gradient and Laplacian of ln(J) and ln(D).
        do 50 i=1,nelec
c         div_v(i)=div_vj(i)+div_vd(i)
          do 50 k=1,ndim
            if(ipr.ge.2) write(6,'(''ifr,vd,vj'',i2,9d12.5)') ifr,vd(k,i),vj(k,i)
   50       velocity(k,i)=vj(k,i)+vd(k,i)
        d2psi=d2d+d2j

        if(ipr.ge.2) write(6,'(''coord'',9d12.5)') ((coord(k,ielec),k=1,ndim),ielec=1,nelec)
        if(ipr.ge.2) write(6,'(''ifr,psid,psij'',i2,9d12.5)') ifr,psid,psij
        if(ipr.ge.2) write(6,'(''ifr,d2d,d2j'',i2,9d12.5)') ifr,d2d,d2j

c get pseudo-potential contribution
c nonloc_pot must be called after determinant because psid is needed in nonloc_pot
c This call computes the entire pot_en if iperiodic=0 but only the nonlocal contrib. otherwise.
        if(ipr.ge.3) write(6,'(''deriv_hpsi: pe before nonloc_pot'',9f12.5)') pe
        if(nloc.gt.0) then
          call deriv_nonloc_pot(coord,rshift,rvec_en,r_en,detu,detd,deti_det,slmui,slmdi,vpsp,psid,pe,dpe,ifr)
        else
c Next 2 lines not needed since we have if(nloc.gt.0) when eval. denergy(iparm)
          do 55 iparm=1,nparmcsf+nparmj+nparms
   55       dpe(iparm)=0
        endif

        if(ipr.ge.3) write(6,'(''pe after nonloc_pot'',9f12.5)') pe

c calculate local energy. for convenience "magnetic energy" is included in pe.
        pe=pe+emag

c Laplacian(JD)/JD = Laplacian(ln(JD)) + V^2
        energy=-hb*d2psi+pe
        do 60 i=1,nelec
          do 60 idim=1,ndim
   60       energy=energy-hb*velocity(idim,i)**2
        if(ipr.ge.3) write(6,'(''ifr,energy='',i2,9f16.10)') ifr,energy,psid,psij,pe

c calculate parameter-derivatives of csf_coefs for optimization
        do 65 iparm=1,nparmcsf
c       denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)
          if(nloc.eq.0) then
            denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)
          else
            denergy(iparm)=-hb*(d2deti_det(iparm)-deti_det(iparm)*d2det_det)+dpe(iparm)
          endif
          do 65 i=1,nelec
            do 65 k=1,ndim
   65         denergy(iparm)=denergy(iparm)-hb*(2*(ddeti_det(k,i,iparm)-deti_det(iparm)*vd(k,i))*vj(k,i))

c calculate parameter-derivatives of Jastrow factor for optimization
c if Psi=J*D, J=exp(f), V=grad(Psi)/Psi, E_L = H*Psi/Psi and _i denotes deriv wrt. parameter c_i, then
c E_L,i = -0.5*(\nabla^2*f_i + 2 grad(f_i).V) + (\hat{V}\psi/psi)_i
        do 70 iparm=1,nparmj+nparms
          if(nloc.eq.0) then
            denergy(nparmcsf+iparm)=-hb*d2g(iparm)
          else
            denergy(nparmcsf+iparm)=-hb*d2g(iparm)+dpe(nparmcsf+iparm)
          endif
          do 70 i=1,nelec
            do 70 k=1,ndim
   70         denergy(nparmcsf+iparm)=denergy(nparmcsf+iparm)-2*hb*(g(k,i,iparm)*velocity(k,i))
c    &    +g(2,i,iparm)*velocity(2,i)+g(3,i,iparm)*velocity(3,i))

      endif

      return
      end

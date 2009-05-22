      subroutine deriv_nonloc_pot(x,rshift,rvec_en,r_en,detu,detd,deti_det,slmui,slmdi,vpsp,psid,pe,dpe,ifr)
c Written by Claudia Filippi; modified by Cyrus Umrigar
c Calculates non-local potential and its derivatives wrt Jastrow parameters
c V_1=A/B, where V_1 is the nonlocal part of the potential on electron 1
c A=(1/N_quad) \sum_l (2l+1) (V_l-V_L) \sum_j^{N_quad} P_l(cos \theta_j) Psi(r_1j,...)
c the sum on j is over the N_quad quadrature points.
c B = Psi = JD
c V' = A'/B - A*B'/B^2
c The arrays in the routine are:
c vpot = (1/N_quad) \sum_l (2l+1) (V_l-V_L) \sum_j^{N_quad} P_l(cos \theta_j) Psi(r_j,...)
c vps  = (V_l-V_L)
      use control_mod
      use deriv_orb_mod
      use eloc_mod
      use periodic_jastrow_mod  !WAS
      use atom_mod
      use dets_mod
      use optim_mod
      use const_mod
      use pseudo_mod
      implicit real*8(a-h,o-z)

      common /contrl_per/ iperiodic,ibasis
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
!JT     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase

      common /derivjas/ gvalue(MPARMJ),g(3,MELEC,MPARMJ),d2g(MPARMJ)
     &,go(MELEC,MELEC,MPARMJ)

      dimension x(3,*),rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &,detu(*),detd(*),deti_det(MPARMD),slmui(nupdn_square,*),slmdi(nupdn_square,*)
     &,dvpsp(MPARM),dpe(MPARM)

      do 20 i=1,nelec
        if(nloc.eq.1) then
          call getvps_fahy(r_en,i)
         elseif(nloc.ge.2 .and. nloc.le.5) then
          call getvps_champ(r_en,i)
         elseif(nloc.eq.6) then
          call getvps_gauss(r_en,i)
         else
          stop "nloc is not implemented in deriv_nonloc_pot.f"
        endif
   20 continue

c local component
      if(iperiodic.eq.0) then
        do 30 ic=1,ncent
          do 30 i=1,nelec
   30       pe=pe+vps(i,ic,lpotp1(iwctype(ic)))
      endif

!     total ee and local en potential
      eloc_pot_loc = pe                      !JT
      call object_modified_by_index (eloc_pot_loc_index)  !JT

c non-local component and its derivative (division by the Jastrow already in nonloc)
      call deriv_nonloc(x,rshift,rvec_en,r_en,detu,detd,slmui,slmdi,vpsp,dvpsp)
      do 35 iparm=1,nparmcsf
  35    dpe(iparm)=(dvpsp(iparm)-vpsp*deti_det(iparm))/psid
      if(ipr.ge.4) write(6,'(''dvpsp2(iparm)'',30d12.4)') (dvpsp(iparm),iparm=1,nparmcsf)
      if(ipr.ge.4) write(6,'(''dpe(iparm)'',30d12.4)') (dpe(iparm),iparm=1,nparmcsf)
      if(ipr.ge.4) write(6,'(''deti_det(iparm)'',30d12.4)') (deti_det(iparm),iparm=1,nparmcsf)
      do 40 iparm=1,nparmj
  40    dpe(nparmcsf+iparm)=(dvpsp(nparmcsf+iparm)-vpsp*gvalue(iparm))/psid
      pe=pe+vpsp/psid

      eloc_pot_nloc = vpsp/psid
      call object_modified_by_index (eloc_pot_nloc_index)  !JT

!WAS
      psid_pjas = psid
      call object_modified_by_index (psid_pjas_index)
!WAS

      if(ipr.ge.4) write(6,'(''pe,vpsp/psid,vpsp,psid,detu(1),detd(1)2='',2f9.4,9d12.4)')
     &pe,vpsp/psid,vpsp,psid,detu(1),detd(1),r_en(1,1)

      return
      end

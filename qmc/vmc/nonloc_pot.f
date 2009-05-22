      subroutine nonloc_pot(x,rshift,rvec_en,r_en,detu,detd,slmui,slmdi,vpsp,psid,pe)
c Written by Claudia Filippi, modified by Cyrus Umrigar
c Calculates the local and nonlocal components of the pseudopotential
c pe_en(loc) is computed in distances and pe_en(nonloc) here in nonloc_pot if nloc !=0 and iperiodic!=0.

      use control_mod
      use deriv_orb_mod
      use eloc_mod
      use periodic_jastrow_mod  !WAS
      use atom_mod
      use dets_mod

      use const_mod
      use pseudo_mod
      implicit real*8(a-h,o-z)

      common /contrl_per/ iperiodic,ibasis
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
!JT     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc

      dimension x(3,*),rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &,detu(*),detd(*),slmui(nupdn_square,*),slmdi(nupdn_square,*)

c Calculate local and nonlocal pseudopotentials for all electrons, nuclei and l-components
c and store in vps
      do 20 i=1,nelec
        if(nloc.eq.1) then
          call getvps_fahy(r_en,i)
         elseif(nloc.ge.2 .and. nloc.le.5) then
          call getvps_champ(r_en,i)
         elseif(nloc.eq.6) then
          call getvps_gauss(r_en,i)
         elseif(nloc.ge.7) then
          stop 'nloc >= 7'
        endif
   20 continue

      if(ipr.ge.4) write(6,'(''r_en='',9f9.5)') ((r_en(iel,ic),iel=1,nelec),ic=1,ncent)

c local component
      if(iperiodic.eq.0) then
        if(ipr.ge.4) write(6,'(''nonloc_pot: pe before local psp'',f9.5)') pe
        do 30 ic=1,ncent
          do 30 i=1,nelec
            if(ipr.ge.4) write(6,'(''i,ic,vps,pe'',2i3,9f9.5)')i,ic,vps(i,ic,lpotp1(iwctype(ic))),pe+vps(i,ic,lpotp1(iwctype(ic)))
   30       pe=pe+vps(i,ic,lpotp1(iwctype(ic)))
        if(ipr.ge.4) write(6,'(''nonloc_pot: pe after local psp'',f12.5)') pe
      endif

!     Total ee and local en potential
      eloc_pot_loc = pe                      !JT
      call object_modified_by_index (eloc_pot_loc_index)  !JT

c non-local component (division by the Jastrow already in nonloc)
      call nonloc(x,rshift,rvec_en,r_en,detu,detd,slmui,slmdi,vpsp)
      pe=pe+vpsp/psid
c     write(6,'(''pe='',9d12.5)') pe,vpsp,psid,detu(1),detd(1)
      if(ipr.ge.4) write(6,'(''pe,vpsp/psid,vpsp,psid,detu(1),detd(1),r_en(1,1)='',2f9.4,9d12.4)')
     &pe,vpsp/psid,vpsp,psid,detu(1),detd(1),r_en(1,1)

      return
      end

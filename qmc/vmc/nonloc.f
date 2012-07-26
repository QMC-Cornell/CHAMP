      subroutine nonloc(x,rshift,rvec_en,r_en,detu,detd,slmui,slmdi,vpsp)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use constants_mod
      use control_mod
      use deriv_orb_mod
      use periodic_jastrow_mod !WAS
      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use pseudo_mod
      use contrl_per_mod
      use periodic_mod
      use qua_mod
      use contrldmc_mod, only : taunow
      use config_dmc_mod, only : xoldw, voldw, psidow, psijow
      use average_mod, only : current_walker
c Temporary:
      use div_v_dmc_mod
      use delocc_mod

      implicit real*8(a-h,o-z)

      dimension x(3,*),rshift(3,nelec,*),rvec_en(3,nelec,*),r_en(nelec,*)
     &,detu(*),detd(*),slmui(nupdn_square,*),slmdi(nupdn_square,*)
      dimension rr_en(nelec,ncent),rr_en2(nelec,ncent),rr_en_sav(ncent),rr_en2_sav(ncent)
     &,xsav(3),rshift_sav(3,ncent),rvec_en_sav(3,ncent),r_en_sav(ncent),vpot(MPS_L)
      dimension xnew(3),vnew(3,nelec) ! local variables for tmoves
      dimension x_tmove_sav(3,nquad*ncent), vpsp_tmove(nquad*ncent), vpsp_tmove_heatbath(nquad*ncent)

      do 10 i=1,nelec
        do 10 ic=1,ncent
          call scale_dist(r_en(i,ic),rr_en(i,ic),1)
   10     call scale_dist(r_en(i,ic),rr_en2(i,ic),3)

c     write(6,'(''l_do_tmoves,x='',l2,30f9.4)') l_do_tmoves, ((x(k,i),k=1,ndim),i=1,nelec)
c     write(6,'(''l_do_tmoves, r_en='',l2,30f9.4)') l_do_tmoves, ((r_en(i,ic),i=1,nelec),ic=1,ncent)
c     write(6,'(''rvec_en='',60f9.4)') (((rvec_en(k,i,ic),k=1,ndim),i=1,nelec),ic=1,ncent)

      if (l_opt_orb_energy) then
        call object_provide_by_index (param_orb_nb_index)
        call object_alloc ('vpot_ex', vpot_ex, MPS_L, param_orb_nb)
        call object_alloc ('vpsp_ex', vpsp_ex, param_orb_nb)
        vpsp_ex = 0.d0
      endif

c     do 3 i=1,nelec
c       do 3 j=1,i-1
c         if(abs(r_en(i,1)-r_en(j,1)).lt.1.d-6) write(6,'(''1 current_walker,i,j,r_en(i,1),r_en(j,1)'',3i3,9d16.8)')
c    &    current_walker,i,j,r_en(i,1),r_en(j,1)
c   3 continue

      vpsp=0
      do 150 i=1,nelec

c Save position ith electron and its distances etc. from all nuclei
        do 11 k=1,ndim
   11     xsav(k)=x(k,i)
        do 12 jc=1,ncent
          r_en_sav(jc)=r_en(i,jc)
          rr_en_sav(jc)=rr_en(i,jc)
          rr_en2_sav(jc)=rr_en2(i,jc)
          do 12 k=1,ndim
            rshift_sav(k,jc)=rshift(k,i,jc)
   12           rvec_en_sav(k,jc)=rvec_en(k,i,jc)

        ntmove_pts=0
        do 100 ic=1,ncent
          ict=iwctype(ic)

c vps was calculated by calling getvps_xx from nonloc_pot
          iskip=1
          do 15 l=1,npotd(ict)
   15       if(l.ne.lpotp1(ict) .and. dabs(vps(i,ic,l)).gt.1.d-4) iskip=0

          if(iskip.eq.0) then

            if(l_do_tmoves) call rotqua ! If we are doing t-moves then the grid needs to be rotated more often to prevent electrons from being along rays

            ri=one/r_en(i,ic)

            do 20 l=1,npotd(ict)
   20         if(l.ne.lpotp1(ict)) vpot(l)=0

            if (l_opt_orb_energy) then  !JT
               vpot_ex = 0.d0           !JT
            endif                       !JT

            do 60 iq=1,nquad
              ntmove_pts=ntmove_pts+1
              costh=rvec_en_sav(1,ic)*xq(iq)+rvec_en_sav(2,ic)*yq(iq)+rvec_en_sav(3,ic)*zq(iq)
              costh=costh*ri

              if(iperiodic.eq.0) then
                x(1,i)=r_en(i,ic)*xq(iq)+cent(1,ic)
                x(2,i)=r_en(i,ic)*yq(iq)+cent(2,ic)
                x(3,i)=r_en(i,ic)*zq(iq)+cent(3,ic)
              else
                x(1,i)=r_en(i,ic)*xq(iq)+cent(1,ic)+rshift(1,i,ic)
                x(2,i)=r_en(i,ic)*yq(iq)+cent(2,ic)+rshift(2,i,ic)
                x(3,i)=r_en(i,ic)*zq(iq)+cent(3,ic)+rshift(3,i,ic)
              endif

c Since we are rotating on sphere around nucleus ic, that elec-nucl distance does not change but distances to other nuclei do
              do 40 jc=1,ncent
                do 38 k=1,ndim
   38             rvec_en(k,i,jc)=x(k,i)-cent(k,jc)

                if(jc.ne.ic) then
                  if(iperiodic.eq.0) then
                    r_en(i,jc)=0
                    do 39 k=1,ndim
   39                 r_en(i,jc)=r_en(i,jc)+rvec_en(k,i,jc)**2
                    r_en(i,jc)=dsqrt(r_en(i,jc))
                   else
                    call find_image4(rshift(1,i,jc),rvec_en(1,i,jc),r_en(i,jc),rlatt,rlatt_inv)
                  endif

                  call scale_dist(r_en(i,jc),rr_en(i,jc),1)
                  call scale_dist(r_en(i,jc),rr_en2(i,jc),3)
                endif
   40         continue

              iel=i

              electron = iel !JT
              call object_modified_by_index (electron_index) !JT

              call nonlocd(iel,x(1,i),rvec_en,r_en,detu,detd,slmui,slmdi,deter)
              call nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,value)

              if (do_pjas) then ! periodic Jastrow implemented by WAS
                 call nonloc_pjas (iel, x(:,1:nelec), value)
              endif

              if(ipr.ge.4) then
                write(6,'(''rr_en,rr_en2'',2d14.6)') rr_en(1,1),rr_en2(1,1)
                write(6,'(''ic,i,iq,deter,value'',3i3,2d14.6)') ic,i,iq,deter,value
              endif

              if(current_walker.eq.0) current_walker=1
              if(ipr.ge.1)
     &          write(6,'(''outside tmoves: vps(i,ic,l),wq(iq),yl0(l,costh),deter,psidow(current_walker,1),exp(value),taunow'',
     &          9d12.4)') vps(i,ic,l),wq(iq),yl0(l,costh),deter,psidow(current_walker,1),exp(value),taunow

              if(l_do_tmoves) then
                vpsp_tmove(ntmove_pts)=0
                x_tmove_sav(1:3,ntmove_pts)=x(1:3,i)
                if(ipr.ge.1) write(6,'(''iw,ielec='',9i5)') current_walker, i
              endif
              do 50 l=1,npotd(ict)
                if(l.ne.lpotp1(ict)) then
                  if(l_do_tmoves) then
                    if(ipr.ge.1)
     &              write(6,'(''vps(i,ic,l),wq(iq),yl0(l,costh),deter,psidow(current_walker,1),exp(value),taunow'',9d12.4)')
     &              vps(i,ic,l),wq(iq),yl0(l,costh),deter,psidow(current_walker,1),exp(value),taunow
                    vpsp_tmove(ntmove_pts)=vpsp_tmove(ntmove_pts)
     &              -vps(i,ic,l)*wq(iq)*yl0(l,costh)*(deter/psidow(current_walker,1))*exp(value)*taunow
                    call systemflush(6)
                  endif
                  vpot(l)=vpot(l)+wq(iq)*yl0(l,costh)*deter*exp(value)
                  if(ipr.ge.1) write(6,'(''l,yl0(l,costh),deter,exp(value),yl0(l,costh)*deter*exp(value),vpot(l)'',i3,9f20.15)')
     &            l,yl0(l,costh),deter,exp(value),yl0(l,costh)*deter*exp(value),vpot(l)

! JT              For singly-excited wave functions
                  if (l_opt_orb_energy) then
                     call object_provide_by_index (psid_ex_in_x_index)
                     do iex = 1, param_orb_nb
                      vpot_ex(l,iex)=vpot_ex(l,iex)+wq(iq)*yl0(l,costh)*psid_ex_in_x(iex)*exp(value)
                     enddo
                  endif

                endif
   50         continue ! npotd(ict)

   60       continue ! nquad

            do 68 k=1,ndim
   68         x(k,i)=xsav(k)
            do 70 jc=1,ncent
              r_en(i,jc)=r_en_sav(jc)
              rr_en(i,jc)=rr_en_sav(jc)
              rr_en2(i,jc)=rr_en2_sav(jc)
              do 70 k=1,ndim
                rshift(k,i,jc)=rshift_sav(k,jc)
   70           rvec_en(k,i,jc)=rvec_en_sav(k,jc)

            do 80 l=1,npotd(ict)
              if(l.ne.lpotp1(ict)) then
                vpsp=vpsp+vps(i,ic,l)*vpot(l)
                if(ipr.ge.4) write(6,'(''nonloc: i,ic,l,vps(i,ic,l),vpot(l),vpsp'',3i5,9d12.4)') i,ic,l,vps(i,ic,l),vpot(l),vpsp

! JT            For singly-excited wave functions
                if (l_opt_orb_energy) then
                   do iex = 1, param_orb_nb
                    vpsp_ex(iex)=vpsp_ex(iex)+vps(i,ic,l)*vpot_ex(l,iex)
                   enddo
                endif

              endif
   80       continue

          endif
  100   continue ! ncent

c       if(l_do_tmoves) write(6,'(''i, current_walker, x_tmove_sav='',2i3, 99d12.4)') i, current_walker,
c    &  ((x_tmove_sav(k,ii),k=1,3),ii=1,ntmove_pts)

c Warning: For the moment correlated sampling for forces does not work with t-moves.
c Decide which if any tmove to perform for electron i
c Note we are using the same array, vpsp_tmove, to store the cumulative values to be used for the heat bath.
        if(l_do_tmoves .and. ntmove_pts.ne.0) then
          vpsp_tmove_sum=1
          do itmove_pts=1,ntmove_pts
            vpsp_tmove_sum=vpsp_tmove_sum+max(vpsp_tmove(itmove_pts),0.d0)
            if(itmove_pts.eq.1) then
              vpsp_tmove_heatbath(itmove_pts)=1
            else
              vpsp_tmove_heatbath(itmove_pts)=vpsp_tmove_heatbath(itmove_pts-1)+max(vpsp_tmove(itmove_pts-1),0.d0)
            endif
          enddo
          if(ipr.ge.1)
     &    write(6,'(''vpsp_tmove_sum, vpsp_tmove_heatbath'',100f10.6)') vpsp_tmove_sum, (vpsp_tmove_heatbath(ii),ii=1,ntmove_pts)

          if(abs(vpsp_tmove_sum-1).ge.1.d-6) random=rannyu(0)
          if(vpsp_tmove_heatbath(1).gt.random*vpsp_tmove_sum) then ! No t-move (vpsp_tmove_heatbath(1)=1)
            iwhich_tmove=0
          else
            iwhich_tmove=ntmove_pts ! Last quadrature point is chosen if we make it through loop without satisfying the if within the loop.
            do itmove_pts=2,ntmove_pts
              if(vpsp_tmove_heatbath(itmove_pts).gt.random*vpsp_tmove_sum) then
                iwhich_tmove=itmove_pts-1
                exit
              endif
            enddo
!           x(1:3,i)=x_tmove_sav(1:3,iwhich_tmove)
            xnew(1:3)=x_tmove_sav(1:3,iwhich_tmove)
            xoldw(1:3,i,current_walker,1)=x_tmove_sav(1:3,iwhich_tmove)
c           write(6,'(''1 xoldw='',99d12.4)') ((xoldw(k,ii,current_walker,1),k=1,3),ii=1,nelec)
c           call systemflush(6)
!           call hpsiedmc(i,iw,xnew,psidn,psijn,vnew)
            call hpsiedmc(i,current_walker,xnew,psidow(current_walker,1),psijow(current_walker,1),voldw(1,1,current_walker,1))
c           write(6,'(''2 xoldw='',99d12.4)') ((xoldw(k,ii,current_walker,1),k=1,3),ii=1,nelec)
c           write(6,'(''2 xnew='',99d12.4)') (xnew(k),k=1,3)
            call systemflush(6)
            call jassav(i)
            if(ibasis.eq.3) then                        ! complex calculations
              call cdetsav(i)
            else
              call detsav(i)
            endif

c           write(6,'(''3 xoldw='',99d12.4)') ((xoldw(k,ii,current_walker,1),k=1,3),ii=1,nelec)
c           if(l_do_tmoves) write(6,'(''3 xnew='',99d12.4)') (xnew(k),k=1,3)
          endif

        endif ! l_do_tmoves
        call systemflush(6)

  150 continue ! nelec

      if (l_do_tmoves) then !FP
c       if(ifr.eq.1) then
c Primary configuration
c         drifdifr=one
c         if(nforce.gt.1) call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)

c         call hpsi(xnew,psidow(current_walker,1),psijow(current_walker,1),voldw(1,1,current_walker,1),
c    &    div_vow(1,current_walker),d2n,pen,pein,enew,denergy,1)

          psi_det = psidow(current_walker,1)
          psi_jas = exp(psijow(current_walker,1))
          call object_modified_by_index (voldw_index)
          call object_modified_by_index (psi_det_index)
          call object_modified_by_index (psi_jas_index)

          if(ibasis.eq.3) then                     !complex basis set
            call cwalksav_det(current_walker)
          else
            call walksav_det(current_walker)
          endif
          call walksav_jas(current_walker)
c       else
c Warning to be done
c       endif
      endif

      call object_modified_by_index (vpsp_ex_index) ! JT

c     write(6,'(''x='',30f9.4)') ((x(k,i),k=1,ndim),i=1,nelec)
c     write(6,'(''r_en='',30f9.4)') ((r_en(i,ic),i=1,nelec),ic=1,ncent)
c     write(6,'(''rvec_en='',60f9.4)') (((rvec_en(k,i,ic),k=1,ndim),i=1,nelec),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------

      function yl0(l,costh)
c (2L+1)*P_L(costh)
c This is not quite Y_L0 but sqrt(4pi/(2L+1)) Y_L0
c Note that the associated P_L^m and the unassociated P_L Legendre polynomials are the same for m=0.
c l is actually L+1.

      implicit real*8(a-h,o-z)

      if(l.eq.1) then
        yl0=1.d0
       elseif(l.eq.2) then
        yl0=3.d0*costh
       elseif(l.eq.3) then
        yl0=2.5d0*(3*costh*costh-1)
       elseif(l.eq.4) then
        yl0=3.5d0*costh*(5*costh*costh-3)
       elseif(l.eq.5) then
        yl0=1.125d0*(35*costh**4-30*costh**2+3)
       else
        stop 'yl0 implemented to l=4 only (Warning: l is l+1)'
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine nonlocd(iel,x,rvec_en,r_en,detu,detd,slmui,slmdi,determ)
c Written by Claudia Filippi, modified by Cyrus Umrigar
      use all_tools_mod
      use control_mod
      use eloc_mod
      use dorb_mod
      use slatn_mod
      use orbe_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use contr2_mod
      use contrl_opt2_mod
      use wfsec_mod
      use contrl_per_mod
      use contr3_mod
      use phifun_mod
      use const_mod
      use slatn2_mod
      implicit real*8(a-h,o-z)

      dimension x(3),rvec_en(3,nelec,*),r_en(nelec,*)
     &,detu(*),detd(*),slmui(nupdn_square,*),slmdi(nupdn_square,*)
      dimension ratio(ndet)

c     determ=0

c get orbitals for all electron iel
c Note x is 3-dim but rvec_en,r_en have info about all electrons.
c So, iel is passed to select elements of rvec_en,r_en,phin and for IO.
      if(iperiodic.eq.0) then

        if(inum_orb.eq.0) then
          call orbitals_loc_anae(iel,rvec_en,r_en,orbe)
         else
          call orbitals_loc_nume(x,orbe)
        endif

       else

        if(inum_orb.eq.0) then
          call orbitals_pwe(iel,x,orbe)
         else
          call orbitals_period_nume(x,orbe)
        endif

      endif

      call object_modified_by_index (orbe_index) !JT

      if(iel.le.nup) then

        ikel=nup*(iel-1)

        do idet=1,ndetup
           ratio(idet)=0
           do j=1,nup
              ratio(idet)=ratio(idet)+slmui(j+ikel,idet)*orbe(iworbdup(j,idet))
           enddo
        enddo

        if(.not. l_opt_exp) then
           do idet=1,ndetup
              detn(idet)=detu(idet)*ratio(idet)
           enddo
        else
           do idet=1,ndetup
              detn(idet)=detu(idet)*ratio(idet)
              do i=1,nup
                 if(i.ne.iel) then
                    ik=nup*(i-1)
                    sum=0
                    do j=1,nup
                       sum=sum+slmui(j+ik,idet)*orbe(iworbdup(j,idet))
                    enddo
                    sum=sum/ratio(idet)
                    do j=1,nup
                       slmin(j+ik,idet)=slmui(j+ik,idet)-slmui(j+ikel,idet)*sum
                    enddo
                 endif
              enddo
              do j=1,nup
                 slmin(j+ikel,idet)=slmui(j+ikel,idet)/ratio(idet)
              enddo
           enddo
        endif

      else

        ikel=ndn*(iel-nup-1)

        do idet=1,ndetdn
           ratio(idet)=0
           do j=1,ndn
              ratio(idet)=ratio(idet)+slmdi(j+ikel,idet)*orbe(iworbddn(j,idet))
           enddo
        enddo
        if(.not. l_opt_exp) then
           do idet=1,ndetdn
              detn(idet)=detd(idet)*ratio(idet)
           enddo
        else
           do idet=1,ndetdn
              detn(idet)=detd(idet)*ratio(idet)
              do i=1,ndn
                 if(i+nup.ne.iel) then
                    ik=ndn*(i-1)
                    sum=0
                    do j=1,ndn
                       sum=sum+slmdi(j+ik,idet)*orbe(iworbddn(j,idet))
                    enddo
                    sum=sum/ratio(idet)
                    do j=1,ndn
                       slmin(j+ik,idet)=slmdi(j+ik,idet)-slmdi(j+ikel,idet)*sum
                    enddo
                 endif
              enddo
              do j=1,ndn
                 slmin(j+ikel,idet)=slmdi(j+ikel,idet)/ratio(idet)
              enddo
           enddo
        endif

      endif

      call object_modified_by_index (detn_index) !JT
      call object_modified_by_index (slmin_index) !fp

      determ=0
      do 115 icsf=1,ncsf
        do 115 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
          if(iel.le.nup) then
            term=detn(iwdetup(idet))*detd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
           else
            term=detu(iwdetup(idet))*detn(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
          endif
  115     determ=determ+term

c Derivatives wrt to csf_coefs for optimizing them
      if(index(mode,'fit').ne.0 .or. igradhess.ge.1 .or. l_opt_csf) then
        do 140 iparm=1,nparmcsf
          icsf=iwcsf(iparm)
          deti_new(iparm)=0
          do 140 idet_in_csf=1,ndet_in_csf(icsf)
            idet=iwdet_in_csf(idet_in_csf,icsf)
            if(iel.le.nup) then
              term=detn(iwdetup(idet))*detd(iwdetdn(idet))*cdet_in_csf(idet_in_csf,icsf)
             else
              term=detu(iwdetup(idet))*detn(iwdetdn(idet))*cdet_in_csf(idet_in_csf,icsf)
            endif
  140       deti_new(iparm)=deti_new(iparm)+term
      endif

      return
      end
c-----------------------------------------------------------------------

!WAS  subroutine nonlocj(iel,x,rshift,rr_en,rr_en2,value)
      subroutine nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,value)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use control_mod
      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use contrl_per_mod
      use jaspar_mod
      use bparm_mod
      use periodic_mod
      use jaso_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*),rshift(3,nelec,ncent),rr_en(nelec,ncent),rr_en2(nelec,ncent),fsn(nelec,nelec),dx(3)

      dimension r_en(nelec,ncent)

      fsumn=0

      if(nelec.lt.2) goto 47

      do 45 jj=1,nelec

        if(jj.eq.iel) goto 45
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif

        sspinn=1
        ipar=0
        isb=1
        if(i.le.nup .or. j.gt.nup) then
          if(nspin2b.eq.2) then
            isb=2
           elseif(nocuspb.eq.0) then
            sspinn=half
          endif
          ipar=1
        endif

        do 10 k=1,ndim
   10     dx(k)=x(k,jj)-x(k,iel)

        if(iperiodic.eq.0) then
          rij=0
          do 20 k=1,ndim
   20       rij=rij+dx(k)**2
          rij=dsqrt(rij)
         else
          call find_image3(dx,rij,rlatt_sim,rlatt_sim_inv)
        endif

c e-e terms
        call scale_dist(rij,u,2)

        fsn(i,j)=psibnl(u,isb,ipar)

c e-e-n terms
c The scaling is switched in psinl, so do not do it here.
c     if(isc.ge.12) call scale_dist(rij,u,3)
      call scale_dist(rij,u,4)

        do 40 ic=1,ncent
          it=iwctype(ic)
!WAS 40   fsn(i,j)=fsn(i,j) + psinl(u,rshift(1,i,ic),rshift(1,j,ic),rr_en2(i,ic),rr_en2(j,ic),it)
   40     fsn(i,j)=fsn(i,j) + psinl(u,rshift(1,i,ic),rshift(1,j,ic),r_en(i,ic),r_en(j,ic),rr_en2(i,ic),rr_en2(j,ic),it)
        fsumn=fsumn+fsn(i,j)-fso(i,j)

   45 continue

c e-n terms
   47 fsn(iel,iel)=0
      do 50 ic=1,ncent
        it=iwctype(ic)
   50   fsn(iel,iel)=fsn(iel,iel)+psianl(rr_en(iel,ic),it)

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)

      value=fsumn

      return
      end

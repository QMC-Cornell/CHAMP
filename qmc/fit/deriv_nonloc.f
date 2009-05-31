      subroutine deriv_nonloc(x,rshift,rvec_en,r_en,detu,detd,slmui,slmdi,vpsp,dvpsp)
c Written by Claudia Filippi, modified by Cyrus Umrigar
      use all_tools_mod
      use constants_mod
      use control_mod
      use deriv_orb_mod
      use periodic_jastrow_mod  !WAS
      use atom_mod
      use dets_mod
      use optim_mod
      use const_mod
      use dim_mod
      use pseudo_mod
      use contrl_per_mod
      use periodic_mod
      use qua_mod
      use slatn2_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*),rshift(3,nelec,ncent),rvec_en(3,nelec,ncent),r_en(nelec,ncent)
     &,detu(*),detd(*),slmui(nupdn_square,*),slmdi(nupdn_square,*)
      dimension rr_en(nelec,ncent),rr_en2(nelec,ncent),rr_en_sav(ncent),rr_en2_sav(ncent)
     &,xsav(3),rshift_sav(3,ncent),rvec_en_sav(3,ncent),r_en_sav(ncent)
     &,vpot(MPS_L),dvpot(MPS_L,MPARM),dvpsp(MPARM),gn(MPARMJ)
      dimension dk_en(nelec,ncent),dk_en2(nelec,ncent),dk_en_sav(ncent),dk_en2_sav(ncent)

      do 10 i=1,nelec
        do 10 ic=1,ncent
          call scale_dist(r_en(i,ic),rr_en(i,ic),1)
          call scale_dist(r_en(i,ic),rr_en2(i,ic),3)
c note: dk2,dr and dr2 are unused variables here
          if(nparms.eq.1) then
            call deriv_scale(r_en(i,ic),dk_en(i,ic),dk2,dr,dr2,1,0)
            call deriv_scale(r_en(i,ic),dk_en2(i,ic),dk2,dr,dr2,3,0)
          endif
   10 continue

! JT beg
      if (l_opt_orb_energy) then
       call object_provide ('param_orb_nb')
       call object_alloc ('vpot_ex', vpot_ex, MPS_L, param_orb_nb)
       call object_alloc ('vpsp_ex', vpsp_ex, param_orb_nb)
       vpsp_ex = 0.d0
      endif
! JT end

! WAS
      if (l_opt_pjas) then
       call object_alloc ('dvpot_pjas', dvpot_pjas, MPS_L, param_pjas_nb)
       call object_alloc ('dvpsp_pjas', dvpsp_pjas, param_pjas_nb)
       do i=1, param_pjas_nb
          dvpsp_pjas (i) = 0.d0
       enddo
      endif
! WAS

      vpsp=0
      do 11 iparm=1,nparmcsf+nparmj+nparms
   11   dvpsp(iparm)=0

      do 100 ic=1,ncent

        ict=iwctype(ic)
        do 100 i=1,nelec

c vps was calculated by calling getvps_tm from deriv_nonloc_pot
          iskip=1
          do 15 l=1,npotd(ict)
   15       if(l.ne.lpotp1(ict) .and. dabs(vps(i,ic,l)).gt.1.d-4) iskip=0

          if(iskip.eq.0) then

            ri=one/r_en(i,ic)

            do 20 l=1,npotd(ict)
              if(l.ne.lpotp1(ict)) then
                vpot(l)=0
                do 18 iparm=1,nparmcsf+nparmj+nparms
   18             dvpot(l,iparm)=0
              endif
   20       continue

            if (l_opt_orb_energy) then  !JT
               vpot_ex = 0.d0           !JT
            endif                       !JT

! WAS
            if (l_opt_pjas) then
               dvpot_pjas = 0.d0
            endif
! WAS

            do 28 k=1,ndim
   28         xsav(k)=x(k,i)
            do 30 jc=1,ncent
              r_en_sav(jc)=r_en(i,jc)
              rr_en_sav(jc)=rr_en(i,jc)
              rr_en2_sav(jc)=rr_en2(i,jc)
              if(nparms.eq.1) then
                dk_en_sav(jc)=dk_en(i,jc)
                dk_en2_sav(jc)=dk_en2(i,jc)
              endif
              do 30 k=1,ndim
                rshift_sav(k,jc)=rshift(k,i,jc)
   30           rvec_en_sav(k,jc)=rvec_en(k,i,jc)

            do 60 iq=1,nquad
              costh=rvec_en_sav(1,ic)*xq(iq)
     &             +rvec_en_sav(2,ic)*yq(iq)
     &             +rvec_en_sav(3,ic)*zq(iq)
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
                  if(nparms.eq.1) then
                    call deriv_scale(r_en(i,jc),dk_en(i,jc),dk2,dr,dr2,1,0)
                    call deriv_scale(r_en(i,jc),dk_en2(i,jc),dk2,dr,dr2,3,0)
                  endif

                endif
   40         continue

              iel=i

              electron = iel !JT
              call object_modified_by_index (electron_index) !JT

              call nonlocd(iel,x(1,i),rvec_en,r_en,detu,detd,slmui,slmdi,deter)
c             call deriv_nonlocj(iel,x,rshift,rr_en,rr_en2,dk_en,dk_en2,value,gn)
c WAS
              call deriv_nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,dk_en,dk_en2,value,gn)

cWAS
              if (l_opt_pjas) then
                 call deriv_nonloc_pjas ( iel, x(:,1:nelec), value )
              endif
cWAS

              do 50 l=1,npotd(ict)
                if(l.ne.lpotp1(ict)) then
                  ppsi_csf=wq(iq)*yl0(l,costh)*exp(value)
                  ppsi_jas=wq(iq)*yl0(l,costh)*deter*exp(value)
                  vpot(l)=vpot(l)+ppsi_jas
                  do 42 iparm=1,nparmcsf
                    if(ipr.ge.4) write(6,'(''l,iparm,ppsi_jas,deti_new(iparm)'',2i5,9d12.4)')
     &              l,iparm,ppsi_jas,deti_new(iparm)
   42               dvpot(l,iparm)=dvpot(l,iparm)+ppsi_csf*deti_new(iparm)
                  do 45 iparm=1,nparmj+nparms
   45               dvpot(l,nparmcsf+iparm)=dvpot(l,nparmcsf+iparm)+ppsi_jas*gn(iparm)

! JT beg
!             For singly-excited wave functions
              if (l_opt_orb_energy) then
                 call object_provide ('psid_ex_in_x')
                 do iex = 1, param_orb_nb
                  vpot_ex(l,iex)=vpot_ex(l,iex)+
     >                wq(iq)*yl0(l,costh)*psid_ex_in_x(iex)*exp(value)
                 enddo
              endif
! JT end

!WAS
              if (l_opt_pjas) then
                 do iex = 1, param_pjas_nb
                    dvpot_pjas(l,iex)=dvpot_pjas(l,iex)+ ppsi_jas * gn_pjas (iex)
                 enddo
              endif
!WAS

                endif
   50         continue
   60       continue

            do 68 k=1,ndim
   68         x(k,i)=xsav(k)
            do 70 jc=1,ncent
              r_en(i,jc)=r_en_sav(jc)
              rr_en(i,jc)=rr_en_sav(jc)
              rr_en2(i,jc)=rr_en2_sav(jc)
              if(nparms.eq.1) then
                dk_en(i,jc)=dk_en_sav(jc)
                dk_en2(i,jc)=dk_en2_sav(jc)
              endif

              do 70 k=1,ndim
                rshift(k,i,jc)=rshift_sav(k,jc)
   70           rvec_en(k,i,jc)=rvec_en_sav(k,jc)

            do 80 l=1,npotd(ict)
              if(l.ne.lpotp1(ict)) then
                if(ipr.ge.4) write(6,'(''i,ic,l,vps(i,ic,l),vpot(l)2'',3i5,9d12.4)') i,ic,l,vps(i,ic,l),vpot(l)
                vpsp=vpsp+vps(i,ic,l)*vpot(l)
                do 75 iparm=1,nparmcsf+nparmj+nparms
                if(ipr.ge.4) write(6,'(''i,ic,l,iparm,dvpot(l,iparm)'',4i5,9d12.4)') i,ic,l,iparm,dvpot(l,iparm)
   75             dvpsp(iparm)=dvpsp(iparm)+vps(i,ic,l)*dvpot(l,iparm)
              endif

! JT beg
!             For singly-excited wave functions
              if (l_opt_orb_energy) then
                 do iex = 1, param_orb_nb
                  vpsp_ex(iex)=vpsp_ex(iex)+vps(i,ic,l)*vpot_ex(l,iex)
                 enddo
              endif
! JT end

!     WAS
!     For periodic jastrow
              if (l_opt_pjas) then
                 do iex = 1, param_pjas_nb
                    dvpsp_pjas(iex)=dvpsp_pjas(iex)+vps(i,ic,l)*dvpot_pjas(l,iex)
                 enddo
!     write(87,'(I4,100e10.3)') l,dvpot_pjas(l,1:param_pjas_nb)
              endif
!!!

   80       continue
          if(ipr.ge.4) write(6,'(''dvpsp(iparm)'',40d12.4)') (dvpsp(iparm),iparm=1,nparmcsf+nparmj)
          endif
  100 continue

      if (l_opt_orb_energy) then
       call object_modified_by_index (vpsp_ex_index) !JT
      endif

!WAS
      if (l_opt_pjas) then
         call object_modified_by_index (dvpsp_pjas_index)
      endif
!WAS

      return
      end
c-----------------------------------------------------------------------


!     subroutine deriv_nonlocj(iel,x,rshift,rr_en,rr_en2,dk_en,dk_en2,value,gn) ! WAS
      subroutine deriv_nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,dk_en,dk_en2,value,gn)
c Written by Claudia Filippi, modified by Cyrus Umrigar
      use constants_mod
      use control_mod
      use atom_mod
      use dets_mod
      use optim_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use contrl_per_mod
      use derivjas_mod
      use jaspar_mod
      use bparm_mod
      use pointer_mod
      use periodic_mod
      use jaso_mod
      implicit real*8(a-h,o-z)


      dimension x(3,*),rshift(3,nelec,ncent),rr_en(nelec,ncent),rr_en2(nelec,ncent)
     &,gn(*),fsn(nelec,nelec),dx(3)
      dimension dk_en(nelec,ncent),dk_en2(nelec,ncent)

!WAS
      dimension r_en(nelec,ncent)
!!

      fsumn=0
      do 5 iparm=1,nparmj+nparms
    5   gn(iparm)=gvalue(iparm)

      if(nelec.lt.2) goto 47

      ipara=nparma(1)
      if(ijas.ge.4.and.ijas.le.6) then
        do 7 it=2,nctype
    7     ipara=ipara+nparma(it)
      endif

c     do 15 i=1,nelec
c       if(i.ne.iel) then
c         rvec_ee(1)=x(1,i)-x(1,iel)
c         rvec_ee(2)=x(2,i)-x(2,iel)
c         rvec_ee(3)=x(3,i)-x(3,iel)
c         rij=rvec_ee(1)**2+rvec_ee(2)**2+rvec_ee(3)**2
c         rij=dsqrt(rij)
c         call scale_dist(rij,u(i),1)
c       endif
c  15 continue

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
        if(i.le.nup .or. j.gt.nup) then
          if(nspin2.ge.2) then
            is=2
            isb=is
            if(nspin2.eq.3 .and. j.gt.nup) then
             is=3
             isb=is
            endif
           else
            is=1
            isb=is
            if(nspin2b.eq.2) then
              isb=2
             elseif(nocuspb.eq.0) then
              sspinn=half
            endif
          endif
          ipar=1
         else
          is=1
          isb=1
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
        if(nparms.eq.1) call deriv_scale(rij,dk,dk2,dr,dr2,2,0)

        gns=0
        iparm0=ipara+nparms
        if(isb.eq.2) iparm0=iparm0+nparmb(1)
        fsn(i,j)=deriv_psibnl(u,dk,gn(iparm0+1),gns,isb,ipar)
c       fsn(i,j)=fsn(i,j) + deriv_psibnl(u(jj),gn(iparm0+1),isb,ipar)

        do 25 jparm=1,nparmb(isb)
          iparm=iparm0+jparm
   25     gn(iparm)=gn(iparm)-go(i,j,iparm)
        if(nparms.eq.1) gn(1)=gn(1)+gns-go(i,j,1)

c e-e-n terms
c The scaling is switched in deriv_psinl, so do not do it here.
c     if(isc.ge.12) call scale_dist(rij,u,3)
        call scale_dist(rij,u,4)
        if(nparms.eq.1) call deriv_scale(rij,dk,dk2,dr,dr2,4,0)
        gns=0
        do 40 ic=1,ncent
          it=iwctype(ic)
c         ri=r_en(i,ic)
c         rj=r_en(j,ic)
          iparm0=npoint(it)+nparms
   40     fsn(i,j)=fsn(i,j) +
c     &    deriv_psinl(u,dk,rshift(1,i,ic),rshift(1,j,ic),rr_en2(i,ic),rr_en2(j,ic)
c     &               ,dk_en2(i,ic),dk_en2(j,ic),gn(iparm0+1),gns,it)
c    &    deriv_psinl(u(jj),rr_en2(i,ic),rr_en2(j,ic),gn(iparm0+1),it)
cc WAS
cc
     &    deriv_psinl(u,dk,rshift(1,i,ic),rshift(1,j,ic),r_en(i,ic),r_en(j,ic),
     &         rr_en2(i,ic),rr_en2(j,ic)
     &               ,dk_en2(i,ic),dk_en2(j,ic),gn(iparm0+1),gns,it)
CCC

        do 42 it=1,nctype
          iparm0=npoint(it)+nparms
          do 42 jparm=1,nparmc(it)
            iparm=iparm0+jparm
   42       gn(iparm)=gn(iparm)-go(i,j,iparm)
        if(nparms.eq.1) gn(1)=gn(1)+gns

        fsumn=fsumn+fsn(i,j)-fso(i,j)
   45 continue

c e-n terms
   47 fsn(iel,iel)=0

      if(ijas.eq.3) then
        if(nparms.ne.0) stop 'ijas.eq.3 not available for nparms.eq.1'
        do 50 ic=1,ncent
   50     fsn(iel,iel)=fsn(iel,iel)+deriv_psianl(rr_en(iel,ic),dk_en(iel,ic),gn(1),gns,1)
        do 51 iparm=1,nparma(1)
   51     gn(iparm)=gn(iparm)-go(iel,iel,iparm)
       elseif(ijas.ge.4.and.ijas.le.6) then
        gns=0
        do 55 ic=1,ncent
          it=iwctype(ic)
          iparm0=npointa(it)+nparms
   55     fsn(iel,iel)=fsn(iel,iel)+
     &                 deriv_psianl(rr_en(iel,ic),dk_en(iel,ic),gn(iparm0+1),gns,it)
        do 56 it=1,nctype
          iparm0=npointa(it)+nparms
          do 56 jparm=1,nparma(it)
            iparm=iparm0+jparm
   56       gn(iparm)=gn(iparm)-go(iel,iel,iparm)
        if(nparms.eq.1) gn(1)=gn(1)+gns-go(iel,iel,1)
      endif

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)
      value=fsumn

      return
      end

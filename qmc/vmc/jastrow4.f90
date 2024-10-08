      subroutine jastrow4(x,v,d2,div_vj,value)
! Written by Cyrus Umrigar
! Jastrow 4,5 must be used with one of isc=2,4,6,7,8,10,12,14,16,17
! Jastrow 6   must be used with one of isc=6,7

! When iperiodic=1, nloc=-4 (infinite quantum wire), then the "en" distance
!  used in en and een terms for the first center (ic = 1) is the distance
!  between the electron and the middle (y-axis) of the wire, rather than
!   the distance between (0,0) and an electron.  ACM, July 2010

      use constants_mod
      use control_mod
      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use wfsec_mod
      use contrl_per_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      use bparm_mod
      use distance_mod
      use jaso_mod
      use pseudo_mod
      use contrldmc_mod, only: l_tau_diffusion, l_psi_approx
      use fragments_mod, only: l_fragments
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-12)

! added WAS
      common /jas_c_cut/ cutjasc,icutjasc

      common /focktmp/ fc,fcu,fcuu,fcs,fcss,fct,fctt,fcst,fcus,fcut

      dimension x(3,*),v(3,*),div_vj(*)
      dimension uu(-2:max(nord,nordb,nordc)),ss(-2:max(nord,norda,nordc)),tt(-2:max(nord,norda,nordc)),rri(-2:max(nord,norda,nordc)) &
     &,rrj(-2:max(nord,norda,nordc))

      ndim1=ndim-1

      fsum=0
      if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) lapjo(1:nelec)=0d0

      do 5 i=-2,-1
        uu(i)=0
        ss(i)=0
        tt(i)=0
        rri(i)=0
    5   rrj(i)=0
      uu(0)=1
      ss(0)=2
      tt(0)=1
      rri(0)=1
      rrj(0)=1

      if(nelec.lt.2) goto 65

! e-e and e-e-n terms
      ij=0
      do 60 i=2,nelec
      im1=i-1
      do 60 j=1,im1
      ij=ij+1

      fso(i,j)=0
      d2ijo(i,j)=0
      if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then
        lapjijo(i,j)=0d0
        lapjijo(j,i)=0d0
      endif
      do 10 k=1,ndim
        fijo(k,i,j)=0
   10   fijo(k,j,i)=0

      sspinn=1
      isb=1
      ipar=0
      if(i.le.nup .or. j.gt.nup) then
        if(nspin2b.eq.2) then
          isb=2
         elseif(nocuspb.eq.0) then
          if(ndim.eq.3) then
            sspinn=half
           elseif(ndim.eq.2) then
            sspinn=third
          endif
        endif
        ipar=1
      endif

      rij=r_ee(ij)

      call scale_dist2(rij,uu(1),dd1,dd2,2)
!     write(6,'(''rij,u in ee'',2f9.5)') rij,uu(1)

! Check rij after scaling because uu(1) used in e-e-n terms too
      if(rij.gt.cutjas_ee) goto 30

      top=sspinn*b(1,isb,iwf)*uu(1)
      topu=sspinn*b(1,isb,iwf)
      topuu=0

      if(ijas.eq.4.or.ijas.eq.5) then
        bot=1+b(2,isb,iwf)*uu(1)
        botu=b(2,isb,iwf)
       elseif(ijas.eq.6) then
        bot=1+b(2,isb,iwf)*(1-uu(1))
        botu=-b(2,isb,iwf)
      endif
      botuu=0
      bot2=bot*bot

!      feeu=topu/bot-botu*top/bot2
!      feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
!      feeuu=feeuu/bot
! simpler expressions are :
      fee=top/bot
      feeu=topu/bot2
      feeuu=-2*feeu*botu/bot

      if(isc.eq.8 .or. isc.eq.10) then
        fee=fee/scalek(iwf)
        feeu=feeu/scalek(iwf)
        feeuu=feeuu/scalek(iwf)
      endif
      fee=fee-asymp_jasb(ipar+1,iwf)

      do 20 iord=2,nordb
        uu(iord)=uu(1)*uu(iord-1)
        if(ijas.eq.4) then
          fee=fee+b(iord+1,isb,iwf)*uu(iord)
          feeu=feeu+b(iord+1,isb,iwf)*iord*uu(iord-1)
          feeuu=feeuu+b(iord+1,isb,iwf)*iord*(iord-1)*uu(iord-2)
         elseif(ijas.eq.5.or.ijas.eq.6) then
          fee=fee+sspinn*b(iord+1,isb,iwf)*uu(iord)
          feeu=feeu+sspinn*b(iord+1,isb,iwf)*iord*uu(iord-1)
          feeuu=feeuu+sspinn*b(iord+1,isb,iwf)*iord*(iord-1)*uu(iord-2)
        endif
   20 continue

      feeuu=feeuu*dd1*dd1+feeu*dd2
      feeu=feeu*dd1/rij

      fso(i,j)=fso(i,j)+fee
      if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then
        lapjijo(i,j)=lapjijo(i,j)+feeuu+ndim1*feeu
        lapjijo(j,i)=lapjijo(j,i)+feeuu+ndim1*feeu
      endif
      do 21 k=1,ndim
        fijo(k,i,j)= fijo(k,i,j) + feeu*rvec_ee(k,ij)
   21   fijo(k,j,i)= fijo(k,j,i) - feeu*rvec_ee(k,ij)
      d2ijo(i,j)=d2ijo(i,j)+2.*(feeuu+ndim1*feeu)

! There are no C terms to order 1.
   30 if(nordc.le.1) goto 58

!     if(isc.ge.12) call scale_dist2(rij,uu(1),dd1,dd2,3)
      call scale_dist2(rij,uu(1),dd1,dd2,4)
      if(ijas.eq.4.or.ijas.eq.5) then
        call switch_scale2(uu(1),dd1,dd2,4)
        do 35 iord=2,nordc
   35     uu(iord)=uu(1)*uu(iord-1)
      endif
!     write(6,'(''rij,u in een'',2f12.9)') rij,uu(1)

      do 57 ic=1,ncent
        it=iwctype(ic)
!       if we have an infinite wire, then for the en and een terms, for ic=1,
!         we only use the distance from the electron to the y-axis of the wire
        if((iperiodic.eq.1).and.(nloc.eq.-4).and.(ic.eq.1)) then
           ri=abs(rvec_en(2,i,ic))  ! y-component of rvec_en is dist to wire center
           rj=abs(rvec_en(2,j,ic))
        else
           ri=r_en(i,ic)
           rj=r_en(j,ic)
        endif

!       write(6,'(''ri,rj'',9f9.4)') ri,rj
!       write(6,'(''rshift(k,i,ic)'',9f9.4)') (rshift(k,i,ic),k=1,ndim),(rshift(k,j,ic),k=1,ndim),(rshift(k,i,ic)-rshift(k,j,ic),k=1,ndim)

        if(ri.gt.cutjas_en .or. rj.gt.cutjas_en) goto 57
        do 37 k=1,ndim
   37     if(abs(rshift(k,i,ic)-rshift(k,j,ic)).gt.eps) goto 57

        call scale_dist2(ri,rri(1),dd7,dd9,3)
        call scale_dist2(rj,rrj(1),dd8,dd10,3)

        if(ijas.eq.4.or.ijas.eq.5) then
          call switch_scale2(rri(1),dd7,dd9,3)
          call switch_scale2(rrj(1),dd8,dd10,3)
        endif
!       write(6,'(''ri,rri in een'',2f12.9)') ri,rri(1)

! Moved to line 252
!       if(icutjasc .gt. 0 .or. iperiodic .ne. 0) then
!          call f_een_cuts (cutjas_en, ri, rj, fcuti, fcutj, fcut, dfcuti, dfcutj,d2fcuti,d2fcutj)
!       endif

        s=ri+rj
        t=ri-rj
!       u2mt2=rij*rij-t*t
        u2pst=rij*rij+s*t
        u2mst=rij*rij-s*t
!       s2mu2=s*s-rij*rij
!       s2mt2=s*s-t*t

        do 50 iord=1,nordc
          rri(iord)=rri(1)*rri(iord-1)
          rrj(iord)=rrj(1)*rrj(iord-1)
          ss(iord)=rri(iord)+rrj(iord)
   50     tt(iord)=rri(iord)*rrj(iord)

        fc=0
        fu=0
        fuu=0
        fi=0
        fii=0
        fj=0
        fjj=0
        fui=0
        fuj=0
        ll=0
        do 55 n=2,nordc
          do 55 k=n-1,0,-1
            if(k.eq.0) then
              l_hi=n-k-2
             else
              l_hi=n-k
            endif
            do 55 l=l_hi,0,-1
              m=(n-k-l)/2
              if(2*m.eq.n-k-l) then
                ll=ll+1
                fc=fc+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
                fu=fu+c(ll,it,iwf)*k*uu(k-1)*ss(l)*tt(m)
                fuu=fuu+c(ll,it,iwf)*k*(k-1)*uu(k-2)*ss(l)*tt(m)
                fi=fi+c(ll,it,iwf)*uu(k) &
     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                fii=fii+c(ll,it,iwf)*uu(k) &
     &          *((l+m)*(l+m-1)*rri(l+m-2)*rrj(m) &
     &          +m*(m-1)*rri(m-2)*rrj(l+m))
                fj=fj+c(ll,it,iwf)*uu(k) &
     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
                fjj=fjj+c(ll,it,iwf)*uu(k) &
     &          *((l+m)*(l+m-1)*rrj(l+m-2)*rri(m) &
     &          +m*(m-1)*rrj(m-2)*rri(l+m))
                fui=fui+c(ll,it,iwf)*k*uu(k-1) &
     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                fuj=fuj+c(ll,it,iwf)*k*uu(k-1) &
     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
              endif
!             write(6,'(''rij,ri,rj'',9f10.5)') rij,ri,rj,uu(1),rri(1),rrj(1)
   55   continue

        if(ifock.gt.0) call fock(uu(1),ss(1),tt(1),rri(1),rrj(1),it)

        fuu=fuu*dd1*dd1+fu*dd2
        fu=fu*dd1/rij

        fui=fui*dd1*dd7
        fuj=fuj*dd1*dd8

        fii=fii*dd7*dd7+fi*dd9
        fjj=fjj*dd8*dd8+fj*dd10
        fi=fi*dd7/ri
        fj=fj*dd8/rj
! een for periodic systems         WAS
        if(icutjasc .gt. 0 .or. iperiodic .ne. 0) then
           call f_een_cuts (cutjas_en, ri, rj, fcuti, fcutj, fcut, dfcuti, dfcutj,d2fcuti,d2fcutj)

           fuu = fuu * fcut
           fii = fii * fcut +(2 * fi * ri * dfcuti + fc * d2fcuti)*fcutj
           fi = fi * fcut + (fc * fcutj *  dfcuti)/ri
           fui = fui * fcut + (fu * fcutj *  dfcuti*rij)

           fjj = fjj * fcut + (2 * fj * dfcutj *rj + fc * d2fcutj)*fcuti
           fj = fj * fcut + (fc * fcuti *  dfcutj)/rj
           fuj = fuj * fcut + (fu * fcuti *  dfcutj * rij)
           fc = fc * fcut
           fu = fu * fcut
        endif
! end WAS

        fso(i,j)=fso(i,j) + fc

!       if we have an infinite wire, then for the en and een terms, we only
!          use the distance from the electron to the center of the wire
!          (for the first center, anyway)
!          thus, derivatives only depend on en distance in the y-direction
        if((iperiodic.eq.1).and.(nloc.eq.-4).and.(ic.eq.1)) then
           fijo(1,i,j)=fijo(1,i,j)                     +fu*rvec_ee(1,ij)
           fijo(2,i,j)=fijo(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
           fijo(3,i,j)=fijo(3,i,j)                     +fu*rvec_ee(3,ij)
           fijo(1,j,i)=fijo(1,j,i)                     -fu*rvec_ee(1,ij)
           fijo(2,j,i)=fijo(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
           fijo(3,j,i)=fijo(3,j,i)                     -fu*rvec_ee(3,ij)
        else
           fijo(1,i,j)=fijo(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
           fijo(2,i,j)=fijo(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
           fijo(3,i,j)=fijo(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
           fijo(1,j,i)=fijo(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
           fijo(2,j,i)=fijo(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
           fijo(3,j,i)=fijo(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)
        endif
!       write(6,'(''i,j,fijo2='',2i5,9d12.4)') i,j,(fijo(k,i,j),k=1,ndim)

!       d2ijo(i,j)=d2ijo(i,j) + 2*(fuu + 2*fu) + fui*u2pst/(ri*rij)
!    &  + fuj*u2mst/(rj*rij) + fii + 2*fi + fjj + 2*fj

!      This equation (and the subsequent u,s,t notation) comes from
!        eqn 5 of Pekeris, Phys Rev., 112, 1649 (1958), which
!        is derived in eqns. 3-5 of Hylleraas, Z. Physik 54, 347 (1929)
!      The expression for wires, where /psi = /psi(y_1, y_2, r_12) rather than
!        /psi(r_1, r_2, r_12) follows from the same type of calculation (ACM)
        if((iperiodic.eq.1).and.(nloc.eq.-4).and.(ic.eq.1)) then ! ri = yi
          d2ijo(i,j)=d2ijo(i,j) + ndim1*2.*fu &
     &         + 2.*fuu + fii + fjj + 2.*fui*(ri-rj)/rij + 2.*fuj*(rj-ri)/rij
          if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then
            lapjijo(i,j)=lapjijo(i,j) + ndim1*fu + fuu + fii + 2.*fui*(ri-rj)/rij
            lapjijo(j,i)=lapjijo(j,i) + ndim1*fu + fuu + fjj + 2.*fuj*(rj-ri)/rij
          endif
        else  ! Pekeris, Phys Rev 112, 1649 (1958) eqn 5:
          d2ijo(i,j)=d2ijo(i,j) + ndim1*(2*fu+fi+fj) &
     &          + 2*fuu + fii +  fjj + fui*u2pst/(ri*rij) + fuj*u2mst/(rj*rij)
          if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then
            lapjijo(i,j)=lapjijo(i,j) + ndim1*(fu+fi) + fuu + fii + fui*u2pst/(ri*rij)
            lapjijo(j,i)=lapjijo(j,i) + ndim1*(fu+fj) + fuu + fjj + fuj*u2mst/(rj*rij)
          endif
        endif


   57 continue

   58 fsum=fsum+fso(i,j)
      v(1,i)=v(1,i)+fijo(1,i,j)
      v(2,i)=v(2,i)+fijo(2,i,j)
      v(3,i)=v(3,i)+fijo(3,i,j)
      v(1,j)=v(1,j)+fijo(1,j,i)
      v(2,j)=v(2,j)+fijo(2,j,i)
      v(3,j)=v(3,j)+fijo(3,j,i)
      div_vj(i)=div_vj(i)+d2ijo(i,j)/2
      div_vj(j)=div_vj(j)+d2ijo(i,j)/2
      if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then
        lapjo(i)=lapjo(i)+lapjijo(i,j)
        lapjo(j)=lapjo(j)+lapjijo(j,i)
      endif
   60 d2=d2+d2ijo(i,j)

! e-n terms
   65 do 90 i=1,nelec

        fso(i,i)=0
        fijo(1,i,i)=0
        fijo(2,i,i)=0
        fijo(3,i,i)=0
        d2ijo(i,i)=0
        if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then
          lapjijo(i,i)=0d0
        endif

        do 80 ic=1,ncent
          it=iwctype(ic)
!       if we have an infinite wire, then for the en and een terms, for ic=1,
!         we only use the distance from the electron to the y-axis of the wire
          if((iperiodic.eq.1).and.(nloc.eq.-4).and.(ic.eq.1)) then
             ri=abs(rvec_en(2,i,ic)) ! y-component of rvec_en is dist to wire center
          else
             ri=r_en(i,ic)
          endif

          if(ri.gt.cutjas_en) goto 80

          call scale_dist2(ri,rri(1),dd7,dd9,1)
!         write(6,'(''ri,rri in en'',2f9.5)') ri,rri(1)

          top=a4(1,it,iwf)*rri(1)
          topi=a4(1,it,iwf)
          topii=0

          if(ijas.eq.4.or.ijas.eq.5) then
            bot=1+a4(2,it,iwf)*rri(1)
            boti=a4(2,it,iwf)
           elseif(ijas.eq.6) then
            bot=1+a4(2,it,iwf)*(1-rri(1))
            boti=-a4(2,it,iwf)
          endif
          botii=0
          bot2=bot*bot

!          feni=topi/bot-boti*top/bot2
!          fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
!          fenii=fenii/bot
! simpler expressions are :
          fen=top/bot
          feni=topi/bot2
          fenii=-2*feni*boti/bot

          if(isc.eq.8 .or. isc.eq.10) then
            fen=fen/scalek(iwf)
            feni=feni/scalek(iwf)
            fenii=fenii/scalek(iwf)
          endif
          fen=fen-asymp_jasa(it,iwf)

          do 70 iord=2,norda
            rri(iord)=rri(1)**iord
            fen=fen+a4(iord+1,it,iwf)*rri(iord)
            feni=feni+a4(iord+1,it,iwf)*iord*rri(iord-1)
   70       fenii=fenii+a4(iord+1,it,iwf)*iord*(iord-1)*rri(iord-2)

          fenii=fenii*dd7*dd7+feni*dd9
          feni=feni*dd7/ri

          fso(i,i)=fso(i,i)+fen

!       if we have an infinite wire, then for the en and een terms, we only
!          use the distance from the electron to the center of the wire
!          (for the first center, anyway)
!          thus, derivatives only depend on en distance in the y-direction
          if((iperiodic.eq.1).and.(nloc.eq.-4).and.(ic.eq.1)) then
             fijo(2,i,i)=fijo(2,i,i) + feni*rvec_en(2,i,ic)
             d2ijo(i,i) = d2ijo(i,i) + fenii  ! check this...
             if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then 
               lapjijo(i,i) = lapjijo(i,i) + fenii
             endif
          else
             fijo(1,i,i)=fijo(1,i,i) + feni*rvec_en(1,i,ic)
             fijo(2,i,i)=fijo(2,i,i) + feni*rvec_en(2,i,ic)
             fijo(3,i,i)=fijo(3,i,i) + feni*rvec_en(3,i,ic)
             d2ijo(i,i) = d2ijo(i,i) + fenii + ndim1*feni
             if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then
               lapjijo(i,i) = lapjijo(i,i) + fenii + ndim1*feni
             endif
          endif

!         write(6,'(''fijo='',9d12.4)') (fijo(k,i,i),k=1,ndim),feni,rvec_en(1,i,ic)

!   Place this line in the if..else block above (ACM):
!          d2ijo(i,i) = d2ijo(i,i) + fenii + ndim1*feni
   80   continue

        fsum=fsum+fso(i,i)
        v(1,i)=v(1,i)+fijo(1,i,i)
        v(2,i)=v(2,i)+fijo(2,i,i)
        v(3,i)=v(3,i)+fijo(3,i,i)
!       write(6,'(''v='',9d12.4)') (v(k,i),k=1,ndim)
        div_vj(i)=div_vj(i)+d2ijo(i,i)
        if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then
          lapjo(i)=lapjo(i)+lapjijo(i,i)
        endif
   90   d2=d2+d2ijo(i,i)

! Warning: c1_jas6 below needs changing now that we have different
! c1_jas6_en and c1_jas6_ee.  Since I no longer use ijaas=6, I do not bother.
      if(ijas.eq.6) then
        stop 'ijas6 needs updating'
        term=1/(c1_jas6*scalek(iwf))
        fsum=term*fsum
        d2=term*d2
        do 100 i=1,nelec
          div_vj(i)=term*div_vj(i)
          do 95 k=1,ndim
   95       v(k,i)=term*v(k,i)
          do 100 j=1,nelec
            d2ijo(i,j)=term*d2ijo(i,j)
            do 100 k=1,ndim
  100         fijo(k,i,j)=term*fijo(k,i,j)
      endif

      if (l_fragments.OR.l_tau_diffusion.OR.l_psi_approx) then
        do i=1,nelec
          do k=1,ndim
            lapjo(i)=lapjo(i)+v(k,i)**2
          enddo
        enddo
      endif

      fsumo=fsum
      d2o=d2
      do 110 i=1,nelec
        fjo(1,i)=v(1,i)
        fjo(2,i)=v(2,i)
  110   fjo(3,i)=v(3,i)

      value=fsum

      return
      end

!-----------------------------------------------------------------------
      integer function nterms4(nord)
! Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)

      i=0
      do 20 n=2,nord
        do 20 k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do 20 l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              i=i+1
            endif
   20 continue
      nterms4=i
!     write(6,'(''nterms4='',i5)') nterms4
      return
      end

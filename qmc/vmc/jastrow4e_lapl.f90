      subroutine jastrow4e_lapl(iel,x,v,value)
! Written by Tyler Anderson, by modifying jastrow4.f90 and jastrow4e.f90
! Jastrow 4,5 must be used with one of isc=2,4,6,7,8,10,12,14,16,17
! Jastrow 6   must be used with one of isc=6,7

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
      use jasn_mod
      use pseudo_mod
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-12)

! added WAS
      common /jas_c_cut/ cutjasc,icutjasc

      common /focktmp/ fc,fcu,fcuu,fcs,fcss,fct,fctt,fcst,fcus,fcut

      dimension x(3,*),v(3,*)
      dimension uu(-2:max(nord,nordb,nordc)),ss(-2:max(nord,norda,nordc)),tt(-2:max(nord,norda,nordc)),rri(-2:max(nord,norda,nordc)) &
     &,rrj(-2:max(nord,norda,nordc))

      ndim1=ndim-1

      fsum=0

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

      do i=1,nelec
        do k=1,ndim
          fjn(k,i)=fjo(k,i)
        enddo
        lapjn(i)=lapjo(i)
      enddo
      fsumn=fsumo

      if(nelec.lt.2) goto 65

! e-e and e-e-n terms
      do jj=1,nelec

      if(jj.eq.iel) cycle
      if(jj.lt.iel) then
        i=iel
        j=jj
      else
        i=jj
        j=iel
      endif
      ij=((i-1)*(i-2))/2+j

      do k=1,ndim
        fijn(k,i,j)=0d0
        fijn(k,j,i)=0d0
      enddo
      lapjijn(i,j)=0d0
      lapjijn(j,i)=0d0
      fsn(i,j)=0d0

      sspinn=1
      ipar=0
      isb=1
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
!     topuu=0

      if(ijas.eq.4.or.ijas.eq.5) then
        bot=1+b(2,isb,iwf)*uu(1)
        botu=b(2,isb,iwf)
       elseif(ijas.eq.6) then
        bot=1+b(2,isb,iwf)*(1-uu(1))
        botu=-b(2,isb,iwf)
      endif
!      botuu=0
      bot2=bot*bot

!      feeu=topu/bot-botu*top/bot2
!      feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
!      feeuu=feeuu/bot
      fee=top/bot
      feeu=topu/bot2
      feeuu=-2*feeu*botu/bot
      if(isc.eq.8 .or. isc.eq.10) then
        fee=fee/scalek(iwf)
        feeu=feeu/scalek(iwf)
        feeuu=feeuu/scalek(iwf)
      endif
      fee=fee-asymp_jasb(ipar+1,iwf)

      do iord=2,nordb
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
      enddo !iord=2,nordb

      feeuu=feeuu*dd1*dd1+feeu*dd2
      feeu=feeu*dd1/rij

      fsn(i,j)=fsn(i,j)+fee
      lapjijn(i,j)=lapjijn(i,j)+feeuu+ndim1*feeu
      lapjijn(j,i)=lapjijn(j,i)+feeuu+ndim1*feeu
      do k=1,ndim
        fijn(k,i,j)= fijn(k,i,j) + feeu*rvec_ee(k,ij)
        fijn(k,j,i)= fijn(k,j,i) - feeu*rvec_ee(k,ij)
      enddo

! There are no C terms to order 1.
   30 if(nordc.le.1) goto 58

!     if(isc.ge.12) call scale_dist2(rij,uu(1),dd1,dd2,3)
      call scale_dist2(rij,uu(1),dd1,dd2,4)
      if(ijas.eq.4.or.ijas.eq.5) then
        call switch_scale2(uu(1),dd1,dd2,4)
        do iord=2,nordc
          uu(iord)=uu(1)*uu(iord-1)
        enddo
      endif

!       Can't call this yet! ri, rj not defined! (ACM)
!       if(icutjasc .gt. 0 .or. iperiodic .ne. 0) then
!         call f_een_cuts (cutjas_en, ri, rj, fcuti, fcutj, fcut, dfcuti, dfcutj    ,d2fcuti,d2fcutj)
!       endif      

      do ic=1,ncent
        it=iwctype(ic)
                                                                             
        if((iperiodic.eq.1).and.(nloc.eq.-4).and.(ic.eq.1)) then
           ri=abs(rvec_en(2,i,ic))  
           rj=abs(rvec_en(2,j,ic))
        else
           ri=r_en(i,ic)
           rj=r_en(j,ic)
        endif

        if(ri.gt.cutjas_en .or. rj.gt.cutjas_en) cycle
        do k=1,ndim
          if(abs(rshift(k,i,ic)-rshift(k,j,ic)).gt.eps) cycle
        enddo

        call scale_dist2(ri,rri(1),dd7,dd9,3)
        call scale_dist2(rj,rrj(1),dd8,dd10,3)

        if(ijas.eq.4.or.ijas.eq.5) then
          call switch_scale2(rri(1),dd7,dd9,3)
          call switch_scale2(rrj(1),dd8,dd10,3)
        endif

! Moved back here, from above the do 50 loop (ACM)
!   Don't think we should call it here either, since we call it just after the c    ontinue at 40
!         if(icutjasc .gt. 0 .or. iperiodic .ne. 0) then
!            call f_een_cuts (cutjas_en, ri, rj, fcuti, fcutj, fcut, dfcuti, dfc    utj, d2fcuti, d2fcutj)
!         endif

        s=ri+rj
        t=ri-rj
        u2pst=rij*rij+s*t
        u2mst=rij*rij-s*t

        do iord=1,nordc
          rri(iord)=rri(1)*rri(iord-1)
          rrj(iord)=rrj(1)*rrj(iord-1)
          ss(iord)=rri(iord)+rrj(iord)
          tt(iord)=rri(iord)*rrj(iord)
        enddo

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
        do n=2,nordc
          do k=n-1,0,-1
            if(k.eq.0) then
              l_hi=n-k-2
             else
              l_hi=n-k
            endif
            do l=l_hi,0,-1
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
            enddo !l=l_hi,0,-1
          enddo !k=n-1,0,-1
        enddo !n=2,nordc

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
        fsn(i,j)=fsn(i,j) + fc

        if((iperiodic.eq.1).and.(nloc.eq.-4).and.(ic.eq.1)) then  ! infinite wires
          fijn(1,i,j)=fijn(1,i,j) + fu*rvec_ee(1,ij)
          fijn(2,i,j)=fijn(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
          fijn(3,i,j)=fijn(3,i,j) + fu*rvec_ee(3,ij)
          fijn(1,j,i)=fijn(1,j,i) - fu*rvec_ee(1,ij)
          fijn(2,j,i)=fijn(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
          fijn(3,j,i)=fijn(3,j,i) - fu*rvec_ee(3,ij)
!         d2ijn(i,j)=d2ijn(i,j) + ndim1*2.*fu &
!     &       + 2.*fuu + fii + fjj + 2.*fui*(ri-rj)/rij + 2.*fuj*(rj-ri)/rij
          lapjijn(i,j)=lapjijn(i,j) + ndim1*fu + fuu + fii + 2.*fui*(ri-rj)/rij
          lapjijn(j,i)=lapjijn(j,i) + ndim1*fu + fuu + fjj + 2.*fuj*(rj-ri)/rij
        else
          fijn(1,i,j)=fijn(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
          fijn(2,i,j)=fijn(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
          fijn(3,i,j)=fijn(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
          fijn(1,j,i)=fijn(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
          fijn(2,j,i)=fijn(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
          fijn(3,j,i)=fijn(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)
!         Pekeris, Phys Rev 112, 1649 (1958) eqn 5:
!         d2ijn(i,j)=d2ijn(i,j) + ndim1*(2*fu+fi+fj) &
!     &        + 2*fuu + fii +  fjj + fui*u2pst/(ri*rij) + fuj*u2mst/(rj*rij)
          lapjijn(i,j)=lapjijn(i,j) + ndim1*(fu+fi) + fuu + fii + fui*u2pst/(ri*rij)
          lapjijn(j,i)=lapjijn(j,i) + ndim1*(fu+fj) + fuu + fjj + fuj*u2mst/(rj*rij)
        endif

      enddo !ic=1,ncent

   58 fsumn=fsumn+fsn(i,j)-fso(i,j)
      fjn(1,i)=fjn(1,i)+fijn(1,i,j)-fijo(1,i,j)
      fjn(2,i)=fjn(2,i)+fijn(2,i,j)-fijo(2,i,j)
      fjn(3,i)=fjn(3,i)+fijn(3,i,j)-fijo(3,i,j)
      fjn(1,j)=fjn(1,j)+fijn(1,j,i)-fijo(1,j,i)
      fjn(2,j)=fjn(2,j)+fijn(2,j,i)-fijo(2,j,i)
      fjn(3,j)=fjn(3,j)+fijn(3,j,i)-fijo(3,j,i)
      lapjn(i)=lapjn(i)+lapjijn(i,j)-lapjijo(i,j)
      lapjn(j)=lapjn(j)+lapjijn(j,i)-lapjijo(j,i)

      enddo ! jj=1,nelec

! e-n terms
   65 fijn(1,iel,iel)=0d0
      fijn(2,iel,iel)=0d0
      fijn(3,iel,iel)=0d0
      fsn(iel,iel)=0d0
      lapjijn(iel,iel)=0d0

      do 80 ic=1,ncent
        it=iwctype(ic)

        if((iperiodic.eq.1).and.(nloc.eq.-4).and.(ic.eq.1)) then
           ri=abs(rvec_en(2,iel,ic))
        else
           ri=r_en(iel,ic)
        endif
        if(ri.gt.cutjas_en) goto 80

        call scale_dist2(ri,rri(1),dd7,dd9,1)

        top=a4(1,it,iwf)*rri(1)
        topi=a4(1,it,iwf)
!        topii=0

        if(ijas.eq.4.or.ijas.eq.5) then
          bot=a4(2,it,iwf)*rri(1)
          boti=a4(2,it,iwf)
         elseif(ijas.eq.6) then
          bot=1+a4(2,it,iwf)*(1-rri(1))
          boti=-a4(2,it,iwf)
        endif
!        botii=0

        bot=1+bot
        bot2=bot*bot
!       feni=topi/bot-boti*top/bot2
!       fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
!       fenii=fenii/bot
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

        do iord=2,norda
          rri(iord)=rri(1)**iord
          fen=fen+a4(iord+1,it,iwf)*rri(iord)
          feni=feni+a4(iord+1,it,iwf)*iord*rri(iord-1)
          fenii=fenii+a4(iord+1,it,iwf)*iord*(iord-1)*rri(iord-2)
        enddo
        
        fenii=fenii*dd7*dd7+feni*dd9
        feni=feni*dd7/ri

        fsn(iel,iel)=fsn(iel,iel)+fen

        if((iperiodic.eq.1).and.(nloc.eq.-4).and.(ic.eq.1)) then !infinte writes
          fijn(2,iel,iel)=fijn(2,iel,iel) + feni*rvec_en(2,iel,ic)
          lapjijn(iel,iel) = lapjijn(iel,iel) + fenii  ! check this...
        else
          fijn(1,iel,iel)=fijn(1,iel,iel) + feni*rvec_en(1,iel,ic)
          fijn(2,iel,iel)=fijn(2,iel,iel) + feni*rvec_en(2,iel,ic)
          fijn(3,iel,iel)=fijn(3,iel,iel) + feni*rvec_en(3,iel,ic)
          lapjijn(iel,iel) = lapjijn(iel,iel) + fenii + ndim1*feni
        endif

   80   continue

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)
      fjn(1,iel)=fjn(1,iel)+fijn(1,iel,iel)-fijo(1,iel,iel)
      fjn(2,iel)=fjn(2,iel)+fijn(2,iel,iel)-fijo(2,iel,iel)
      fjn(3,iel)=fjn(3,iel)+fijn(3,iel,iel)-fijo(3,iel,iel)
      lapjn(iel)=lapjn(iel)+lapjijn(iel,iel)-lapjijo(iel,iel)

      do i=1,nelec
        do k=1,ndim
          lapjn(i)=lapjn(i) + fjn(k,i)**2 - fjo(k,i)**2
        enddo
      enddo
                       
      do i=1,nelec
        v(1,i)=fjn(1,i)
        v(2,i)=fjn(2,i)
        v(3,i)=fjn(3,i)
      enddo
                     
      value=fsumn     

      return
      end

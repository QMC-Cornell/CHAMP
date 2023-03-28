      subroutine plot_jas
! Written by Cyrus Umrigar
! subroutine to plot e-n and e-e Jastrow

      use atom_mod
      use contr2_mod
      use jaspar6_mod
      implicit real*8(a-h,o-z)


      dimension fen(ncent),fee(2)

      write(6,'('' r    f_en(1-nctype)    f_ee(1,2)'')') ; call systemflush(6)
      np=101
! Warning: I should really be plotting f_en to cutjas_en and f_ee to cutjas_ee but
! for the moment plot both to larger distance.
      if(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17) then
        rmax=cutjas_ee
       else
        rmax=10.d0
      endif
      dr=rmax/(np-1)
      r=-dr
      do 10 i=1,np
        r=r+dr
        call jastrow4ab(r,fen,fee)
   10   write(6,'(f8.4,9f12.6)') r,(fen(it),it=1,nctype),(fee(isb),isb=1,2) ; call systemflush(6)

      return
      end
!-----------------------------------------------------------------------
      subroutine jastrow4ab(r,fen,fee)
! Written by Cyrus Umrigar
! Jastrow 4,5 must be used with one of isc=2,4,6,7
! Jastrow 6   must be used with one of isc=6,7

      use constants_mod
      use atom_mod
      use const_mod
      use contr2_mod
      use wfsec_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      use bparm_mod
      implicit real*8(a-h,o-z)

      dimension uu(-2:max(nord,nordb,nordc)),ss(-2:max(nord,norda,nordc)),tt(-2:max(nord,norda,nordc)),rri(-2:max(nord,norda,nordc)) &
     &,rrj(-2:max(nord,norda,nordc)),fen(ncent),fee(2)

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

! ipar=0 antiparallel
! ipar=1 parallel
      do 20 ipar=0,1

      if(ipar.eq.0) then
        sspinn=1
        isb=1
       else
        if(nspin2b.eq.2) then
          isb=2
         elseif(nocuspb.eq.0) then
          sspinn=half
        endif
      endif

      call scale_dist(r,uu(1),2)

      top=sspinn*b(1,isb,iwf)*uu(1)

      if(ijas.eq.4.or.ijas.eq.5) then
        bot=1+b(2,isb,iwf)*uu(1)
       elseif(ijas.eq.6) then
        bot=1+b(2,isb,iwf)*(1-uu(1))
      endif

      fee(ipar+1)=top/bot-asymp_jasb(ipar+1,iwf)

      do 20 iord=2,nordb
        uu(iord)=uu(1)*uu(iord-1)
        if(ijas.eq.4) then
          fee(ipar+1)=fee(ipar+1)+b(iord+1,isb,iwf)*uu(iord)
         elseif(ijas.eq.5.or.ijas.eq.6) then
          fee(ipar+1)=fee(ipar+1)+sspinn*b(iord+1,isb,iwf)*uu(iord)
        endif
   20 continue

! e-n terms
      do 70 it=1,nctype

          call scale_dist(r,rri(1),1)

          top=a4(1,it,iwf)*rri(1)

          if(ijas.eq.4.or.ijas.eq.5) then
            bot=1+a4(2,it,iwf)*rri(1)
           elseif(ijas.eq.6) then
            bot=1+a4(2,it,iwf)*(1-rri(1))
          endif
          bot2=bot*bot

          fen(it)=top/bot-asymp_jasa(it,iwf)

        do 70 iord=2,norda
          rri(iord)=rri(1)**iord
   70     fen(it)=fen(it)+a4(iord+1,it,iwf)*rri(iord)

      return
      end

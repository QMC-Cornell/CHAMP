      subroutine jastrow4e(iel,x,v,value)
c Written by Cyrus Umrigar and Claudia Filippi
c Jastrow 4,5 must be used with one of isc=2,4,6,7,12,14,16,17
c Jastrow 6   must be used with one of isc=6,7

      use control_mod
      use atom_mod
      use dets_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'force.h'

!JT      parameter (half=.5d0,third=1.d0/3.d0,eps=1.d-12)
      parameter (eps=1.d-12)

      common /dim/ ndim
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
!!! added WAS
      common /jas_c_cut/ cutjasc,icutjasc
      common /contrl_per/ iperiodic,ibasis
!!!
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /bparm/ nspin2b,nocuspb

      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
      common /focktmp/ fc,fcu,fcuu,fcs,fcss,fct,fctt,fcst,fcus,fcut

      common /jasn/ fsn(MELEC,MELEC),fijn(3,MELEC,MELEC)
     &,d2ijn(MELEC,MELEC),d2n,fsumn,fjn(3,MELEC)
      common /jaso/ fso(MELEC,MELEC),fijo(3,MELEC,MELEC)
     &,d2ijo(MELEC,MELEC),d2o,fsumo,fjo(3,MELEC)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      dimension x(3,*),v(3,*)
      dimension uu(-2:MORDJ),ss(-2:MORDJ),tt(-2:MORDJ),rri(-2:MORDJ)
     &,rrj(-2:MORDJ)

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

      do 10 i=1,nelec
        do 10 k=1,ndim
   10   fjn(k,i)=fjo(k,i)
      fsumn=fsumo

      if(nelec.lt.2) goto 65

      do 60 jj=1,nelec

      if(jj.eq.iel) goto 60
      if(jj.lt.iel) then
        i=iel
        j=jj
       else
        i=jj
        j=iel
      endif
      ij=((i-1)*(i-2))/2+j

      fijn(1,i,j)=0
      fijn(2,i,j)=0
      fijn(3,i,j)=0
      fijn(1,j,i)=0
      fijn(2,j,i)=0
      fijn(3,j,i)=0
      fsn(i,j)=0

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

      call scale_dist1(rij,uu(1),dd1,2)

      if(rij.gt.cutjas_ee) goto 22




      top=sspinn*b(1,isb,iwf)*uu(1)
      topu=sspinn*b(1,isb,iwf)
c     topuu=0

      bot=1+b(2,isb,iwf)*uu(1)
      botu=b(2,isb,iwf)
c     botuu=0
      bot2=bot*bot

c      fee=top/bot-asymp_jasb(ipar+1,iwf)
c      feeu=topu/bot-botu*top/bot2
      fee=top/bot
      feeu=topu/bot2
      if(isc.eq.8 .or. isc.eq.10) then
        fee=fee/scalek(iwf)
        feeu=feeu/scalek(iwf)
      endif
      fee=fee-asymp_jasb(ipar+1,iwf)

      do 20 iord=2,nordb
        uu(iord)=uu(1)*uu(iord-1)
        fee=fee+b(iord+1,isb,iwf)*uu(iord)
   20   feeu=feeu+b(iord+1,isb,iwf)*iord*uu(iord-1)

      feeu=feeu*dd1/rij

      fsn(i,j)=fsn(i,j) + fee

      fijn(1,i,j)=fijn(1,i,j) + feeu*rvec_ee(1,ij)
      fijn(2,i,j)=fijn(2,i,j) + feeu*rvec_ee(2,ij)
      fijn(3,i,j)=fijn(3,i,j) + feeu*rvec_ee(3,ij)
      fijn(1,j,i)=fijn(1,j,i) - feeu*rvec_ee(1,ij)
      fijn(2,j,i)=fijn(2,j,i) - feeu*rvec_ee(2,ij)
      fijn(3,j,i)=fijn(3,j,i) - feeu*rvec_ee(3,ij)

c There are no C terms to order 1.
   22 if(nordc.le.1) goto 55

c     if(isc.ge.12) call scale_dist1(rij,uu(1),dd1,3)
      call scale_dist1(rij,uu(1),dd1,4)
      if(ijas.eq.4.or.ijas.eq.5) then
        call switch_scale1(uu(1),dd1,4)
        do 25 iord=2,nordc
   25     uu(iord)=uu(1)*uu(iord-1)
      endif

      if(icutjasc .gt. 0 .or. iperiodic .ne. 0) then
c     call f_een_cuts (cutjasc, ri, rj, fcuti, fcutj, fcut,
         call f_een_cuts (cutjas_en, ri, rj, fcuti, fcutj, fcut,
     +        dfcuti, dfcutj,d2fcuti,d2fcutj)
      endif

      do 50 ic=1,ncent
        it=iwctype(ic)

        ri=r_en(i,ic)
        rj=r_en(j,ic)

        if(ri.gt.cutjas_en .or. rj.gt.cutjas_en) goto 50
        do 27 k=1,ndim
   27     if(abs(rshift(k,i,ic)-rshift(k,j,ic)).gt.eps) goto 50

        call scale_dist1(ri,rri(1),dd7,3)
        call scale_dist1(rj,rrj(1),dd8,3)

        if(ijas.eq.4.or.ijas.eq.5) then
          call switch_scale1(rri(1),dd7,3)
          call switch_scale1(rrj(1),dd8,3)
        endif

        if(icutjasc .gt. 0 .or. iperiodic .ne. 0) then
c     call f_een_cuts (cutjasc, ri, rj, fcuti, fcutj, fcut,
           call f_een_cuts (cutjas_en, ri, rj, fcuti, fcutj, fcut,
     +          dfcuti, dfcutj,d2fcuti,d2fcutj)
        endif


        do 30 iord=1,nordc
          rri(iord)=rri(1)*rri(iord-1)
          rrj(iord)=rrj(1)*rrj(iord-1)
          ss(iord)=rri(iord)+rrj(iord)
   30     tt(iord)=rri(iord)*rrj(iord)

    	fc=0
  	fu=0
  	fi=0
  	fj=0
	ll=0
        do 40 n=2,nordc
          do 40 k=n-1,0,-1
            if(k.eq.0) then
              l_hi=n-k-2
             else
              l_hi=n-k
            endif
            do 40 l=l_hi,0,-1
              m=(n-k-l)/2
              if(2*m.eq.n-k-l) then
                ll=ll+1
                fc=fc+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
                fu=fu+c(ll,it,iwf)*k*uu(k-1)*ss(l)*tt(m)
                fi=fi+c(ll,it,iwf)*uu(k)
     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                fj=fj+c(ll,it,iwf)*uu(k)
     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
              endif
   40   continue

        if(ifock.gt.0) call fock(uu(1),ss(1),tt(1),rri(1),rrj(1),it)

        fu=fu*dd1/rij
        fi=fi*dd7/ri
        fj=fj*dd8/rj
!!!  een for periodic systems         WAS
        if(icutjasc .gt.  0 .or. iperiodic .ne. 0) then
           fi = fi * fcut + (fc * fcutj *  dfcuti)/ri
           fj = fj * fcut + (fc * fcuti *  dfcutj)/rj
           fu = fu * fcut
           fc = fc * fcut
        endif
!!! end WAS
        fsn(i,j)=fsn(i,j) + fc

        fijn(1,i,j)=fijn(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        fijn(2,i,j)=fijn(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        fijn(3,i,j)=fijn(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        fijn(1,j,i)=fijn(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        fijn(2,j,i)=fijn(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        fijn(3,j,i)=fijn(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)

   50 continue

   55 fsumn=fsumn+fsn(i,j)-fso(i,j)
      fjn(1,i)=fjn(1,i)+fijn(1,i,j)-fijo(1,i,j)
      fjn(2,i)=fjn(2,i)+fijn(2,i,j)-fijo(2,i,j)
      fjn(3,i)=fjn(3,i)+fijn(3,i,j)-fijo(3,i,j)
      fjn(1,j)=fjn(1,j)+fijn(1,j,i)-fijo(1,j,i)
      fjn(2,j)=fjn(2,j)+fijn(2,j,i)-fijo(2,j,i)
      fjn(3,j)=fjn(3,j)+fijn(3,j,i)-fijo(3,j,i)
   60 continue

c e-n terms

   65 fijn(1,iel,iel)=0
      fijn(2,iel,iel)=0
      fijn(3,iel,iel)=0
      fsn(iel,iel)=0

      do 80 ic=1,ncent
        it=iwctype(ic)

        ri=r_en(iel,ic)
        if(ri.gt.cutjas_en) goto 80

        call scale_dist1(ri,rri(1),dd7,1)

        top=a4(1,it,iwf)*rri(1)
        topi=a4(1,it,iwf)

        bot=a4(2,it,iwf)*rri(1)
        boti=a4(2,it,iwf)

        bot=1+bot
        bot2=bot*bot
c       fen=top/bot-asymp_jasa(it,iwf)
c       feni=topi/bot-boti*top/bot2
        fen=top/bot
        feni=topi/bot2

        if(isc.eq.8 .or. isc.eq.10) then
          fen=fen/scalek(iwf)
          feni=feni/scalek(iwf)
        endif
        fen=fen-asymp_jasa(it,iwf)

        do 70 iord=2,norda
          rri(iord)=rri(1)**iord
          fen=fen+a4(iord+1,it,iwf)*rri(iord)
   70     feni=feni+a4(iord+1,it,iwf)*iord*rri(iord-1)

        feni=feni*dd7/ri

        fsn(iel,iel)=fsn(iel,iel)+fen

        fijn(1,iel,iel)=fijn(1,iel,iel) + feni*rvec_en(1,iel,ic)
        fijn(2,iel,iel)=fijn(2,iel,iel) + feni*rvec_en(2,iel,ic)
        fijn(3,iel,iel)=fijn(3,iel,iel) + feni*rvec_en(3,iel,ic)
   80 continue

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)
      fjn(1,iel)=fjn(1,iel)+fijn(1,iel,iel)-fijo(1,iel,iel)
      fjn(2,iel)=fjn(2,iel)+fijn(2,iel,iel)-fijo(2,iel,iel)
      fjn(3,iel)=fjn(3,iel)+fijn(3,iel,iel)-fijo(3,iel,iel)

      do 110 i=1,nelec
        v(1,i)=fjn(1,i)
        v(2,i)=fjn(2,i)
  110   v(3,i)=fjn(3,i)

      value=fsumn

      return
      end

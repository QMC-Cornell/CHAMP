      subroutine deriv_jastrow4(x,v,d2,value)
      use all_tools_mod
      use control_mod

c Written by Cyrus Umrigar and Claudia Filippi
c Jastrow 4,5 must be used with one of isc=2,4,6,7,12,14,16,17
c Jastrow 6   must be used with one of isc=6,7

      implicit real*8(a-h,o-z)

      parameter(NEQSX=6*MORDJ,MTERMS=55)
!JT      parameter (zero=0.d0,one=1.d0,two=2.d0)
!JT      parameter (half=.5d0,third=1.d0/3.d0,eps=1.d-12)
      parameter (eps=1.d-12)

      common /dim/ ndim
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
!!!   added WAS
      common /jas_c_cut/ cutjasc,icutjasc
!JT      common /contrl_per/ iperiodic,ibasis
!!!
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
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

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      common /jaso/ fso(MELEC,MELEC),fijo(3,MELEC,MELEC)
     &,d2ijo(MELEC,MELEC),d2o,fsumo,fjo(3,MELEC)

      common /contrl_opt2/ igradhess,iadd_diag_opt
      common /gradhessj_nonlin/ d1d2a(MCTYPE),d2d2a(MCTYPE),d1d2b(2),d2d2b(2),didk(MPARMJ)

c     common /jasnonloc/ fso(MELEC,MELEC),fsumo

      common /optim/ lo(MORB),npoint(MORB),
     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
     &imnbas(MCENT),
     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
     &necn,nebase
      common /pointer/ npointa(MPARMJ*NCTYP3X)
      common /derivjas/ gvalue(MPARMJ),g(3,MELEC,MPARMJ),d2g(MPARMJ)
     &,go(MELEC,MELEC,MPARMJ)

      common /cuspmat4/ d(NEQSX,MTERMS),iwc4(NEQSX),nterms
      common /vardep/ nvdepend(NEQSX,MCTYPE),iwdepend(NEQSX,MPARMJ,MCTYPE)
     &,cdep(NEQSX,MPARMJ,MCTYPE)

      dimension x(3,*),v(3,*)
      dimension uu(-3:MORDJ),ss(-3:MORDJ),tt(-3:MORDJ),rri(-3:MORDJ)
     &,rrj(-3:MORDJ)

      ndim1=ndim-1

      fsum=0

      do 2 iparm=1,nparmj+nparms
        gvalue(iparm)=0
        d2g(iparm)=0
        didk(iparm)=0
        do 2 i=1,nelec
          g(1,i,iparm)=0
          g(2,i,iparm)=0
   2      g(3,i,iparm)=0

      do 3 isb=1,nspin2b
        d1d2b(isb)=0
   3    d2d2b(isb)=0
      do 4 it=1,nctype
        d1d2a(it)=0
   4    d2d2a(it)=0

      do 5 i=-3,-1
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

      if(igradhess.gt.0) then
        iderivk=2
      else
        iderivk=1
      endif

      if(nelec.lt.2) goto 65

c e-e and e-e-n terms
      ij=0
      do 60 i=2,nelec
      im1=i-1
      do 60 j=1,im1
      ij=ij+1

      fso(i,j)=0
      d2ijo(i,j)=0
      do 10 k=1,ndim
        fijo(k,i,j)=0
   10   fijo(k,j,i)=0

      do 11 iparm=1,nparmj+nparms
   11   go(i,j,iparm)=0

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

      if(rij.gt.cutjas_ee) goto 30

      call scale_dist2(rij,uu(1),dd1,dd2,2)
      if(nparms.eq.1) call deriv_scale(rij,dk,dk2,dr,dr2,2,iderivk)


      top=sspinn*b(1,isb,iwf)*uu(1)
      topu=sspinn*b(1,isb,iwf)
      topuu=0

      bot=1+b(2,isb,iwf)*uu(1)
      botu=b(2,isb,iwf)
      botuu=0
      bot2=bot*bot


c      feeu=topu/bot-botu*top/bot2
c      feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
c      feeuu=feeuu/bot
c simpler expressions are :
      fee=top/bot
      feeu=topu/bot2
      feeuu=-2*feeu*botu/bot
      tempuu=feeuu
      if(isc.eq.8 .or. isc.eq.10) then
        fee=fee/scalek(iwf)
        feeu=feeu/scalek(iwf)
        feeuu=feeuu/scalek(iwf)
      endif
      fee=fee-asymp_jasb(ipar+1,iwf)


      do 20 iord=2,nordb
        uu(iord)=uu(1)*uu(iord-1)
        fee=fee+b(iord+1,isb,iwf)*uu(iord)
        feeu=feeu+b(iord+1,isb,iwf)*iord*uu(iord-1)
   20   feeuu=feeuu+b(iord+1,isb,iwf)*iord*(iord-1)*uu(iord-2)

c for scale derivatives we also need feeuuu:
      if(nparms.eq.1) then
        if(isc.eq.8 .or. isc.eq.10) then
          stop 'scalek opt not fully implemented for isc=8,10'
        endif
        feeuuu=-3*tempuu*botu/bot
        do 21 iord=3,nordb
   21     feeuuu=feeuuu+b(iord+1,isb,iwf)*iord*(iord-1)*(iord-2)*uu(iord-3)
      endif

c      feeuu=feeuu*dd1*dd1+feeu*dd2
c      feeu=feeu*dd1/rij
c we will need feeuu and feeu later
      tempuu=feeuu*dd1*dd1+feeu*dd2
      tempu=feeu*dd1/rij

      fso(i,j)=fso(i,j) + fee
      do 22 k=1,ndim
        fijo(k,i,j)= fijo(k,i,j) + tempu*rvec_ee(k,ij)
   22   fijo(k,j,i)= fijo(k,j,i) - tempu*rvec_ee(k,ij)
      d2ijo(i,j)=d2ijo(i,j)+2*(tempuu+ndim1*tempu)

c     v(1,i)=v(1,i) + tempu*rvec_ee(1,ij)
c     v(2,i)=v(2,i) + tempu*rvec_ee(2,ij)
c     v(3,i)=v(3,i) + tempu*rvec_ee(3,ij)
c     v(1,j)=v(1,j) - tempu*rvec_ee(1,ij)
c     v(2,j)=v(2,j) - tempu*rvec_ee(2,ij)
c     v(3,j)=v(3,j) - tempu*rvec_ee(3,ij)

c     d2 =d2 + two*(tempuu+ndim1*feeu)

c derivatives of wave function wrt b(1),b(2) and rest of b(i)

      iparma=nparma(1)
      do 23 it=2,nctype
   23  iparma=iparma+nparma(it)

      iparm0=iparma+nparms
      if(isb.eq.2) iparm0=iparm0+nparmb(1)
      do 25 jparm=1,nparmb(isb)
        iparm=iparm0+jparm

        if(iwjasb(jparm,isb).eq.1) then
          top=sspinn*uu(1)
          topu=sspinn
          topuu=zero

          bot=one+b(2,isb,iwf)*uu(1)
          botu=b(2,isb,iwf)
          botuu=zero
          bot2=bot*bot

          botasymp=1+b(2,isb,iwf)*asymp_r_ee(iwf)

          gee=top/bot-sspinn*asymp_r_ee(iwf)/botasymp
c         geeu=topu/bot-botu*top/bot2
c         geeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
c         geeuu=geeuu/bot
          geeu=topu/bot2
          geeuu=-2*botu*geeu/bot
          if(isc.eq.8 .or. isc.eq.10) then
            gee=gee/scalek(iwf)
            geeu=geeu/scalek(iwf)
            geeuu=geeuu/scalek(iwf)
          endif

          if(igradhess.gt.0 .and. nparms.eq.1) then
            if(isc.eq.8 .or. isc.eq.10) then
              stop 'scalek opt not fully implemented for isc=8,10'
            endif
            asympterm=sspinn/(botasymp*botasymp)*dasymp_r_ee(iwf)
c            asympterm=0
            didk(iparm)=didk(iparm)+geeu*dk-asympterm
          endif

         elseif(iwjasb(jparm,isb).eq.2) then
          top=-sspinn*b(1,isb,iwf)*uu(1)*uu(1)
          topu=-sspinn*2*b(1,isb,iwf)*uu(1)
          topuu=-sspinn*2*b(1,isb,iwf)

          bot0=one+b(2,isb,iwf)*uu(1)
          bot=bot0*bot0
          botu=2*bot0*b(2,isb,iwf)
          botuu=2*b(2,isb,iwf)*b(2,isb,iwf)
          bot2=bot*bot

          botasymp=1+b(2,isb,iwf)*asymp_r_ee(iwf)
          botasymp2=botasymp*botasymp

          gee=top/bot+sspinn*b(1,isb,iwf)*asymp_r_ee(iwf)**2/botasymp2
c          geeu=topu/bot-botu*top/bot2
c          geeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
c          geeuu=geeuu/bot
          geeu=topu/(bot*bot0)
          geeuu=topuu*(1-2*b(2,isb,iwf)*uu(1))/bot2
          if(isc.eq.8 .or. isc.eq.10) then
            gee=gee/scalek(iwf)
            geeu=geeu/scalek(iwf)
            geeuu=geeuu/scalek(iwf)
          endif

          if(igradhess.gt.0 .or. l_opt_jas_2nd_deriv) then
            if(isc.eq.8 .or. isc.eq.10) then
              d2d2b(isb)=d2d2b(isb)-
     &2*(top*uu(1)/(bot*bot0)+sspinn*b(1,isb,iwf)*asymp_r_ee(iwf)**3/(1+b(2,isb,iwf)*asymp_r_ee(iwf))**3)/scalek(iwf)
              d1d2b(isb)=d1d2b(isb)-
     &(sspinn*uu(1)*uu(1)/bot-sspinn*asymp_r_ee(iwf)**2/(1+b(2,isb,iwf)*asymp_r_ee(iwf))**2)/scalek(iwf)
            else
              d2d2b(isb)=d2d2b(isb)-
     &2*top*uu(1)/bot/bot0-2*sspinn*b(1,isb,iwf)*asymp_r_ee(iwf)**3/(1+b(2,isb,iwf)*asymp_r_ee(iwf))**3
              d1d2b(isb)=d1d2b(isb)-
     &sspinn*uu(1)*uu(1)/bot+sspinn*asymp_r_ee(iwf)**2/(1+b(2,isb,iwf)*asymp_r_ee(iwf))**2
            endif
            if(nparms.eq.1) then
              if(isc.eq.8 .or. isc.eq.10) then
                stop 'scalek opt not fully implemented for isc=8,10'
              endif
              asympterm=topuu*asymp_r_ee(iwf)/(botasymp*botasymp2)*dasymp_r_ee(iwf)
              didk(iparm)=didk(iparm)+geeu*dk-asympterm
            endif
          endif

         else
          iord=iwjasb(jparm,isb)-1
          if(ijas.eq.4) then
            asympiord=asymp_r_ee(iwf)**iord
            gee=uu(iord)-asympiord
            geeu=iord*uu(iord-1)
            geeuu=iord*(iord-1)*uu(iord-2)
            if(nparms.eq.1 .and. igradhess.gt.0) then
              asympterm=asympiord/asymp_r_ee(iwf)*dasymp_r_ee(iwf)
c              asympterm=0
              didk(iparm)=didk(iparm)+iord*(uu(iord-1)*dk-asympterm)
            endif
           elseif(ijas.eq.5) then
            gee=sspinn*(uu(iord)-asymp_r_ee(iwf)**iord)
            geeu=sspinn*iord*uu(iord-1)
            geeuu=sspinn*iord*(iord-1)*uu(iord-2)
          endif

        endif

        geeuu=geeuu*dd1*dd1+geeu*dd2
        geeu=geeu*dd1/rij

        go(i,j,iparm)=go(i,j,iparm)+gee
        gvalue(iparm)=gvalue(iparm)+gee

        g(1,i,iparm)=g(1,i,iparm) + geeu*rvec_ee(1,ij)
        g(2,i,iparm)=g(2,i,iparm) + geeu*rvec_ee(2,ij)
        g(3,i,iparm)=g(3,i,iparm) + geeu*rvec_ee(3,ij)
        g(1,j,iparm)=g(1,j,iparm) - geeu*rvec_ee(1,ij)
        g(2,j,iparm)=g(2,j,iparm) - geeu*rvec_ee(2,ij)
        g(3,j,iparm)=g(3,j,iparm) - geeu*rvec_ee(3,ij)

        d2g(iparm)=d2g(iparm) + two*(geeuu+ndim1*geeu)

   25 continue

c derivatives (go,gvalue and g) wrt scalek parameter
      if(nparms.eq.1) then
        if(isc.eq.8 .or. isc.eq.10) then
          stop 'scalek opt not fully implemented for isc=8,10'
        endif
        iparm=1
c        call deriv_scale(rij,dk,dk2,dr,dr2,2,iderivk)

        gee=feeu*dk-dasymp_jasb(ipar+1,iwf)*dasymp_r_ee(iwf)
        geeu=(feeuu*dk*dd1+feeu*dr)/rij
        geeuu=feeuu*dk*dd2+(feeuuu*dk*dd1+2*feeuu*dr)*dd1+feeu*dr2

        go(i,j,iparm)=go(i,j,iparm)+gee
        gvalue(iparm)=gvalue(iparm)+gee

        g(1,i,iparm)=g(1,i,iparm) + geeu*rvec_ee(1,ij)
        g(2,i,iparm)=g(2,i,iparm) + geeu*rvec_ee(2,ij)
        g(3,i,iparm)=g(3,i,iparm) + geeu*rvec_ee(3,ij)
        g(1,j,iparm)=g(1,j,iparm) - geeu*rvec_ee(1,ij)
        g(2,j,iparm)=g(2,j,iparm) - geeu*rvec_ee(2,ij)
        g(3,j,iparm)=g(3,j,iparm) - geeu*rvec_ee(3,ij)

        d2g(iparm)=d2g(iparm) + two*(geeuu+ndim1*geeu)

        if(igradhess.gt.0) then
          didk(1)=didk(1)+feeuu*dk*dk+feeu*dk2
     &-d2asymp_jasb(ipar+1,iwf)*dasymp_r_ee(iwf)*dasymp_r_ee(iwf)-dasymp_jasb(ipar+1,iwf)*d2asymp_r_ee(iwf)
        endif

      endif

c There are no C terms to order 1.
   30 if(nordc.le.1) goto 58

c     if(isc.ge.12) call scale_dist2(rij,uu(1),dd1,dd2,3)
      call scale_dist2(rij,uu(1),dd1,dd2,4)
      if(ijas.eq.4.or.ijas.eq.5) then
        call switch_scale2(uu(1),dd1,dd2,4)
        do 35 iord=2,nordc
   35     uu(iord)=uu(1)*uu(iord-1)
      endif

      do 57 ic=1,ncent
        it=iwctype(ic)

        ri=r_en(i,ic)
        rj=r_en(j,ic)

        if(ri.gt.cutjas_en .or. rj.gt.cutjas_en) goto 57
        do 37 k=1,ndim
   37     if(abs(rshift(k,i,ic)-rshift(k,j,ic)).gt.eps) goto 57

        call scale_dist2(ri,rri(1),dd7,dd9,3)
        call scale_dist2(rj,rrj(1),dd8,dd10,3)

        if(ijas.eq.4.or.ijas.eq.5) then
          call switch_scale2(rri(1),dd7,dd9,3)
          call switch_scale2(rrj(1),dd8,dd10,3)
        endif

!!!! WAS
        if(icutjasc .gt. 0 .or. iperiodic .ne. 0) then
c           call f_een_cuts (cutjasc, ri, rj, fcuti, fcutj, fcut,
           call f_een_cuts (cutjas_en, ri, rj, fcuti, fcutj, fcut,
     +          dfcuti, dfcutj, d2fcuti, d2fcutj)
        endif
!!!

        if(nparms.eq.1) then
          call deriv_scale(rij,dkij,dk2ij,drij,dr2ij,4,iderivk)
          call deriv_scale(ri,dki,dk2i,dri,dr2i,3,iderivk)
          call deriv_scale(rj,dkj,dk2j,drj,dr2j,3,iderivk)
          call switch_dscale(uu(1),dd1,dd2,dkij,dk2ij,drij,dr2ij,4,iderivk)
          call switch_dscale(rri(1),dd7,dd9,dki,dk2i,dri,dr2i,3,iderivk)
          call switch_dscale(rrj(1),dd8,dd10,dkj,dk2j,drj,dr2j,3,iderivk)
        endif

        s=ri+rj
        t=ri-rj
c       u2mt2=rij*rij-t*t
        u2pst=rij*rij+s*t
        u2mst=rij*rij-s*t
c       s2mu2=s*s-rij*rij
c       s2mt2=s*s-t*t

        do 50 iord=1,nordc
c         rri(iord)=rri(1)**iord
c         rrj(iord)=rrj(1)**iord
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
        fij=0
        fuuu=0
        fuui=0
        fuuj=0
        fuii=0
        fuij=0
        fujj=0
        fiii=0
        fiij=0
        fijj=0
        fjjj=0

        ll=0
        jj=1
        jparm=1
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
                p1=uu(k)
                p1u=k*uu(k-1)
                p1uu=k*(k-1)*uu(k-2)
                p2=ss(l)*tt(m)
                p2i=((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                p2ii=((l+m)*(l+m-1)*rri(l+m-2)*rrj(m)+m*(m-1)*rri(m-2)*rrj(l+m))
                p2j=((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
                p2jj=((l+m)*(l+m-1)*rrj(l+m-2)*rri(m)+m*(m-1)*rrj(m-2)*rri(l+m))
                pc=p1*p2
                pu=p1u*p2
                puu=p1uu*p2
                ppi=p1*p2i
                pii=p1*p2ii
                pj=p1*p2j
                pjj=p1*p2jj
                pui=p1u*p2i
                puj=p1u*p2j

c                pc=uu(k)*ss(l)*tt(m)
c                pu=k*uu(k-1)*ss(l)*tt(m)
c                puu=k*(k-1)*uu(k-2)*ss(l)*tt(m)
c                ppi=uu(k)
c     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
c                pii=uu(k)
c     &          *((l+m)*(l+m-1)*rri(l+m-2)*rrj(m)
c     &          +m*(m-1)*rri(m-2)*rrj(l+m))
c                pj=uu(k)
c     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
c                pjj=uu(k)
c     &          *((l+m)*(l+m-1)*rrj(l+m-2)*rri(m)
c     &          +m*(m-1)*rrj(m-2)*rri(l+m))
c                pui=k*uu(k-1)
c     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
c                puj=k*uu(k-1)
c     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))

                fc=fc+c(ll,it,iwf)*pc
                fu=fu+c(ll,it,iwf)*pu
                fuu=fuu+c(ll,it,iwf)*puu
                fi=fi+c(ll,it,iwf)*ppi
                fii=fii+c(ll,it,iwf)*pii
                fj=fj+c(ll,it,iwf)*pj
                fjj=fjj+c(ll,it,iwf)*pjj
                fui=fui+c(ll,it,iwf)*pui
                fuj=fuj+c(ll,it,iwf)*puj

c quantities needed for scale derivatives:
                if(nparms.eq.1) then

                  p1uuu=k*(k-1)*(k-2)*uu(k-3)
                  p2ij =((l+m)*rri(l+m-1)*m*rrj(m-1)             + m*(l+m)*rri(m-1)*rrj(l+m-1))
                  p2iii=((l+m)*(l+m-1)*(l+m-2)*rri(l+m-3)*rrj(m) + m*(m-1)*(m-2)*rri(m-3)*rrj(l+m))
                  p2iij=((l+m)*(l+m-1)*m*rri(l+m-2)*rrj(m-1)     + m*(m-1)*(l+m)*rri(m-2)*rrj(l+m-1))
                  p2ijj=((l+m)*rri(l+m-1)*m*(m-1)*rrj(m-2)       + m*(l+m)*(l+m-1)*rri(m-1)*rrj(l+m-2))
                  p2jjj=((l+m)*(l+m-1)*(l+m-2)*rrj(l+m-3)*rri(m) + m*(m-1)*(m-2)*rrj(m-3)*rri(l+m))

                  fij=fij+c(ll,it,iwf)*p1*p2ij
                  fuuu=fuuu+c(ll,it,iwf)*p1uuu*p2
                  fuui=fuui+c(ll,it,iwf)*p1uu*p2i
                  fuuj=fuuj+c(ll,it,iwf)*p1uu*p2j
                  fuii=fuii+c(ll,it,iwf)*p1u*p2ii
                  fuij=fuij+c(ll,it,iwf)*p1u*p2ij
                  fujj=fujj+c(ll,it,iwf)*p1u*p2jj
                  fiii=fiii+c(ll,it,iwf)*p1*p2iii
                  fiij=fiij+c(ll,it,iwf)*p1*p2iij
                  fijj=fijj+c(ll,it,iwf)*p1*p2ijj
                  fjjj=fjjj+c(ll,it,iwf)*p1*p2jjj

                endif

c derivatives of wave function wrt c-parameters
c ideriv = 0 parameter is not varied and is not dependent
c        = 1 parameter is a dependent parameter
c        = 2 parameter is an independent parameter that is varied
                ideriv=0
                if(ll.eq.iwjasc(jparm,it)) then
                  ideriv=2
                 else
                  do 31 id=1,2*(nordc-1)           !   dowhile loop would be more efficient here?
                    if(ll.eq.iwc4(id)) then
                      jj=id
                      if(nvdepend(jj,it).gt.0) ideriv=1
                    endif
   31             continue
                endif

                if(ideriv.gt.0) then
                  gp=pc
                  gu=pu
                  guu=puu
                  gi=ppi
                  gii=pii
                  gj=pj
                  gjj=pjj
                  gui=pui
                  guj=puj

                  guu=guu*dd1*dd1+gu*dd2
                  gu=gu*dd1/rij

                  gui=gui*dd1*dd7
                  guj=guj*dd1*dd8

                  gii=gii*dd7*dd7+gi*dd9
                  gjj=gjj*dd8*dd8+gj*dd10
                  gi=gi*dd7/ri
                  gj=gj*dd8/rj

!!!!  een for periodic systems         WAS
                  if(icutjasc .gt. 0 .or. iperiodic .ne. 0) then

                     guu = guu * fcut
                     gii = gii * fcut +(2 * gi * ri * dfcuti + gp * d2fcuti)*fcutj
                     gi = gi * fcut + (gp * fcutj *  dfcuti)/ri
                     gui = gui * fcut + (gu * fcutj *  dfcuti*rij)

                     gjj = gjj * fcut + (2 * gj * dfcutj *rj + gp * d2fcutj)*fcuti
                     gj = gj * fcut + (gp * fcuti *  dfcutj)/rj
                     guj = guj * fcut + (gu * fcuti *  dfcutj * rij)
                     gp = gp * fcut
                     gu = gu * fcut

                  endif
!!! end WAS

                  if(ideriv.eq.1) then

                    do 33 id=1,nvdepend(jj,it)

                      iparm=npoint(it)+iwdepend(jj,id,it)+nparms
                      cd=cdep(jj,id,it)

                      go(i,j,iparm)=go(i,j,iparm)+cd*gp
                      gvalue(iparm)=gvalue(iparm)+cd*gp

                      g(1,i,iparm)=g(1,i,iparm)+cd*(gi*rvec_en(1,i,ic)+gu*rvec_ee(1,ij))
                      g(2,i,iparm)=g(2,i,iparm)+cd*(gi*rvec_en(2,i,ic)+gu*rvec_ee(2,ij))
                      g(3,i,iparm)=g(3,i,iparm)+cd*(gi*rvec_en(3,i,ic)+gu*rvec_ee(3,ij))
                      g(1,j,iparm)=g(1,j,iparm)+cd*(gj*rvec_en(1,j,ic)-gu*rvec_ee(1,ij))
                      g(2,j,iparm)=g(2,j,iparm)+cd*(gj*rvec_en(2,j,ic)-gu*rvec_ee(2,ij))
                      g(3,j,iparm)=g(3,j,iparm)+cd*(gj*rvec_en(3,j,ic)-gu*rvec_ee(3,ij))

   33                 d2g(iparm)=d2g(iparm) + cd*((ndim-1)*(2*gu+gi+gj)
     &                + 2*guu + gii +  gjj + gui*u2pst/(ri*rij) + guj*u2mst/(rj*rij))

c  33                 d2g(iparm)=d2g(iparm) + cd*(2*(guu + 2*gu)
c    &                + gui*u2pst/(ri*rij) + guj*u2mst/(rj*rij)
c    &                + gii + 2*gi + gjj + 2*gj)

c                   jj=jj+1

                    if(igradhess.gt.0 .and. nparms.eq.1) then
                      didk(iparm)=didk(iparm)+cd*(pu*dkij+ppi*dki+pj*dkj)
                    endif

                   elseif(ideriv.eq.2) then

                    iparm=npoint(it)+jparm+nparms

                    go(i,j,iparm)=go(i,j,iparm)+gp
                    gvalue(iparm)=gvalue(iparm)+gp

                    g(1,i,iparm)=g(1,i,iparm)+gi*rvec_en(1,i,ic)+gu*rvec_ee(1,ij)
                    g(2,i,iparm)=g(2,i,iparm)+gi*rvec_en(2,i,ic)+gu*rvec_ee(2,ij)
                    g(3,i,iparm)=g(3,i,iparm)+gi*rvec_en(3,i,ic)+gu*rvec_ee(3,ij)
                    g(1,j,iparm)=g(1,j,iparm)+gj*rvec_en(1,j,ic)-gu*rvec_ee(1,ij)
                    g(2,j,iparm)=g(2,j,iparm)+gj*rvec_en(2,j,ic)-gu*rvec_ee(2,ij)
                    g(3,j,iparm)=g(3,j,iparm)+gj*rvec_en(3,j,ic)-gu*rvec_ee(3,ij)
                    d2g(iparm)=d2g(iparm) + (ndim-1)*(2*gu+gi+gj)
     &              + 2*guu + gii +  gjj + gui*u2pst/(ri*rij) + guj*u2mst/(rj*rij)

c                   d2g(iparm)=d2g(iparm) + 2*(guu + 2*gu)
c    &              + gui*u2pst/(ri*rij) + guj*u2mst/(rj*rij)
c    &              + gii + 2*gi + gjj + 2*gj

                    jparm=jparm+1

                    if(igradhess.gt.0 .and. nparms.eq.1) then
                      didk(iparm)=didk(iparm)+pu*dkij+ppi*dki+pj*dkj
                    endif

                  endif
c               write(6,'(''i,j,iparm,go(i,j,iparm),gvalue(iparm),(g(k,i,iparm),g(k,j,iparm),k=1,3)'',3i5,12d12.4)')
c    &          i,j,iparm,go(i,j,iparm),gvalue(iparm),(g(kk,i,iparm),g(kk,j,iparm),kk=1,3)

                endif
              endif
   55   continue

c derivatives (go,gvalue and g) wrt scalek parameter
        if(nparms.eq.1) then

          iparm=1

          gp=fu*dkij + fi*dki + fj*dkj
          gu=(dkij*fuu*dd1 + fu*drij + dki*fui*dd1 + dkj*fuj*dd1)/rij
          gi=(dkij*fui*dd7 + dki*fii*dd7 + fi*dri + dkj*fij*dd7)/ri
          gj=(dkij*fuj*dd8 + dkj*fjj*dd8 + fj*drj + dki*fij*dd8)/rj

          go(i,j,iparm)=go(i,j,iparm) + gp
          gvalue(iparm)=gvalue(iparm) + gp

          g(1,i,iparm)=g(1,i,iparm) + gi*rvec_en(1,i,ic) + gu*rvec_ee(1,ij)
          g(2,i,iparm)=g(2,i,iparm) + gi*rvec_en(2,i,ic) + gu*rvec_ee(2,ij)
          g(3,i,iparm)=g(3,i,iparm) + gi*rvec_en(3,i,ic) + gu*rvec_ee(3,ij)
          g(1,j,iparm)=g(1,j,iparm) + gj*rvec_en(1,j,ic) - gu*rvec_ee(1,ij)
          g(2,j,iparm)=g(2,j,iparm) + gj*rvec_en(2,j,ic) - gu*rvec_ee(2,ij)
          g(3,j,iparm)=g(3,j,iparm) + gj*rvec_en(3,j,ic) - gu*rvec_ee(3,ij)

          gui=dkij*fuui*dd1*dd7 + drij*fui*dd7 + dri*fui*dd1 + dki*fuii*dd1*dd7 + dkj*fuij*dd1*dd7
          guj=dkij*fuuj*dd1*dd8 + drij*fuj*dd8 + drj*fuj*dd1 + dkj*fujj*dd1*dd8 + dki*fuij*dd1*dd8
          gtu=dd2*(fuu*dkij + fui*dki + fuj*dkj) + 2*fuu*dd1*drij + fu*dr2ij
     &      + dd1*dd1*(fuuu*dkij + fuui*dki + fuuj*dkj)
          gti=dd9*(fui*dkij + fii*dki + fij*dkj) + fi*dr2i + 2*fii*dd7*dri
     &      + dd7*dd7*(fuii*dkij + fiii*dki + fiij*dkj)
          gtj=dd10*(fuj*dkij + fjj*dkj + fij*dki) + fj*dr2j + 2*fjj*dd8*drj
     &       + dd8*dd8*(fujj*dkij + fjjj*dkj + fijj*dki)


          d2g(iparm)=d2g(iparm) + (2*gu+gi+gj)*ndim1 + gui*u2pst/(ri*rij) + guj*u2mst/(rj*rij)
     &              + 2*gtu + gti + gtj

          if(igradhess.gt.0 .and. nparms.eq.1) then
            didk(iparm)=didk(iparm)+(fuu*dkij+fui*dki+fuj*dkj)*dkij+fu*dk2ij
     &                             +(fui*dkij+fii*dki+fij*dkj)*dki+fi*dk2i
     &                             +(fuj*dkij+fij*dki+fjj*dkj)*dkj+fj*dk2j
          endif

        endif

        if(ifock.gt.0) call fock(uu(1),ss(1),tt(1),rri(1),rrj(1),it)

        fuu=fuu*dd1*dd1+fu*dd2
        fu=fu*dd1/rij

        fui=fui*dd1*dd7
        fuj=fuj*dd1*dd8

        fii=fii*dd7*dd7+fi*dd9
        fjj=fjj*dd8*dd8+fj*dd10
        fi=fi*dd7/ri
        fj=fj*dd8/rj


!!!!  een for periodic systems         WAS
       if(icutjasc .gt. 0 .or. iperiodic .ne. 0) then

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
!!!! end WAS

        fso(i,j)=fso(i,j) + fc

        fijo(1,i,j)=fijo(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        fijo(2,i,j)=fijo(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        fijo(3,i,j)=fijo(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        fijo(1,j,i)=fijo(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        fijo(2,j,i)=fijo(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        fijo(3,j,i)=fijo(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)
c       write(6,'(''i,j,fijo2='',2i5,9d12.4)') i,j,(fijo(k,i,j),k=1,ndim)

c       d2ijo(i,j)=d2ijo(i,j) + 2*(fuu + 2*fu) + fui*u2pst/(ri*rij)
c    &  + fuj*u2mst/(rj*rij) + fii + 2*fi + fjj + 2*fj
        d2ijo(i,j)=d2ijo(i,j) + ndim1*(2*fu+fi+fj)
     &  + 2*fuu + fii +  fjj + fui*u2pst/(ri*rij) + fuj*u2mst/(rj*rij)

c       v(1,i)=v(1,i) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
c       v(2,i)=v(2,i) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
c       v(3,i)=v(3,i) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
c       v(1,j)=v(1,j) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
c       v(2,j)=v(2,j) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
c       v(3,j)=v(3,j) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)

c       d2 = d2 + ndim1*(2*fu+fi+fj)
c    &  + 2*fuu + fii +  fjj + fui*u2pst/(ri*rij) + fuj*u2mst/(rj*rij)

cc      d2 = d2 + 2*(fuu + 2*fu) + fui*u2pst/(ri*rij)
cc   &  + fuj*u2mst/(rj*rij) + fii + 2*fi + fjj + 2*fj

   57 continue

   58 fsum=fsum+fso(i,j)
      v(1,i)=v(1,i)+fijo(1,i,j)
      v(2,i)=v(2,i)+fijo(2,i,j)
      v(3,i)=v(3,i)+fijo(3,i,j)
      v(1,j)=v(1,j)+fijo(1,j,i)
      v(2,j)=v(2,j)+fijo(2,j,i)
      v(3,j)=v(3,j)+fijo(3,j,i)
c  write(6,'(''v='',2i2,9d12.4)') i,j,(v(k,i),v(k,j),k=1,3)
c     div_vj(i)=div_vj(i)+d2ijo(i,j)/2
c     div_vj(j)=div_vj(j)+d2ijo(i,j)/2
   60 d2=d2+d2ijo(i,j)
c     write(6,'(''fsum,d2='',9d12.4)') fsum,d2

c e-n terms
   65 do 90 i=1,nelec

        fso(i,i)=0
        fijo(1,i,i)=0
        fijo(2,i,i)=0
        fijo(3,i,i)=0
        d2ijo(i,i)=0

        do 66 iparm=1,nparmj+nparms
   66     go(i,i,iparm)=zero

        do 80 ic=1,ncent
          it=iwctype(ic)

          ri=r_en(i,ic)
          if(ri.gt.cutjas_en) goto 80

          call scale_dist2(ri,rri(1),dd7,dd9,1)
          if(nparms.eq.1) call deriv_scale(ri,dk,dk2,dr,dr2,1,iderivk)

          top=a4(1,it,iwf)*rri(1)
          topi=a4(1,it,iwf)
          topii=0

          bot=1+a4(2,it,iwf)*rri(1)
          boti=a4(2,it,iwf)
          botii=0
          bot2=bot*bot

c          feni=topi/bot-boti*top/bot2
c          fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
c          fenii=fenii/bot
c simpler expressions are :
          fen=top/bot
          feni=topi/bot2
          fenii=-2*feni*boti/bot
          tempii=fenii
          if(isc.eq.8 .or. isc.eq.10) then
            fen=fen/scalek(iwf)
            feni=feni/scalek(iwf)
            fenii=fenii/scalek(iwf)
          endif
          fen=fen-asymp_jasa(it,iwf)

          do 70 iord=2,norda
            rri(iord)=rri(1)*rri(iord-1)
            fen=fen+a4(iord+1,it,iwf)*rri(iord)
            feni=feni+a4(iord+1,it,iwf)*iord*rri(iord-1)
   70       fenii=fenii+a4(iord+1,it,iwf)*iord*(iord-1)*rri(iord-2)

c for scale derivatives we also need feniii:
          if(nparms.eq.1) then
            if(isc.eq.8 .or. isc.eq.10) then
              stop 'scalek opt not fully implemented for isc=8,10'
            endif
            feniii=-3*tempii*boti/bot
            do 73 iord=3,norda
   73         feniii=feniii+a4(iord+1,it,iwf)*iord*(iord-1)*(iord-2)*rri(iord-3)
          endif

          tempii=fenii*dd7*dd7+feni*dd9
          tempi=feni*dd7/ri

          fso(i,i)=fso(i,i)+fen

          fijo(1,i,i)=fijo(1,i,i) + tempi*rvec_en(1,i,ic)
          fijo(2,i,i)=fijo(2,i,i) + tempi*rvec_en(2,i,ic)
          fijo(3,i,i)=fijo(3,i,i) + tempi*rvec_en(3,i,ic)
c         write(6,'(''fijo='',9d12.4)') (fijo(k,i,i),k=1,ndim),feni,rvec_en(1,i,ic)

          d2ijo(i,i) = d2ijo(i,i) + tempii + ndim1*tempi

c         v(1,i)=v(1,i) + feni*rvec_en(1,i,ic)
c         v(2,i)=v(2,i) + feni*rvec_en(2,i,ic)
c         v(3,i)=v(3,i) + feni*rvec_en(3,i,ic)

c         d2 = d2 + fenii + ndim1*feni

          do 78 jparm=1,nparma(it)
            iparm=npointa(it)+jparm+nparms

            if(iwjasa(jparm,it).eq.1) then
              top=rri(1)
              topi=one
              topii=zero

              bot=one+a4(2,it,iwf)*rri(1)
              boti=a4(2,it,iwf)
              botii=zero
              bot2=bot*bot

              botasymp=1+a4(2,it,iwf)*asymp_r_en(iwf)

              gen=top/bot-asymp_r_en(iwf)/botasymp
c              geni=topi/bot-boti*top/bot2
c              genii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
c              genii=genii/bot
              geni=topi/bot2
              genii=-2*boti*geni/bot
              if(isc.eq.8 .or. isc.eq.10) then
                gen=gen/scalek(iwf)
                geni=geni/scalek(iwf)
                genii=genii/scalek(iwf)
              endif

              if(igradhess.gt.0 .and. nparms.eq.1) then
                if(isc.eq.8 .or. isc.eq.10) then
                  stop 'scalek opt not fully implemented for isc=8,10'
                endif
                asympterm=dasymp_r_en(iwf)/(botasymp*botasymp)
c                asympterm=0
                didk(iparm)=didk(iparm)+geni*dk-asympterm
              endif

             elseif(iwjasa(jparm,it).eq.2) then
              top=-a4(1,it,iwf)*rri(1)*rri(1)
              topi=-2*a4(1,it,iwf)*rri(1)
              topii=-2*a4(1,it,iwf)

              bot0=one+a4(2,it,iwf)*rri(1)
              bot=bot0*bot0
              boti=2*bot0*a4(2,it,iwf)
              botii=2*a4(2,it,iwf)*a4(2,it,iwf)
              bot2=bot*bot

              botasymp=1+a4(2,it,iwf)*asymp_r_en(iwf)
              botasymp2=botasymp*botasymp

              gen=top/bot+a4(1,it,iwf)*asymp_r_en(iwf)**2/botasymp2
c             geni=topi/bot-boti*top/bot2
c             genii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
c             genii=genii/bot
              geni=topi/(bot*bot0)
              genii=topii*(1-2*a4(2,it,iwf)*rri(1))/bot2
              if(isc.eq.8 .or. isc.eq.10) then
                gen=gen/scalek(iwf)
                geni=geni/scalek(iwf)
                genii=genii/scalek(iwf)
              endif

              if(igradhess.gt.0 .or. l_opt_jas_2nd_deriv) then
                d2d2a(it)=d2d2a(it)-2*top*rri(1)/bot/bot0-2*a4(1,it,iwf)*asymp_r_en(iwf)**3/(1+a4(2,it,iwf)*asymp_r_en(iwf))**3
                d1d2a(it)=d1d2a(it)-rri(1)*rri(1)/bot+asymp_r_en(iwf)**2/(1+a4(2,it,iwf)*asymp_r_en(iwf))**2
                if(isc.eq.8 .or. isc.eq.10) then
                  d2d2a(it)=d2d2a(it)/scalek(iwf)
                  d1d2a(it)=d1d2a(it)/scalek(iwf)
                endif
                if(ipr.ge.1)
     & write(6,'(''it,rri(1),a4(1,it,iwf),asymp_r_en(iwf),a4(2,it,iwf),top,bot,bot0,d1d2a(it),d2d2a(it)'',i2,20d12.4)')
     & it,rri(1),a4(1,it,iwf),asymp_r_en(iwf),a4(2,it,iwf),top,bot,bot0,d1d2a(it),d2d2a(it)
                if(nparms.eq.1) then
                  if(isc.eq.8 .or. isc.eq.10) then
                    stop 'scalek opt not fully implemented for isc=8,10'
                  endif
                  asympterm=topii*asymp_r_en(iwf)/(botasymp*botasymp2)*dasymp_r_en(iwf)
c                  asympterm=0
                  didk(iparm)=didk(iparm)+geni*dk-asympterm
                endif
              endif

             else
              iord=iwjasa(jparm,it)-1
              asympiord=asymp_r_en(iwf)**iord
              gen=rri(iord)-asympiord
              geni=iord*rri(iord-1)
              genii=iord*(iord-1)*rri(iord-2)
              if(nparms.eq.1 .and. igradhess.gt.0) then
                asympterm=asympiord/asymp_r_en(iwf)*dasymp_r_en(iwf)
c                asympterm=0
                didk(iparm)=didk(iparm)+iord*(rri(iord-1)*dk-asympterm)
              endif
            endif

            genii=genii*dd7*dd7+geni*dd9
            geni=geni*dd7/r_en(i,ic)

            go(i,i,iparm)=go(i,i,iparm)+gen
            gvalue(iparm)=gvalue(iparm)+gen

            g(1,i,iparm)=g(1,i,iparm)+geni*rvec_en(1,i,ic)
            g(2,i,iparm)=g(2,i,iparm)+geni*rvec_en(2,i,ic)
            g(3,i,iparm)=g(3,i,iparm)+geni*rvec_en(3,i,ic)

   78       d2g(iparm)=d2g(iparm)+genii+ndim1*geni

c derivatives (go,gvalue and g) wrt scalek parameter
          if(nparms.eq.1) then
            if(isc.eq.8 .or. isc.eq.10) then
              stop 'scalek opt not fully implemented for isc=8,10'
            endif

            iparm=1
c            call deriv_scale(ri,dk,dk2,dr,dr2,1,iderivk)

            gen=feni*dk-dasymp_jasa(it,iwf)*dasymp_r_en(iwf)
            geni=(fenii*dk*dd7+feni*dr)/ri
            genii=fenii*dk*dd9+(feniii*dk*dd7+2*fenii*dr)*dd7+feni*dr2

            go(i,i,iparm)=go(i,i,iparm)+gen
            gvalue(iparm)=gvalue(iparm)+gen

            g(1,i,iparm)=g(1,i,iparm) + geni*rvec_en(1,i,ic)
            g(2,i,iparm)=g(2,i,iparm) + geni*rvec_en(2,i,ic)
            g(3,i,iparm)=g(3,i,iparm) + geni*rvec_en(3,i,ic)

            d2g(iparm)=d2g(iparm) + genii+ndim1*geni

            if(igradhess.gt.0) then
              didk(1)=didk(1)+fenii*dk*dk+feni*dk2
     &-d2asymp_jasa(it,iwf)*dasymp_r_en(iwf)*dasymp_r_en(iwf)-dasymp_jasa(it,iwf)*d2asymp_r_en(iwf)
            endif

          endif

   80     continue

        fsum=fsum+fso(i,i)
        v(1,i)=v(1,i)+fijo(1,i,i)
        v(2,i)=v(2,i)+fijo(2,i,i)
        v(3,i)=v(3,i)+fijo(3,i,i)
c       write(6,'(''v='',9d12.4)') (v(k,i),k=1,ndim)
c       div_vj(i)=div_vj(i)+d2ijo(i,i)
   90   d2=d2+d2ijo(i,i)

      if(ijas.eq.6) then
        term=1/(c1_jas6_en(iwf)*scalek(iwf))
        fsum=term*fsum
        d2=term*d2
        do 100 i=1,nelec
c         div_vj(i)=term*div_vj(i)
          do 95 k=1,ndim
   95       v(k,i)=term*v(k,i)
          do 100 j=1,nelec
            d2ijo(i,j)=term*d2ijo(i,j)
            do 100 k=1,ndim
  100         fijo(k,i,j)=term*fijo(k,i,j)
      endif

      fsumo=fsum
      d2o=d2
      do 110 i=1,nelec
        fjo(1,i)=v(1,i)
        fjo(2,i)=v(2,i)
  110   fjo(3,i)=v(3,i)

      value=fsum

      call object_modified_by_index (gvalue_index) !JT
      call object_modified_by_index (g_index) !JT

      if (l_opt_jas_2nd_deriv) then
       call object_modified_by_index (d2d2a_index) !JT
       call object_modified_by_index (d1d2a_index) !JT
       call object_modified_by_index (d1d2b_index) !JT
       call object_modified_by_index (d2d2b_index) !JT
      endif

      return
      end

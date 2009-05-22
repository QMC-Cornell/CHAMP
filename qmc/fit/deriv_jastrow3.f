      subroutine deriv_jastrow3(x,v,d2,value)
c Written by Claudia Filippi
      use atom_mod
      use dets_mod
      use optim_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use wfsec_mod
      use derivjas_mod
      implicit real*8(a-h,o-z)

!JT      include '../vmc/vmc.h'
!JT      include '../vmc/force.h'
!JT      include 'fit.h'

!JT      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
!JT      parameter (half=.5d0,third=1.d0/3.d0)
      parameter(NEQSX=6*MORDJ)

!JT      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
!JT     &,ifock,i3body,irewgt,iaver,istrch
!JT     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord

      common /cuspmat/ cm(NEQSX,NEQSX),iwc(NEQSX),neqs,ishe
      common /vardep/ nvdepend(NEQSX,MCTYPE),iwdepend(NEQSX,MPARMJ,MCTYPE)
     &,cdep(NEQSX,MPARMJ,MCTYPE)

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /bparm/ nspin2b,nocuspb
      common /focktmp/ fc,fu,fuu,fs,fss,ft,ftt,fst,fus,fut
      common /jaso/ fso(MELEC,MELEC),fijo(3,MELEC,MELEC)
     &,d2ijo(MELEC,MELEC),d2o,fsumo,fjo(3,MELEC)

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase
!JT      common /derivjas/ gvalue(MPARMJ),g(3,MELEC,MPARMJ),d2g(MPARMJ)
!JT     &,go(MELEC,MELEC,MPARMJ)
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      dimension x(3,*),v(3,*)
      dimension u(-2:MORDJ),s(-2:MORDJ),t(-2:MORDJ)

      ndim1=ndim-1

      fsum=zero
      u(-2)=zero
      s(-2)=zero
      t(-2)=zero
      u(-1)=zero
      s(-1)=zero
      t(-1)=zero
      u(0)=one
      s(0)=one
      t(0)=one

      do 10 iparm=1,nparmj
        gvalue(iparm)=zero
        d2g(iparm)=zero
        do 10 i=1,nelec
          g(1,i,iparm)=zero
          g(2,i,iparm)=zero
  10      g(3,i,iparm)=zero

      if(nelec.lt.2) goto 47

      ncentt=ncent

      ij=0
      do 45 i=2,nelec
      im1=i-1
      do 45 j=1,im1
      ij=ij+1

      sspinn=1
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
      endif

      fso(i,j)=zero
      do 20 iparm=1,nparmj
   20   go(i,j,iparm)=zero

      rij=r_ee(ij)

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2) then
          u(1)=(one-dexp(-scalek(iwf)*rij))/scalek(iwf)
          dd1=one-scalek(iwf)*u(1)
          dd2=-scalek(iwf)*dd1
	 elseif(isc.eq.3) then
	  rijs=scalek(iwf)*rij
          rijs2=rijs*rijs
          rijs3=rijs*rijs2
	  exprij=dexp(-rijs-half*rijs2-third*rijs3)
	  u(1)=(one-exprij)/scalek(iwf)
	  dd1=(one+rijs+rijs2)*exprij
	  dd2=-scalek(iwf)*rijs2*(three+two*rijs+rijs2)*exprij
	 elseif(isc.eq.4) then
	  denij=one/(one+scalek(iwf)*rij)
	  u(1)=rij*denij
	  dd1=denij*denij
	  dd2=-two*scalek(iwf)*denij*dd1
 	 elseif(isc.eq.5) then
	  denij=one/(one+(scalek(iwf)*rij)**3)**third
	  u(1)=rij*denij
	  dd1=denij**4
	  dd2=-four*(scalek(iwf)*rij)**2*scalek(iwf)*dd1*dd1/denij
	endif
       else
        u(1)=rij
        dd1=one
        dd2=zero
      endif

      top=sspinn*b(1,isb,iwf)*u(1)
      topu=sspinn*b(1,isb,iwf)
      topuu=zero

      bot=b(2,isb,iwf)*u(1)
      botu=b(2,isb,iwf)
      botuu=zero

      bot=one+bot
      bot2=bot*bot
      fee=top/bot
      feeu=topu/bot-botu*top/bot2
      feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
      feeuu=feeuu/bot

      feeuu=feeuu*dd1*dd1+feeu*dd2
      feeu=feeu*dd1/rij

      fso(i,j)=fso(i,j) + fee

      v(1,i)=v(1,i) + feeu*rvec_ee(1,ij)
      v(2,i)=v(2,i) + feeu*rvec_ee(2,ij)
      v(3,i)=v(3,i) + feeu*rvec_ee(3,ij)
      v(1,j)=v(1,j) - feeu*rvec_ee(1,ij)
      v(2,j)=v(2,j) - feeu*rvec_ee(2,ij)
      v(3,j)=v(3,j) - feeu*rvec_ee(3,ij)

      d2=d2 + two*(feeuu+ndim1*feeu)

c derivatives of wave function wrt b(2)

      iparm0=nparma(1)
      if(isb.eq.2) iparm0=iparm0+nparmb(1)
      do 30 jparm=1,nparmb(isb)
        iparm=iparm0+jparm

        if(iwjasb(jparm,isb).eq.1) then
          top=u(1)
          topu=one
          topuu=zero

          bot=one+b(2,isb,iwf)*u(1)
          botu=b(2,isb,iwf)
          botuu=zero
         elseif(iwjasb(jparm,isb).eq.2) then
          top=-b(1,isb,iwf)*u(1)*u(1)
          topu=-2*b(1,isb,iwf)*u(1)
          topuu=-2*b(1,isb,iwf)

          bot0=one+b(2,isb,iwf)*u(1)
          bot=bot0*bot0
          botu=2*bot0*b(2,isb,iwf)
          botuu=2*b(2,isb,iwf)*b(2,isb,iwf)
        endif

        bot2=bot*bot
        gee=sspinn*top/bot
        geeu=topu/bot-botu*top/bot2
        geeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
        geeuu=geeuu/bot

        geeuu=geeuu*dd1*dd1+geeu*dd2
        geeu=geeu*dd1/rij

        go(i,j,iparm)=go(i,j,iparm)+gee
        gvalue(iparm)=gvalue(iparm)+gee

        g(1,i,iparm)=g(1,i,iparm) + sspinn*geeu*rvec_ee(1,ij)
        g(2,i,iparm)=g(2,i,iparm) + sspinn*geeu*rvec_ee(2,ij)
        g(3,i,iparm)=g(3,i,iparm) + sspinn*geeu*rvec_ee(3,ij)
        g(1,j,iparm)=g(1,j,iparm) - sspinn*geeu*rvec_ee(1,ij)
        g(2,j,iparm)=g(2,j,iparm) - sspinn*geeu*rvec_ee(2,ij)
        g(3,j,iparm)=g(3,j,iparm) - sspinn*geeu*rvec_ee(3,ij)

        d2g(iparm)=d2g(iparm) + two*sspinn*(geeuu+ndim1*geeu)
  30  continue

      if(nord.eq.0) goto 41

      do 40 ic=1,ncentt
        it=iwctype(ic)

        ri=r_en(i,ic)
        rj=r_en(j,ic)

        if(scalek(iwf).ne.zero) then
          scale2=half*scalek(iwf)
	  if(isc.eq.2) then
            rri=(one-dexp(-scalek(iwf)*ri))/scalek(iwf)
            rrj=(one-dexp(-scalek(iwf)*rj))/scalek(iwf)
            dd3=one-scale2*(rri+rrj)
            dd4=-scale2*(rri-rrj)
            dd5=-scale2*dd3
            dd6=-scale2*dd4
	   elseif(isc.eq.3) then
            ris=scalek(iwf)*ri
            ris2=ris*ris
            ris3=ris*ris2
            rjs=scalek(iwf)*rj
            rjs2=rjs*rjs
            rjs3=rjs*rjs2
            expri=dexp(-ris-half*ris2-third*ris3)
            exprj=dexp(-rjs-half*rjs2-third*rjs3)
            rri=(one-expri)/scalek(iwf)
            rrj=(one-exprj)/scalek(iwf)
            driri=(one+ris+ris2)*expri
            drjrj=(one+rjs+rjs2)*exprj
            d2riri=-scalek(iwf)*ris2*(three+two*ris+ris2)*expri
            d2rjrj=-scalek(iwf)*rjs2*(three+two*rjs+rjs2)*exprj
            dd3=half*(driri+drjrj)
            dd4=half*(driri-drjrj)
            dd5=(d2riri+d2rjrj)/4
            dd6=(d2riri-d2rjrj)/4
           elseif(isc.eq.4) then
            deni=one/(one+scalek(iwf)*ri)
            denj=one/(one+scalek(iwf)*rj)
            rri=ri*deni
            rrj=rj*denj
            dd3=half*(deni*deni+denj*denj)
            dd4=half*(deni*deni-denj*denj)
            dd5=-scale2*(denj**3+deni**3)
            dd6=-scale2*(deni**3-denj**3)
           elseif(isc.eq.5) then
            deni=one/(one+(scalek(iwf)*ri)**3)**third
            denj=one/(one+(scalek(iwf)*rj)**3)**third
            rri=ri*deni
            rrj=rj*denj
            driri=deni**4
            drjrj=denj**4
            d2riri=-four*(scalek(iwf)*ri)**2*scalek(iwf)*deni**7
            d2rjrj=-four*(scalek(iwf)*rj)**2*scalek(iwf)*denj**7
            dd3=half*(driri+drjrj)
            dd4=half*(driri-drjrj)
            dd5=(d2riri+d2rjrj)/4
            dd6=(d2riri-d2rjrj)/4
          endif
         else
          rri=ri
          rrj=rj
          dd3=one
          dd4=zero
          dd5=zero
          dd6=zero
        endif

        ss=r_en(i,ic)+r_en(j,ic)
        tt=r_en(i,ic)-r_en(j,ic)
        u2mt2=rij*rij-tt*tt
        s2mu2=ss*ss-rij*rij
        s2mt2=ss*ss-tt*tt

        s(1)=rri+rrj
        t(1)=rri-rrj
        do 32 jp=2,nord
          u(jp)=u(jp-1)*u(1)
          s(jp)=s(jp-1)*s(1)
   32     t(jp)=t(jp-1)*t(1)

    	fc=zero
  	fu=zero
  	fs=zero
  	ft=zero
  	fuu=zero
  	fss=zero
  	ftt=zero
  	fus=zero
	fut=zero
      	fst=zero

        jparm=1
	ll=0
	jj=1
        do 36 jp=1,nord
	  do 36 ju=jp,0,-1
	    pc=u(ju)
	    pcu=ju*u(ju-1)
	    pcuu=ju*(ju-1)*u(ju-2)

            fp=zero
            fps=zero
            fpss=zero
            fpt=zero
            fptt=zero
            fpst=zero

	    jsx=jp-ju
	    do 35 js=jsx,0,-1
	      ll=ll+1
	      jt=jsx-js

              if(mod(jt,2).ne.0.and.nup.eq.ndn) then
                c(ll,it,iwf)=zero
               else
                p=s(js)*t(jt)
                ps=js*s(js-1)*t(jt)
                pss=js*(js-1)*s(js-2)*t(jt)
                pt=s(js)*jt*t(jt-1)
                ptt=s(js)*jt*(jt-1)*t(jt-2)
                pst=js*s(js-1)*jt*t(jt-1)

    	        fp=fp+c(ll,it,iwf)*p
      	        fps=fps+c(ll,it,iwf)*ps
    	        fpss=fpss+c(ll,it,iwf)*pss
  	        fpt=fpt+c(ll,it,iwf)*pt
    	        fptt=fptt+c(ll,it,iwf)*ptt
     	        fpst=fpst+c(ll,it,iwf)*pst

c derivatives of wave function wrt c-parameters
                ideriv=0
                if(ll.eq.iwc(jj)) then
                  if(nvdepend(jj,it).gt.0) then
                    ideriv=1
                   else
                    jj=jj+1
                  endif
                 elseif(ll.eq.iwjasc(jparm,it)) then
                  ideriv=2
                endif

                if(ideriv.gt.0) then

      	          gu=p*pcu
                  gs=ps*pc
                  gt=pt*pc
    	          guu=p*pcuu
                  gsstt=(pss+ptt)*pc
                  gst=pst*pc
        	  gus=ps*pcu
            	  gut=pt*pcu

                  guu=guu*dd1*dd1+gu*dd2
                  gsstt=gsstt*(dd3**2+dd4**2)+4*gst*dd3*dd4
     &                 +2*(gs*dd5+gt*dd6)
                  tus=gus
                  gus=(tus*dd3+gut*dd4)*dd1
                  gut=(tus*dd4+gut*dd3)*dd1
                  gu=gu*dd1
                  ts=gs
                  gs=ts*dd3+gt*dd4
                  gt=ts*dd4+gt*dd3

                  gu=gu/rij
                  gi=(gs+gt)/r_en(i,ic)
                  gj=(gs-gt)/r_en(j,ic)

                  gp=p*pc
                  if(ideriv.eq.1) then

                  do 33 id=1,nvdepend(jj,it)

                  iparm=npoint(it)+iwdepend(jj,id,it)
                  cd=cdep(jj,id,it)

                  go(i,j,iparm)=go(i,j,iparm)+cd*gp
                  gvalue(iparm)=gvalue(iparm)+cd*gp

                  g(1,i,iparm)=g(1,i,iparm)+cd*(gi*rvec_en(1,i,ic)+gu*rvec_ee(1,ij))
                  g(2,i,iparm)=g(2,i,iparm)+cd*(gi*rvec_en(2,i,ic)+gu*rvec_ee(2,ij))
                  g(3,i,iparm)=g(3,i,iparm)+cd*(gi*rvec_en(3,i,ic)+gu*rvec_ee(3,ij))
                  g(1,j,iparm)=g(1,j,iparm)+cd*(gj*rvec_en(1,j,ic)-gu*rvec_ee(1,ij))
                  g(2,j,iparm)=g(2,j,iparm)+cd*(gj*rvec_en(2,j,ic)-gu*rvec_ee(2,ij))
                  g(3,j,iparm)=g(3,j,iparm)+cd*(gj*rvec_en(3,j,ic)-gu*rvec_ee(3,ij))

 33               d2g(iparm)=d2g(iparm) + cd*two*(guu + gsstt  +
     &            (ndim1*gu + two*(ss*u2mt2*gus+tt*s2mu2*gut) /
     &            (rij*s2mt2) + two*ndim1*(ss*gs-tt*gt)/s2mt2))

                  jj=jj+1

                  elseif(ideriv.eq.2) then

                  iparm=npoint(it)+jparm

                  go(i,j,iparm)=go(i,j,iparm)+gp
                  gvalue(iparm)=gvalue(iparm)+gp

                  g(1,i,iparm)=g(1,i,iparm)+gi*rvec_en(1,i,ic)+gu*rvec_ee(1,ij)
                  g(2,i,iparm)=g(2,i,iparm)+gi*rvec_en(2,i,ic)+gu*rvec_ee(2,ij)
                  g(3,i,iparm)=g(3,i,iparm)+gi*rvec_en(3,i,ic)+gu*rvec_ee(3,ij)
                  g(1,j,iparm)=g(1,j,iparm)+gj*rvec_en(1,j,ic)-gu*rvec_ee(1,ij)
                  g(2,j,iparm)=g(2,j,iparm)+gj*rvec_en(2,j,ic)-gu*rvec_ee(2,ij)
                  g(3,j,iparm)=g(3,j,iparm)+gj*rvec_en(3,j,ic)-gu*rvec_ee(3,ij)

                  d2g(iparm)=d2g(iparm) + two*(guu + gsstt  +
     &            (ndim1*gu + two*(ss*u2mt2*gus+tt*s2mu2*gut) /
     &            (rij*s2mt2) + two*ndim1*(ss*gs-tt*gt)/s2mt2))

                  jparm=jparm+1
                  endif

                endif
              endif
   35       continue

            fc=fc+pc*fp
    	    fu=fu+pcu*fp
    	    fuu=fuu+pcuu*fp
    	    fs=fs+pc*fps
    	    fss=fss+pc*fpss
    	    ft=ft+pc*fpt
    	    ftt=ftt+pc*fptt
    	    fus=fus+pcu*fps
	    fut=fut+pcu*fpt
   36	    fst=fst+pc*fpst

        if(ifock.gt.0) call fock(u(1),s(1),t(1),rri,rrj,it)

        fsstt=fss+ftt

        fuu=fuu*dd1*dd1+fu*dd2
        fsstt=fsstt*(dd3**2+dd4**2)+four*fst*dd3*dd4+two*(fs*dd5+ft*dd6)
        tus=fus
        fus=(tus*dd3+fut*dd4)*dd1
        fut=(tus*dd4+fut*dd3)*dd1
        fu=fu*dd1
        ts=fs
        fs=ts*dd3+ft*dd4
        ft=ts*dd4+ft*dd3

        fu=fu/rij
        fi=(fs+ft)/r_en(i,ic)
        fj=(fs-ft)/r_en(j,ic)

        fso(i,j)=fso(i,j) + fc

        v(1,i)=v(1,i) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        v(2,i)=v(2,i) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        v(3,i)=v(3,i) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        v(1,j)=v(1,j) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        v(2,j)=v(2,j) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        v(3,j)=v(3,j) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)

        d2=d2 + two*(fuu + fsstt  + (ndim1*fu +
     &  two*(ss*u2mt2*fus+tt*s2mu2*fut)/(rij*s2mt2) +
     &  two*ndim1*(ss*fs-tt*ft)/s2mt2))

  40  continue

  41  fsum=fsum+fso(i,j)

  45  continue

c e-n terms
  47  do 55 i=1,nelec

        fso(i,i)=zero
        do 48 iparm=1,nparmj
  48      go(i,i,iparm)=zero

        do 50 ic=1,ncentt

          ri=r_en(i,ic)

          if(scalek(iwf).ne.zero) then
            if(isc.eq.2) then
              rri=(one-dexp(-scalek(iwf)*ri))/scalek(iwf)
              dd7=one-scalek(iwf)*rri
              dd9=-scalek(iwf)*dd7
             elseif(isc.eq.3) then
              ris=scalek(iwf)*ri
              ris2=ris*ris
              ris3=ris*ris2
              expri=dexp(-ris-half*ris2-third*ris3)
              rri=(one-expri)/scalek(iwf)
              dd7=(one+ris+ris2)*expri
              dd9=-scalek(iwf)*ris2*(three+two*ris+ris2)*expri
             elseif(isc.eq.4) then
              deni=one/(one+scalek(iwf)*ri)
              rri=ri*deni
              dd7=deni*deni
              dd9=-two*scalek(iwf)*deni*dd7
             elseif(isc.eq.5) then
              deni=one/(one+(scalek(iwf)*ri)**3)**third
              rri=ri*deni
              dd7=deni**4
              dd9=-four*(scalek(iwf)*ri)**2*scalek(iwf)*dd7*dd7/deni
            endif
           else
            rri=ri
            dd7=one
            dd9=zero
          endif

          rri2=rri*rri
c         rri3=rri2*rri
c         rri4=rri3*rri

          top=a(1,iwf)*rri
          topi=a(1,iwf)
          topii=zero

          bot=a(2,iwf)*rri
          boti=a(2,iwf)
          botii=zero

          bot=one+bot
          bot2=bot*bot
          fen=top/bot
          feni=topi/bot-boti*top/bot2
          fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
          fenii=fenii/bot

          fso(i,i)=fso(i,i)+fen
          fenii=fenii*dd7*dd7+feni*dd9
          feni=feni*dd7/r_en(i,ic)

          if(ipr.ge.3) write(6,'(''d:a(1,iwf),a(2,iwf),rri,top,bot,fsum'',9f10.4)')
     &    a(1,iwf),a(2,iwf),rri,top,bot,fsum

          v(1,i)=v(1,i) + feni*rvec_en(1,i,ic)
          v(2,i)=v(2,i) + feni*rvec_en(2,i,ic)
          v(3,i)=v(3,i) + feni*rvec_en(3,i,ic)

          d2 = d2 + fenii + ndim1*feni

          do 49 jparm=1,nparma(1)

            if(iwjasa(jparm,1).eq.1) then
              top=rri
              topi=one
              topii=zero

              bot=one+a(2,iwf)*rri
              boti=a(2,iwf)
              botii=zero
             elseif(iwjasa(jparm,1).eq.2) then
              top=-a(1,iwf)*rri2
              topi=-2*a(1,iwf)*rri
              topii=-2*a(1,iwf)

              bot0=one+a(2,iwf)*rri
              bot=bot0*bot0
              boti=2*bot0*a(2,iwf)
              botii=2*a(2,iwf)*a(2,iwf)
            endif

            bot2=bot*bot
            geni=topi/bot-boti*top/bot2
            genii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
            genii=genii/bot

            genii=genii*dd7*dd7+geni*dd9
            geni=geni*dd7/r_en(i,ic)
            gen=top/bot

            go(i,i,jparm)=go(i,i,jparm)+gen
            gvalue(jparm)=gvalue(jparm)+gen

            g(1,i,jparm)=g(1,i,jparm)+geni*rvec_en(1,i,ic)
            g(2,i,jparm)=g(2,i,jparm)+geni*rvec_en(2,i,ic)
            g(3,i,jparm)=g(3,i,jparm)+geni*rvec_en(3,i,ic)

            d2g(jparm)=d2g(jparm)+genii+two*geni
  49      continue

  50    continue

        fsum=fsum+fso(i,i)

   55 continue

   60 continue

      fsumo=fsum

      value=fsum

      return
      end

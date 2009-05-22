      subroutine jastrow3(x,v,d2,div_vj,value)
c Written by Claudia Filippi and Cyrus Umrigar

      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use wfsec_mod
      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'
!JT      include 'force.h'

!JT      parameter (one=1.d0,two=2.d0,three=3.d0,four=4.d0
!JT     &,half=.5d0,third=1.d0/3.d0)

!JT      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
!JT     &,ifock,i3body,irewgt,iaver,istrch
!JT     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent

!JT      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      common /focktmp/ fc,fcu,fcuu,fcs,fcss,fct,fctt,fcst,fcus,fcut
      common /bparm/ nspin2b,nocuspb
      common /jaso/ fso(MELEC,MELEC),fijo(3,MELEC,MELEC)
     &,d2ijo(MELEC,MELEC),d2o,fsumo,fjo(3,MELEC)

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      dimension x(3,*),v(3,*),div_vj(*)
      dimension u(-2:MORDJ),s(-2:MORDJ),t(-2:MORDJ)

      fsum=0
      u(-2)=0
      s(-2)=0
      t(-2)=0
      u(-1)=0
      s(-1)=0
      t(-1)=0
      u(0)=one
      s(0)=one
      t(0)=one

      if(nelec.lt.2) goto 47

c e-e and e-e-n terms
      ij=0
      do 45 i=2,nelec
      im1=i-1
      do 45 j=1,im1
      ij=ij+1

      sspinn=one
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

      fijo(1,i,j)=0
      fijo(2,i,j)=0
      fijo(3,i,j)=0
      fijo(1,j,i)=0
      fijo(2,j,i)=0
      fijo(3,j,i)=0
      d2ijo(i,j)=0
      fso(i,j)=0

      rij=r_ee(ij)

c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 are supposed to have 0 2nd and 3rd derivatives at 0
c however, isc = 3 has an error.  It should be
c exprij=dexp(-rijs-half*rijs2-(5./6.)*rijs3)
c with other corresponding changes.
      if(scalek(iwf).ne.0) then
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
        dd2=0
      endif

      top=b(1,isb,iwf)*u(1)
      topu=b(1,isb,iwf)
      topuu=0

      bot=b(2,isb,iwf)*u(1)
      botu=b(2,isb,iwf)
      botuu=0

      bot=one+bot
      bot2=bot*bot
      fee=top/bot
      feeu=topu/bot-botu*top/bot2
      feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
      feeuu=feeuu/bot

      feeuu=feeuu*dd1*dd1+feeu*dd2
      feeu=feeu*dd1/rij

      fso(i,j)=fso(i,j) + sspinn*fee

      fijo(1,i,j)=fijo(1,i,j) + sspinn*feeu*rvec_ee(1,ij)
      fijo(2,i,j)=fijo(2,i,j) + sspinn*feeu*rvec_ee(2,ij)
      fijo(3,i,j)=fijo(3,i,j) + sspinn*feeu*rvec_ee(3,ij)
      fijo(1,j,i)=fijo(1,j,i) - sspinn*feeu*rvec_ee(1,ij)
      fijo(2,j,i)=fijo(2,j,i) - sspinn*feeu*rvec_ee(2,ij)
      fijo(3,j,i)=fijo(3,j,i) - sspinn*feeu*rvec_ee(3,ij)

      d2ijo(i,j)=d2ijo(i,j) + two*sspinn*(feeuu+(ndim-1)*feeu)

      if(nord.eq.0) goto 41

      do 40 ic=1,ncent
        it=iwctype(ic)

        ri=r_en(i,ic)
        rj=r_en(j,ic)

        if(scalek(iwf).ne.0) then
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
          dd4=0
          dd5=0
          dd6=0
        endif

        ss=r_en(i,ic)+r_en(j,ic)
        tt=r_en(i,ic)-r_en(j,ic)
        u2mt2=rij*rij-tt*tt
        s2mu2=ss*ss-rij*rij
        s2mt2=ss*ss-tt*tt

        s(1)=rri+rrj
        t(1)=rri-rrj
        do 31 jp=2,nord
          u(jp)=u(jp-1)*u(1)
          s(jp)=s(jp-1)*s(1)
   31     t(jp)=t(jp-1)*t(1)

    	fc=0
  	fcu=0
  	fcs=0
  	fct=0
  	fcuu=0
  	fcss=0
  	fctt=0
  	fcus=0
	fcut=0
      	fcst=0

	ll=0
        do 36 jp=1,nord
	  do 36 ju=jp,0,-1
	    pc=u(ju)
	    pcu=ju*u(ju-1)
	    pcuu=ju*(ju-1)*u(ju-2)

            fp=0
            fps=0
            fpss=0
            fpt=0
            fptt=0
            fpst=0

	    jsx=jp-ju
	    do 35 js=jsx,0,-1
	      ll=ll+1
	      jt=jsx-js

              if(mod(jt,2).ne.0.and.nup.eq.ndn) then
                c(ll,it,iwf)=0
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
              endif
   35       continue

            fc=fc+pc*fp
    	    fcu=fcu+pcu*fp
    	    fcuu=fcuu+pcuu*fp
    	    fcs=fcs+pc*fps
    	    fcss=fcss+pc*fpss
    	    fct=fct+pc*fpt
    	    fctt=fctt+pc*fptt
    	    fcus=fcus+pcu*fps
	    fcut=fcut+pcu*fpt
   36	    fcst=fcst+pc*fpst

c     write(6,'(''u,s,t,fc,fcu,fcs,fct'',9f9.5)') u(1),s(1),t(1),fc,fcu,fcs,fct,fcss,fctt,fcus,fcut,fcst
        if(ifock.gt.0) call fock(u(1),s(1),t(1),rri,rrj,it)

        top=fc

        fu=fcu
        fs=fcs
        ft=fct

        fuu=fcuu
        fsstt=fcss+fctt
        fus=fcus
        fut=fcut
        fst=fcst

        fuu=fuu*dd1*dd1+fu*dd2
        fsstt=fsstt*(dd3**2+dd4**2)+four*fst*dd3*dd4+two*(fs*dd5+ft*dd6)
        gus=fus
        fus=(fus*dd3+fut*dd4)*dd1
        fut=(gus*dd4+fut*dd3)*dd1
        fu=fu*dd1
        gs=fs
        fs=fs*dd3+ft*dd4
        ft=gs*dd4+ft*dd3

        fu=fu/rij
        fi=(fs+ft)/r_en(i,ic)
        fj=(fs-ft)/r_en(j,ic)

        fso(i,j)=fso(i,j) + top

        fijo(1,i,j)=fijo(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        fijo(2,i,j)=fijo(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        fijo(3,i,j)=fijo(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        fijo(1,j,i)=fijo(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        fijo(2,j,i)=fijo(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        fijo(3,j,i)=fijo(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)

c Warning: The generalization to arbitrary dimension is to be checked but is most likely OK.
c We never use this for anything except 3d anyway.
c       d2ijo(i,j)=d2ijo(i,j) + two*(fuu + fsstt  + two*(fu +
c    &  (ss*u2mt2*fus+tt*s2mu2*fut)/(rij*s2mt2) +
c    &  two*(ss*fs-tt*ft)/s2mt2))
        d2ijo(i,j)=d2ijo(i,j) + 2*((ndim-1)*(fu + 2*(ss*fs-tt*ft)/s2mt2)
     &  + fuu + fsstt + 2*(ss*u2mt2*fus+tt*s2mu2*fut)/(rij*s2mt2))

  40  continue

  41  v(1,i)=v(1,i)+fijo(1,i,j)
      v(2,i)=v(2,i)+fijo(2,i,j)
      v(3,i)=v(3,i)+fijo(3,i,j)
      v(1,j)=v(1,j)+fijo(1,j,i)
      v(2,j)=v(2,j)+fijo(2,j,i)
      v(3,j)=v(3,j)+fijo(3,j,i)
      div_vj(i)=div_vj(i)+d2ijo(i,j)/2
      div_vj(j)=div_vj(j)+d2ijo(i,j)/2
      d2=d2+d2ijo(i,j)
      fsum=fsum+fso(i,j)

  45  continue

c     write(6,'(''fu,fuu,fus,fut,fsstt,d2,f(1,1)='',9f9.5)')
c    &fu,fuu,fus,fut,fsstt,d2,v(1,1)

c e-n terms
  47  do 55 i=1,nelec

        fijo(1,i,i)=0
        fijo(2,i,i)=0
        fijo(3,i,i)=0
        d2ijo(i,i)=0
        fso(i,i)=0

        do 50 ic=1,ncent

          ri=r_en(i,ic)

          if(scalek(iwf).ne.0) then
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
            dd9=0
          endif

          top=a(1,iwf)*rri
          topi=a(1,iwf)
          topii=0

          bot=a(2,iwf)*rri
          boti=a(2,iwf)
          botii=0

          bot=one+bot
          bot2=bot*bot
          fen=top/bot
          feni=topi/bot-boti*top/bot2
          fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
          fenii=fenii/bot

          fso(i,i)=fso(i,i)+fen
          fenii=fenii*dd7*dd7+feni*dd9
          feni=feni*dd7/r_en(i,ic)

          fijo(1,i,i)=fijo(1,i,i) + feni*rvec_en(1,i,ic)
          fijo(2,i,i)=fijo(2,i,i) + feni*rvec_en(2,i,ic)
          fijo(3,i,i)=fijo(3,i,i) + feni*rvec_en(3,i,ic)

   50     d2ijo(i,i) = d2ijo(i,i) + fenii + (ndim-1)*feni

      fsum=fsum+fso(i,i)
      v(1,i)=v(1,i)+fijo(1,i,i)
      v(2,i)=v(2,i)+fijo(2,i,i)
      v(3,i)=v(3,i)+fijo(3,i,i)
      div_vj(i)=div_vj(i)+d2ijo(i,i)
      d2=d2+d2ijo(i,i)

   55 continue

   60 continue

      fsumo=fsum
      d2o=d2
      do 70 i=1,nelec
        fjo(1,i)=v(1,i)
        fjo(2,i)=v(2,i)
   70   fjo(3,i)=v(3,i)

      value=fsum

      return
      end

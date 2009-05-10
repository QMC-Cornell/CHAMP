      subroutine jastrow3e(iel,x,v,value)
c Written by Claudia Filippi by modifying jastrow3

      use atom_mod
      use dets_mod
      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'
!JT      include 'force.h'
!JT      include 'pseudo.h'

!JT      parameter (one=1.d0)
!JT      parameter (half=.5d0,third=1.d0/3.d0)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /focktmp/ fc,fcu,fcuu,fcs,fcss,fct,fctt,fcst,fcus,fcut
      common /bparm/ nspin2b,nocuspb

      common /jasn/ fsn(MELEC,MELEC),fijn(3,MELEC,MELEC)
     &,d2ijn(MELEC,MELEC),d2n,fsumn,fjn(3,MELEC)

      common /jaso/ fso(MELEC,MELEC),fijo(3,MELEC,MELEC)
     &,d2ijo(MELEC,MELEC),d2o,fsumo,fjo(3,MELEC)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      dimension x(3,*),v(3,*)
      dimension u(-2:MORDJ),s(-2:MORDJ),t(-2:MORDJ)

      u(-2)=0
      s(-2)=0
      t(-2)=0
      u(-1)=0
      s(-1)=0
      t(-1)=0
      u(0)=1
      s(0)=1
      t(0)=1

      do 1 i=1,nelec
        fjn(1,i)=fjo(1,i)
        fjn(2,i)=fjo(2,i)
    1   fjn(3,i)=fjo(3,i)
      fsumn=fsumo

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
      ij=((i-1)*(i-2))/2+j

      fijn(1,i,j)=0
      fijn(2,i,j)=0
      fijn(3,i,j)=0
      fijn(1,j,i)=0
      fijn(2,j,i)=0
      fijn(3,j,i)=0
      fsn(i,j)=0

      sspinn=1
      if(i.le.nup .or. j.gt.nup) then
        isb=1
        if(nspin2b.eq.2) then
          isb=2
         elseif(nocuspb.eq.0) then
          if(ndim.eq.3) then
            sspinn=half
           elseif(ndim.eq.2) then
            sspinn=third
          endif
        endif
       else
        isb=1
      endif

      rij=r_ee(ij)
c     riji=one/rij

      if(scalek(iwf).ne.0) then
        if(isc.eq.2) then
          u(1)=(one-dexp(-scalek(iwf)*rij))/scalek(iwf)
          dd1=one-scalek(iwf)*u(1)
	 elseif(isc.eq.3) then
	  rijs=scalek(iwf)*rij
          rijs2=rijs*rijs
          rijs3=rijs*rijs2
	  exprij=dexp(-rijs-half*rijs2-third*rijs3)
	  u(1)=(one-exprij)/scalek(iwf)
	  dd1=(one+rijs+rijs2)*exprij
	 elseif(isc.eq.4) then
	  denij=one/(one+scalek(iwf)*rij)
	  u(1)=rij*denij
	  dd1=denij*denij
 	 elseif(isc.eq.5) then
	  denij=one/(one+(scalek(iwf)*rij)**3)**third
	  u(1)=rij*denij
	  dd1=denij**4
	endif
       else
        u(1)=rij
        dd1=one
      endif

      top=b(1,isb,iwf)*u(1)
      topu=b(1,isb,iwf)
c     topuu=0

      bot=b(2,isb,iwf)*u(1)
      botu=b(2,isb,iwf)
c     botuu=0

      bot=one+bot
      bot2=bot*bot
      fee=top/bot
      feeu=topu/bot-botu*top/bot2
      feeu=feeu*dd1/rij

      fsn(i,j)=fsn(i,j) + sspinn*fee

      fijn(1,i,j)=fijn(1,i,j) + sspinn*feeu*rvec_ee(1,ij)
      fijn(2,i,j)=fijn(2,i,j) + sspinn*feeu*rvec_ee(2,ij)
      fijn(3,i,j)=fijn(3,i,j) + sspinn*feeu*rvec_ee(3,ij)
      fijn(1,j,i)=fijn(1,j,i) - sspinn*feeu*rvec_ee(1,ij)
      fijn(2,j,i)=fijn(2,j,i) - sspinn*feeu*rvec_ee(2,ij)
      fijn(3,j,i)=fijn(3,j,i) - sspinn*feeu*rvec_ee(3,ij)

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
c           d2riri=-scalek(iwf)*ris2*(3+2*ris+ris2)*expri
c           d2rjrj=-scalek(iwf)*rjs2*(3+2*rjs+rjs2)*exprj
            dd3=half*(driri+drjrj)
            dd4=half*(driri-drjrj)
           elseif(isc.eq.4) then
            deni=one/(one+scalek(iwf)*ri)
            denj=one/(one+scalek(iwf)*rj)
            rri=ri*deni
            rrj=rj*denj
            dd3=half*(deni*deni+denj*denj)
            dd4=half*(deni*deni-denj*denj)
           elseif(isc.eq.5) then
            deni=one/(one+(scalek(iwf)*ri)**3)**third
            denj=one/(one+(scalek(iwf)*rj)**3)**third
            rri=ri*deni
            rrj=rj*denj
            driri=deni**4
            drjrj=denj**4
c           d2riri=-four*(scalek(iwf)*ri)**2*scalek(iwf)*deni**7
c           d2rjrj=-four*(scalek(iwf)*rj)**2*scalek(iwf)*denj**7
            dd3=half*(driri+drjrj)
            dd4=half*(driri-drjrj)
          endif
         else
          rri=ri
          rrj=rj
          dd3=one
          dd4=0
        endif

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

	ll=0
        do 36 jp=1,nord
	  do 36 ju=jp,0,-1
	    pc=u(ju)
	    pcu=ju*u(ju-1)

            fp=0
            fps=0
            fpt=0

	    jsx=jp-ju
	    do 34 js=jsx,0,-1
	      ll=ll+1
	      jt=jsx-js

              p=s(js)*t(jt)
              ps=js*s(js-1)*t(jt)
              pt=s(js)*jt*t(jt-1)

    	      fp=fp+c(ll,it,iwf)*p
  	      fps=fps+c(ll,it,iwf)*ps
   34         fpt=fpt+c(ll,it,iwf)*pt

              fc=fc+pc*fp
    	      fcu=fcu+pcu*fp
    	      fcs=fcs+pc*fps
   36         fct=fct+pc*fpt

        if(ifock.gt.0) call fock(u(1),s(1),t(1),rri,rrj,it)

        fu=fcu
        fs=fcs
        ft=fct

        fu=fu*dd1
        gs=fs
        fs=fs*dd3+ft*dd4
        ft=gs*dd4+ft*dd3

        fu=fu/rij
        fi=(fs+ft)/r_en(i,ic)
        fj=(fs-ft)/r_en(j,ic)

        fsn(i,j)=fsn(i,j) + fc

        fijn(1,i,j)=fijn(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        fijn(2,i,j)=fijn(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        fijn(3,i,j)=fijn(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        fijn(1,j,i)=fijn(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        fijn(2,j,i)=fijn(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        fijn(3,j,i)=fijn(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)

  40  continue

  41  fsumn=fsumn+fsn(i,j)-fso(i,j)
      fjn(1,i)=fjn(1,i)+fijn(1,i,j)-fijo(1,i,j)
      fjn(2,i)=fjn(2,i)+fijn(2,i,j)-fijo(2,i,j)
      fjn(3,i)=fjn(3,i)+fijn(3,i,j)-fijo(3,i,j)
      fjn(1,j)=fjn(1,j)+fijn(1,j,i)-fijo(1,j,i)
      fjn(2,j)=fjn(2,j)+fijn(2,j,i)-fijo(2,j,i)
      fjn(3,j)=fjn(3,j)+fijn(3,j,i)-fijo(3,j,i)

  45  continue

  47  fijn(1,iel,iel)=0
      fijn(2,iel,iel)=0
      fijn(3,iel,iel)=0
      fsn(iel,iel)=0
      d2ijn(iel,iel)=0

      do 50 ic=1,ncent
        it=iwctype(ic)

        ri=r_en(iel,ic)

        if(scalek(iwf).ne.0) then
          if(isc.eq.2) then
            rri=(one-dexp(-scalek(iwf)*ri))/scalek(iwf)
            dd7=one-scalek(iwf)*rri
           elseif(isc.eq.3) then
            ris=scalek(iwf)*ri
            ris2=ris*ris
            ris3=ris*ris2
            expri=dexp(-ris-half*ris2-third*ris3)
            rri=(one-expri)/scalek(iwf)
            dd7=(one+ris+ris2)*expri
           elseif(isc.eq.4) then
            deni=one/(one+scalek(iwf)*ri)
            rri=ri*deni
            dd7=deni*deni
           elseif(isc.eq.5) then
            deni=one/(one+(scalek(iwf)*ri)**3)**third
            rri=ri*deni
            dd7=deni**4
          endif
         else
          rri=ri
          dd7=one
        endif

        top=a(1,iwf)*rri
        topi=a(1,iwf)

        bot=a(2,iwf)*rri
        boti=a(2,iwf)

        bot=one+bot
        bot2=bot*bot
        fen=top/bot
        feni=topi/bot-boti*top/bot2

        fsn(iel,iel)=fsn(iel,iel)+fen
        feni=feni*dd7/r_en(iel,ic)

        fijn(1,iel,iel)=fijn(1,iel,iel) + feni*rvec_en(1,iel,ic)
        fijn(2,iel,iel)=fijn(2,iel,iel) + feni*rvec_en(2,iel,ic)
   50   fijn(3,iel,iel)=fijn(3,iel,iel) + feni*rvec_en(3,iel,ic)

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)
      fjn(1,iel)=fjn(1,iel)+fijn(1,iel,iel)-fijo(1,iel,iel)
      fjn(2,iel)=fjn(2,iel)+fijn(2,iel,iel)-fijo(2,iel,iel)
      fjn(3,iel)=fjn(3,iel)+fijn(3,iel,iel)-fijo(3,iel,iel)

      if(ipr.ge.3) write(6,'(''a(1,iwf),a(2,iwf),rri,top,bot,fsumn'',9f10.4)')
     &a(1,iwf),a(2,iwf),rri,top,bot,fsumn

   60 continue

      do 70 i=1,nelec
        v(1,i)=fjn(1,i)
        v(2,i)=fjn(2,i)
   70   v(3,i)=fjn(3,i)

      value=fsumn

      return
      end

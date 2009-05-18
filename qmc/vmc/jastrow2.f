      subroutine jastrow2(x,v,d2,div_vj,value)
c Written by Cyrus Umrigar

      use atom_mod
      use dets_mod
      use const_mod
      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'
!JT      include 'force.h'

!JT      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0
!JT     &,half=.5d0,third=1.d0/3.d0)

      common /dim/ ndim
      common /pars/ a00,a20,a21,eps_fock,c0000,c1110,c2000,
     &   xm1,xm2,xm12,xms,xma,Zfock

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar2/ a1(MPARMJ,3,MWF),a2(MPARMJ,3,MWF)
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /gamder/ dsphi21,dtphi21,duphi21,d2sphi21,d2tphi21,
     &d2uphi21,d2stphi21,d2suphi21,d2utphi21,
     &                dsphi31,dtphi31,duphi31,d2sphi31,d2tphi31,
     &d2uphi31,d2stphi31,d2suphi31,d2utphi31,
     &                dsphi20,dtphi20,duphi20,d2sphi20,d2tphi20,
     &d2uphi20,d2stphi20,d2suphi20,d2utphi20,
     &                dsy21,dty21,duy21,d2sy21,d2ty21,
     &d2uy21,d2sty21,d2suy21,d2uty21
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent

      common /wfsec/ iwftype(MFORCE),iwf,nwftype


      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
      dimension x(3,*),v(3,*),div_vj(*)

      fsum=zero
      if(nelec.lt.2) goto 60

      ij=0
      do 40 i=2,nelec
      im1=i-1
      do 40 j=1,im1
      ij=ij+1
      if(i.le.nup .or. j.gt.nup) then
        if(nspin2.ge.2) then
          sspinn=one
          is=2
          if(nspin2.eq.3 .and. j.gt.nup) is=3
         else
          is=1
          if(ndim.eq.3) then
            sspinn=half
           elseif(ndim.eq.2) then
            sspinn=third
          endif
        endif
       else
        sspinn=one
        is=1
      endif

      rij=r_ee(ij)

      do 40 ic=1,ncent

        wtj=one/ncent
        sspin=sspinn*wtj

        ri=r_en(i,ic)
        rj=r_en(j,ic)

c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 are supposed to have zero 2nd and 3rd derivatives at 0
c however, isc = 3 has an error.  It should be
c exprij=dexp(-rijs-half*rijs2-(5./6.)*rijs3)
c with other corresponding changes.
        if(scalek(iwf).ne.zero) then
          scale2=half*scalek(iwf)
          if(isc.eq.2) then
            u=(one-dexp(-scalek(iwf)*rij))/scalek(iwf)
            rri=(one-dexp(-scalek(iwf)*ri))/scalek(iwf)
            rrj=(one-dexp(-scalek(iwf)*rj))/scalek(iwf)
            dd1=one-scalek(iwf)*u
            dd2=-scalek(iwf)*dd1
            dd3=one-scale2*(rri+rrj)
            dd4=-scale2*(rri-rrj)
            dd5=-scale2*dd3
            dd6=-scale2*dd4
           elseif(isc.eq.3) then
            rijs=scalek(iwf)*rij
            rijs2=rijs*rijs
            rijs3=rijs*rijs2
            ris=scalek(iwf)*ri
            ris2=ris*ris
            ris3=ris*ris2
            rjs=scalek(iwf)*rj
            rjs2=rjs*rjs
            rjs3=rjs*rjs2
            exprij=dexp(-rijs-half*rijs2-third*rijs3)
            expri=dexp(-ris-half*ris2-third*ris3)
            exprj=dexp(-rjs-half*rjs2-third*rjs3)
            u=(one-exprij)/scalek(iwf)
            rri=(one-expri)/scalek(iwf)
            rrj=(one-exprj)/scalek(iwf)
            dd1=(one+rijs+rijs2)*exprij
            dd2=-scalek(iwf)*rijs2*(three+two*rijs+rijs2)*exprij
            driri=(one+ris+ris2)*expri
            drjrj=(one+rjs+rjs2)*exprj
            d2riri=-scalek(iwf)*ris2*(three+two*ris+ris2)*expri
            d2rjrj=-scalek(iwf)*rjs2*(three+two*rjs+rjs2)*exprj
            dd3=half*(driri+drjrj)
            dd4=half*(driri-drjrj)
            dd5=(d2riri+d2rjrj)/4
            dd6=(d2riri-d2rjrj)/4
           elseif(isc.eq.4) then
            denij=one/(one+scalek(iwf)*rij)
            deni=one/(one+scalek(iwf)*ri)
            denj=one/(one+scalek(iwf)*rj)
            u=rij/(one+scalek(iwf)*rij)
            rri=ri*deni
            rrj=rj*denj
            dd1=denij*denij
            dd2=-two*scalek(iwf)*denij*dd1
            dd3=half*(deni*deni+denj*denj)
            dd4=half*(deni*deni-denj*denj)
            dd5=-scale2*(denj**3+deni**3)
            dd6=-scale2*(deni**3-denj**3)
           elseif(isc.eq.5) then
            denij=one/(one+(scalek(iwf)*rij)**3)**third
            deni=one/(one+(scalek(iwf)*ri)**3)**third
            denj=one/(one+(scalek(iwf)*rj)**3)**third
            u=rij*denij
            rri=ri*deni
            rrj=rj*denj
            dd1=denij**4
            dd2=-four*(scalek(iwf)*rij)**2*scalek(iwf)*dd1*dd1/denij
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
          u=rij
          rri=r_en(i,ic)
          rrj=r_en(j,ic)
          dd1=one
          dd2=zero
          dd3=one
          dd4=zero
          dd5=zero
          dd6=zero
        endif

        s=rri+rrj
        t=(rri-rrj)
        u2=u*u
        u3=u*u2
        u4=u2*u2
        u5=u3*u2
        s2=s*s
        s3=s*s2
        s4=s2*s2
        si=one/s
        si2=si*si
        si3=si2*si
        si4=si2*si2
	si5=si3*si2
        t2=t*t
        t4=t2*t2
        us=u*s
        rr2=(s2+t2)/2
        rr=dsqrt(rr2)
        rr3=rr2*rr
        rr4=rr3*rr
        rr5=rr4*rr
        rrlog=dlog(rr2)
        ss=r_en(i,ic)+r_en(j,ic)
        tt=(r_en(i,ic)-r_en(j,ic))
        u2mt2=rij*rij-tt*tt
        s2mu2=ss*ss-rij*rij
        s2mt2=ss*ss-tt*tt

        term38=(rr2-u2)*rrlog
c       term1=u*(12*rr2+u2)/rr3
        terma1=(a1(35,is,iwf)*u/rr - a1(37,is,iwf)*u3/rr3)/2
     &  +((rr2*rrlog+rr2-u2)/rr2)*(a1(38,is,iwf)+a1(39,is,iwf)*u)
        terma2=(a1(35,is,iwf)/rr-3*a1(37,is,iwf)*u2/rr3)/2
     &  - 2*a1(38,is,iwf)*u/rr2
     &  + a1(39,is,iwf)*(rr2*rrlog+rr2-3*u2)/rr2
        termb1=(a2(35,is,iwf)*u/rr - a2(37,is,iwf)*u3/rr3)/2
     &  +((rr2*rrlog+rr2-u2)/rr2)*(a2(38,is,iwf)+a2(39,is,iwf)*u)
        termb2=(a2(35,is,iwf)/rr-3*a2(37,is,iwf)*u2/rr3)/2
     &  - 2*a2(38,is,iwf)*u/rr2
     &  + a2(39,is,iwf)*(rr2*rrlog+rr2-3*u2)/rr2
        st35=-us/(4*rr3)
        st36=t2/(4*rr3)
        st37=three*u3*s/(four*rr5)
        st38=s*(rr2+u2)/rr4
        st39=u*st38

        termsa31=a1(47,is,iwf)*u3+a1(48,is,iwf)*s3+a1(49,is,iwf)*u2*s
     &  +a1(50,is,iwf)*u*s2+ (a1(51,is,iwf)*u+a1(52,is,iwf)*s)*t2
        du_termsa31=a1(50,is,iwf)*s2+a1(51,is,iwf)*t2
     &  +2*a1(49,is,iwf)*s*u+3*a1(47,is,iwf)*u2
        ds_termsa31=3*a1(48,is,iwf)*s2+a1(52,is,iwf)*t2
     &  +2*a1(50,is,iwf)*s*u+a1(49,is,iwf)*u2
        dt_termsa31=2*a1(52,is,iwf)*s + 2*a1(51,is,iwf)*u
        d2u_termsa31=2*a1(49,is,iwf)*s + 6*a1(47,is,iwf)*u
        d2sd2t_termsa31=2*(3*a1(48,is,iwf)*s+a1(50,is,iwf)*u
     &  +a1(52,is,iwf)*s
     &  +a1(51,is,iwf)*u
     &  )
        d2us_termsa31=2*(a1(50,is,iwf)*s+a1(49,is,iwf)*u)
        d2ut_termsa31=2*a1(51,is,iwf)
        d2st_termsa31=2*a1(52,is,iwf)

        termsb31=a2(47,is,iwf)*u3+a2(48,is,iwf)*s3+a2(49,is,iwf)*u2*s
     &  +a2(50,is,iwf)*u*s2
     &  + (a2(51,is,iwf)*u+a2(52,is,iwf)*s)*t2
        du_termsb31=a2(50,is,iwf)*s2+a2(51,is,iwf)*t2
     &  +2*a2(49,is,iwf)*s*u+3*a2(47,is,iwf)*u2
        ds_termsb31=3*a2(48,is,iwf)*s2+a2(52,is,iwf)*t2
     &  +2*a2(50,is,iwf)*s*u+a2(49,is,iwf)*u2
        dt_termsb31=2*a2(52,is,iwf)*s + 2*a2(51,is,iwf)*u
        d2u_termsb31=2*a2(49,is,iwf)*s + 6*a2(47,is,iwf)*u
        d2sd2t_termsb31=2*(3*a2(48,is,iwf)*s+a2(50,is,iwf)*u
     &  +a2(52,is,iwf)*s
     &  +a2(51,is,iwf)*u
     &  )
        d2us_termsb31=2*(a2(50,is,iwf)*s+a2(49,is,iwf)*u)
        d2ut_termsb31=2*a2(51,is,iwf)
        d2st_termsb31=2*a2(52,is,iwf)

        ds_rrlog=s/rr2
        dt_rrlog=one/rr2
        d2st_rrlog=-s/rr4

        if(ifock.ge.1) call psigm(u,rri,rrj,phi21,phi20,phi31,1)

        y21=rr2-u2

        au=a1(1,is,iwf)+2*a1(4,is,iwf)*u+a1(7,is,iwf)*s
     &  +3*a1(10,is,iwf)*u2+2*a1(13,is,iwf)*u*s
     &  +a1(14,is,iwf)*s2+a1(16,is,iwf)*t2+ 4*a1(20,is,iwf)*u3
     &  +3*a1(23,is,iwf)*u2*s
     &  +a1(24,is,iwf)*s3+a1(31,is,iwf)*s*t2+2*u*(a1(32,is,iwf)*s2
     &  +a1(33,is,iwf)*t2)+a1(35,is,iwf)*rr+3*a1(37,is,iwf)*u2/rr
     &  +(-2*a1(38,is,iwf)*u+a1(39,is,iwf)*(rr2-3*u2))*rrlog
     &  +a1(40,is,iwf)*duphi21+a1(44,is,iwf)*duphi20
     &  +a1(45,is,iwf)*duphi31+a21*duy21
     &  +du_termsa31*rrlog
     &  +(3*a1(53,is,iwf)*u2+a1(54,is,iwf)*t2)*si
     &  +(4*a1(55,is,iwf)*u3+2*a1(56,is,iwf)*u*t2)*si2
     &  +(5*a1(58,is,iwf)*u4+3*a1(59,is,iwf)*u2*t2+a1(60,is,iwf)*t4)*si3
     &  +(4*a1(61,is,iwf)*u3+2*a1(62,is,iwf)*u*t2)*si
     &  +(5*a1(64,is,iwf)*u4+3*a1(65,is,iwf)*u2*t2+a1(66,is,iwf)*t4)*si2
     &  +(5*a1(67,is,iwf)*u4+3*a1(68,is,iwf)*u2*t2+a1(69,is,iwf)*t4)*si

        as=a1(2,is,iwf)+2*a1(5,is,iwf)*s+a1(7,is,iwf)*u
     &  +3*a1(11,is,iwf)*s2+a1(13,is,iwf)*u2+2*a1(14,is,iwf)*u*s
     &  +a1(18,is,iwf)*t2+4*a1(21,is,iwf)*s3+a1(23,is,iwf)*u3
     &  +3*a1(24,is,iwf)*u*s2+a1(31,is,iwf)*u*t2+2*s*(a1(32,is,iwf)*u2
     &  +a1(34,is,iwf)*t2)+a1(36,is,iwf)*(s2+t2/2)/rr +s*terma1
     &  +a1(40,is,iwf)*dsphi21+a1(44,is,iwf)*dsphi20
     &  +a1(45,is,iwf)*dsphi31+a21*dsy21
     &  +ds_termsa31*rrlog + termsa31*ds_rrlog
     &  -(a1(53,is,iwf)*u3+a1(54,is,iwf)*u*t2)*si2
     &  -2*(a1(55,is,iwf)*u4+a1(56,is,iwf)*u2*t2+a1(57,is,iwf)*t4)*si3
     &  -3*(a1(58,is,iwf)*u5+a1(59,is,iwf)*u3*t2+a1(60,is,iwf)*u*t4)*si4
     &  -  (a1(61,is,iwf)*u4+a1(62,is,iwf)*u2*t2+a1(63,is,iwf)*t4)*si2
     &  -2*(a1(64,is,iwf)*u5+a1(65,is,iwf)*u3*t2+a1(66,is,iwf)*u*t4)*si3
     &  -  (a1(67,is,iwf)*u5+a1(68,is,iwf)*u3*t2+a1(69,is,iwf)*u*t4)*si2

        at=t*(2*(a1(6,is,iwf)+a1(16,is,iwf)*u+a1(18,is,iwf)*s
     &  +2*a1(22,is,iwf)*t2
     &  +a1(31,is,iwf)*u*s+a1(33,is,iwf)*u2+a1(34,is,iwf)*s2)
     &  +a1(36,is,iwf)*s/(2*rr) +terma1
     &  +dt_termsa31*rrlog + termsa31*dt_rrlog
     &  +2*a1(54,is,iwf)*u*si
     &  +(2*a1(56,is,iwf)*u2+4*a1(57,is,iwf)*t2)*si2
     &  +(2*a1(59,is,iwf)*u3+4*a1(60,is,iwf)*u*t2)*si3
     &  +(2*a1(62,is,iwf)*u2+4*a1(63,is,iwf)*t2)*si
     &  +(2*a1(65,is,iwf)*u3+4*a1(66,is,iwf)*u*t2)*si2
     &  +(2*a1(68,is,iwf)*u3+4*a1(69,is,iwf)*u*t2)*si)
     &  +a1(40,is,iwf)*dtphi21+a1(44,is,iwf)*dtphi20
     &   +a1(45,is,iwf)*dtphi31+a21*dty21

        auu=2*(a1(4,is,iwf)+3*a1(10,is,iwf)*u+a1(13,is,iwf)*s
     &  + 6*a1(20,is,iwf)*u2
     &  +3*a1(23,is,iwf)*u*s+a1(32,is,iwf)*s2+a1(33,is,iwf)*t2)
     &  +6*a1(37,is,iwf)*u/rr-2*(a1(38,is,iwf)+3*a1(39,is,iwf)*u)*rrlog
     &  +a1(40,is,iwf)*d2uphi21+a1(44,is,iwf)*d2uphi20
     &  +a1(45,is,iwf)*d2uphi31+a21*d2uy21+d2u_termsa31*rrlog
     &  +6*a1(53,is,iwf)*u*si
     &  +(12*a1(55,is,iwf)*u2+2*a1(56,is,iwf)*t2)*si2
     &  +(20*a1(58,is,iwf)*u3+6*a1(59,is,iwf)*u*t2)*si3
     &  +(12*a1(61,is,iwf)*u2+2*a1(62,is,iwf)*t2)*si
     &  +(20*a1(64,is,iwf)*u3+6*a1(65,is,iwf)*u*t2)*si2
     &  +(20*a1(67,is,iwf)*u3+6*a1(68,is,iwf)*u*t2)*si

        assatt=2*(a1(5,is,iwf)+a1(6,is,iwf)
     &  + (a1(14,is,iwf)+a1(16,is,iwf))*u
     &  +(3*a1(11,is,iwf)+a1(18,is,iwf))*s
     &  +(a1(32,is,iwf)+a1(33,is,iwf))*u2
     &  +(6*a1(21,is,iwf)+a1(34,is,iwf))*s2
     &  +(6*a1(22,is,iwf)+a1(34,is,iwf))*t2
     &  +(3*a1(24,is,iwf)+a1(31,is,iwf))*u*s)
     &  +((a1(35,is,iwf)*u+3*a1(36,is,iwf)*s)*rr2
     &  +a1(37,is,iwf)*u3)/(2*rr3)
     &  +2*(a1(38,is,iwf)+a1(39,is,iwf)*u)*(two+rrlog)
     &  +a1(40,is,iwf)*(d2sphi21+d2tphi21)
     &  +a1(44,is,iwf)*(d2sphi20+d2tphi20)
     &  +a1(45,is,iwf)*(d2sphi31+d2tphi31)+a21*(d2sy21+d2ty21)
     &  +d2sd2t_termsa31*rrlog
     &  +2*(ds_termsa31*ds_rrlog+dt_termsa31*dt_rrlog*t2)
     &  +(2*a1(53,is,iwf)*u3+4*a1(54,is,iwf)*u*rr2)*si3
     &  +(6*a1(55,is,iwf)*u4+2*a1(56,is,iwf)*u2*(s2+3*t2)
     &  +6*a1(57,is,iwf)*t2*(2*s2+t2))*si4
     &  +(12*a1(58,is,iwf)*u5+2*a1(59,is,iwf)*u3*(s2+6*t2)
     &  +24*a1(60,is,iwf)*u*t2*rr2)*si5
     &  +(2*a1(61,is,iwf)*u4+4*a1(62,is,iwf)*u2*rr2
     &  +2*a1(63,is,iwf)*t2*(6*s2+t2))*si3
     &  +(6*a1(64,is,iwf)*u5+2*a1(65,is,iwf)*u3*(s2+3*t2)
     &  +6*a1(66,is,iwf)*u*t2*(2*s2+t2))*si4
     &  +(2*a1(67,is,iwf)*u5+4*a1(68,is,iwf)*u3*rr2
     &  +2*a1(69,is,iwf)*u*t2*(6*s2+t2))*si3

        aus=a1(7,is,iwf)+2*(a1(13,is,iwf)*u+a1(14,is,iwf)*s)
     &  + 3*(a1(23,is,iwf)*u2+a1(24,is,iwf)*s2)+a1(31,is,iwf)*t2
     &  +4*a1(32,is,iwf)*u*s+ s*terma2
     &  +a1(40,is,iwf)*d2suphi21+a1(44,is,iwf)*d2suphi20
     &  +a1(45,is,iwf)*d2suphi31+a21*d2suy21
     &  +d2us_termsa31*rrlog + du_termsa31*ds_rrlog
     &  -(3*a1(53,is,iwf)*u2+a1(54,is,iwf)*t2)*si2
     &  -(8*a1(55,is,iwf)*u3+4*a1(56,is,iwf)*u*t2)*si3
     &  -(15*a1(58,is,iwf)*u4+9*a1(59,is,iwf)*u2*t2
     &  +3*a1(60,is,iwf)*t4)*si4
     &  -(4*a1(61,is,iwf)*u3+2*a1(62,is,iwf)*u*t2)*si2
     &  -(10*a1(64,is,iwf)*u4+6*a1(65,is,iwf)*u2*t2
     &  +2*a1(66,is,iwf)*t4)*si3
     &  -(5*a1(67,is,iwf)*u4+3*a1(68,is,iwf)*u2*t2+a1(69,is,iwf)*t4)*si2

        aut=t*(2*(a1(16,is,iwf) + a1(31,is,iwf)*s+2*a1(33,is,iwf)*u)
     &  + terma2
     &  +d2ut_termsa31*rrlog + du_termsa31*dt_rrlog
     &  +2*a1(54,is,iwf)*si
     &  +4*a1(56,is,iwf)*u*si2
     &  +(6*a1(59,is,iwf)*u2+4*a1(60,is,iwf)*t2)*si3
     &  +4*a1(62,is,iwf)*u*si
     &  +(6*a1(65,is,iwf)*u2+4*a1(66,is,iwf)*t2)*si2
     &  +(6*a1(68,is,iwf)*u2+4*a1(69,is,iwf)*t2)*si)
     &  +a1(40,is,iwf)*d2utphi21+a1(44,is,iwf)*d2utphi20
     &  +a1(45,is,iwf)*d2utphi31
     &  +a21*d2uty21

        ast=t*(2*((a1(18,is,iwf) + a1(31,is,iwf)*u+2*a1(34,is,iwf)*s))
     &  +a1(35,is,iwf)*st35+a1(36,is,iwf)*st36+a1(37,is,iwf)*st37
     &  +a1(38,is,iwf)*st38
     &  +a1(39,is,iwf)*st39
     &  +d2st_termsa31*rrlog + termsa31*d2st_rrlog
     &  +ds_termsa31*dt_rrlog +dt_termsa31*ds_rrlog
     &  -2*a1(54,is,iwf)*u*si2
     &  -(4*a1(56,is,iwf)*u2+8*a1(57,is,iwf)*t2)*si3
     &  -(6*a1(59,is,iwf)*u3+12*a1(60,is,iwf)*u*t2)*si4
     &  -(2*a1(62,is,iwf)*u2+4*a1(63,is,iwf)*t2)*si2
     &  -(4*a1(65,is,iwf)*u3+8*a1(66,is,iwf)*u*t2)*si3
     &  -(2*a1(68,is,iwf)*u3+4*a1(69,is,iwf)*u*t2)*si2)
     &  +a1(40,is,iwf)*d2stphi21+a1(44,is,iwf)*d2stphi20
     &  +a1(45,is,iwf)*d2stphi31
     &  +a21*d2sty21

        aaa=       a1(1,is,iwf)*u + a1(2,is,iwf)*s
     &  + a1(4,is,iwf)*u2
     &  + a1(5,is,iwf)*s2  + a1(6,is,iwf)*t2 + a1(7,is,iwf)*u*s
     &                + a1(10,is,iwf)*u3 + a1(11,is,iwf)*s3
     &  + a1(13,is,iwf)*u2*s + a1(14,is,iwf)*u*s2
     &  + a1(16,is,iwf)*u*t2
     &                  + a1(18,is,iwf)*s*t2
     &  + a1(20,is,iwf)*u4 + a1(21,is,iwf)*s4 + a1(22,is,iwf)*t4
     &  + a1(23,is,iwf)*u3*s
     &  + a1(24,is,iwf)*u*s3
     &  + a1(31,is,iwf)*u*s*t2 + a1(32,is,iwf)*u2*s2
     &  +t2*(a1(33,is,iwf)*u2+a1(34,is,iwf)*s2)
     &  + rr*(a1(35,is,iwf)*u + a1(36,is,iwf)*s) + a1(37,is,iwf)*u3/rr
     &  + term38*(a1(38,is,iwf) + a1(39,is,iwf)*u)
     &  +a1(40,is,iwf)*phi21+a1(44,is,iwf)*phi20+a1(45,is,iwf)*phi31
     &  +a21*y21
     &  + termsa31*rrlog
     &  + (a1(53,is,iwf)*u3+a1(54,is,iwf)*u*t2)*si
     &  + (a1(55,is,iwf)*u4+a1(56,is,iwf)*u2*t2+a1(57,is,iwf)*t4)*si2
     &  + (a1(58,is,iwf)*u5+a1(59,is,iwf)*u3*t2+a1(60,is,iwf)*u*t4)*si3
     &  + (a1(61,is,iwf)*u4+a1(62,is,iwf)*u2*t2+a1(63,is,iwf)*t4)*si
     &  + (a1(64,is,iwf)*u5+a1(65,is,iwf)*u3*t2+a1(66,is,iwf)*u*t4)*si2
     &  + (a1(67,is,iwf)*u5+a1(68,is,iwf)*u3*t2+a1(69,is,iwf)*u*t4)*si

        bu=a2(1,is,iwf)+2*a2(4,is,iwf)*u+a2(7,is,iwf)
     &  *s+3*a2(10,is,iwf)*u2+2*a2(13,is,iwf)*u*s
     &  +a2(14,is,iwf)*s2+a2(16,is,iwf)*t2+ 4*a2(20,is,iwf)*u3
     &  +3*a2(23,is,iwf)*u2*s
     &  +a2(24,is,iwf)*s3+a2(31,is,iwf)*s*t2+2*u*(a2(32,is,iwf)*s2
     &  +a2(33,is,iwf)*t2)
     &  +a2(35,is,iwf)*rr+3*a2(37,is,iwf)*u2/rr
     &  +(-2*a2(38,is,iwf)*u+a2(39,is,iwf)*(rr2-3*u2))*rrlog
     &  +a2(40,is,iwf)*duphi21+a2(44,is,iwf)*duphi20
     &  +a2(45,is,iwf)*duphi31
     &  +du_termsb31*rrlog
     &  +(3*a2(53,is,iwf)*u2+a2(54,is,iwf)*t2)*si
     &  +(4*a2(55,is,iwf)*u3+2*a2(56,is,iwf)*u*t2)*si2
     &  +(5*a2(58,is,iwf)*u4+3*a2(59,is,iwf)*u2*t2+a2(60,is,iwf)*t4)*si3
     &  +(4*a2(61,is,iwf)*u3+2*a2(62,is,iwf)*u*t2)*si
     &  +(5*a2(64,is,iwf)*u4+3*a2(65,is,iwf)*u2*t2+a2(66,is,iwf)*t4)*si2
     &  +(5*a2(67,is,iwf)*u4+3*a2(68,is,iwf)*u2*t2+a2(69,is,iwf)*t4)*si

        bs=a2(2,is,iwf)+2*a2(5,is,iwf)*s+a2(7,is,iwf)*u
     &  +3*a2(11,is,iwf)*s2+a2(13,is,iwf)*u2
     &  +2*a2(14,is,iwf)*u*s+a2(18,is,iwf)*t2+4*a2(21,is,iwf)*s3
     &  +a2(23,is,iwf)*u3
     &  +3*a2(24,is,iwf)*u*s2+a2(31,is,iwf)*u*t2
     &  +2*s*(a2(32,is,iwf)*u2+a2(34,is,iwf)*t2)
     &  +a2(36,is,iwf)*(s2+t2/2)/rr +s*termb1
     &  +a2(40,is,iwf)*dsphi21+a2(44,is,iwf)*dsphi20+
     &  a2(45,is,iwf)*dsphi31
     &  +ds_termsb31*rrlog + termsb31*ds_rrlog
     &  -(a2(53,is,iwf)*u3+a2(54,is,iwf)*u*t2)*si2
     &  -2*(a2(55,is,iwf)*u4+a2(56,is,iwf)*u2*t2+a2(57,is,iwf)*t4)*si3
     &  -3*(a2(58,is,iwf)*u5+a2(59,is,iwf)*u3*t2+a2(60,is,iwf)*u*t4)*si4
     &  -  (a2(61,is,iwf)*u4+a2(62,is,iwf)*u2*t2+a2(63,is,iwf)*t4)*si2
     &  -2*(a2(64,is,iwf)*u5+a2(65,is,iwf)*u3*t2+a2(66,is,iwf)*u*t4)*si3
     &  -  (a2(67,is,iwf)*u5+a2(68,is,iwf)*u3*t2+a2(69,is,iwf)*u*t4)*si2

        bt=t*(2*(a2(6,is,iwf)+a2(16,is,iwf)*u+a2(18,is,iwf)*s
     &  +2*a2(22,is,iwf)*t2
     &  +a2(31,is,iwf)*u*s+a2(33,is,iwf)*u2+a2(34,is,iwf)*s2)
     &  +a2(36,is,iwf)*s/(2*rr) +termb1
     &  +dt_termsb31*rrlog + termsb31*dt_rrlog
     &  +2*a2(54,is,iwf)*u*si
     &  +(2*a2(56,is,iwf)*u2+4*a2(57,is,iwf)*t2)*si2
     &  +(2*a2(59,is,iwf)*u3+4*a2(60,is,iwf)*u*t2)*si3
     &  +(2*a2(62,is,iwf)*u2+4*a2(63,is,iwf)*t2)*si
     &  +(2*a2(65,is,iwf)*u3+4*a2(66,is,iwf)*u*t2)*si2
     &  +(2*a2(68,is,iwf)*u3+4*a2(69,is,iwf)*u*t2)*si)
     &  +a2(40,is,iwf)*dtphi21+a2(44,is,iwf)*dtphi20
     &  +a2(45,is,iwf)*dtphi31

        buu=2*(a2(4,is,iwf)+3*a2(10,is,iwf)*u+a2(13,is,iwf)*s
     &  + 6*a2(20,is,iwf)*u2
     &  +3*a2(23,is,iwf)*u*s+a2(32,is,iwf)*s2+a2(33,is,iwf)*t2)
     &  +6*a2(37,is,iwf)*u/rr-2*(a2(38,is,iwf)+3*a2(39,is,iwf)*u)*rrlog
     &  +a2(40,is,iwf)*d2uphi21+a2(44,is,iwf)*d2uphi20
     &  +a2(45,is,iwf)*d2uphi31
     &  +d2u_termsb31*rrlog
     &  +6*a2(53,is,iwf)*u*si
     &  +(12*a2(55,is,iwf)*u2+2*a2(56,is,iwf)*t2)*si2
     &  +(20*a2(58,is,iwf)*u3+6*a2(59,is,iwf)*u*t2)*si3
     &  +(12*a2(61,is,iwf)*u2+2*a2(62,is,iwf)*t2)*si
     &  +(20*a2(64,is,iwf)*u3+6*a2(65,is,iwf)*u*t2)*si2
     &  +(20*a2(67,is,iwf)*u3+6*a2(68,is,iwf)*u*t2)*si

        bssbtt=2*(a2(5,is,iwf)+a2(6,is,iwf)
     &  + (a2(14,is,iwf)+a2(16,is,iwf))*u
     &  +(3*a2(11,is,iwf)+a2(18,is,iwf))*s
     &  +(a2(32,is,iwf)+a2(33,is,iwf))*u2
     &  +(6*a2(21,is,iwf)+a2(34,is,iwf))*s2
     &  +(6*a2(22,is,iwf)+a2(34,is,iwf))*t2
     &  +(3*a2(24,is,iwf)+a2(31,is,iwf))*u*s)
     &  +((a2(35,is,iwf)*u+3*a2(36,is,iwf)*s)*rr2
     &  +a2(37,is,iwf)*u3)/(2*rr3)
     &  +2*(a2(38,is,iwf)+a2(39,is,iwf)*u)*(two+rrlog)
     &  +a2(40,is,iwf)*(d2sphi21+d2tphi21)
     &  +a2(44,is,iwf)*(d2sphi20+d2tphi20)
     &  +a2(45,is,iwf)*(d2sphi31+d2tphi31)
     &  +d2sd2t_termsb31*rrlog
     &  +2*(ds_termsb31*ds_rrlog+dt_termsb31*dt_rrlog*t2)
     &  +(2*a2(53,is,iwf)*u3+4*a2(54,is,iwf)*u*rr2)*si3
     &  +(6*a2(55,is,iwf)*u4+2*a2(56,is,iwf)*u2*(s2+3*t2)
     &  +6*a2(57,is,iwf)*t2*(2*s2+t2))
     &  *si4
     &  +(12*a2(58,is,iwf)*u5+2*a2(59,is,iwf)*u3*(s2+6*t2)
     &  +24*a2(60,is,iwf)*u*t2*rr2)
     &  *si5
     &  +(2*a2(61,is,iwf)*u4+4*a2(62,is,iwf)*u2*rr2
     &  +2*a2(63,is,iwf)*t2*(6*s2+t2))*si3
     &  +(6*a2(64,is,iwf)*u5+2*a2(65,is,iwf)*u3*(s2+3*t2)
     &  +6*a2(66,is,iwf)*u*t2*(2*s2+t2
     &  ))*si4
     &  +(2*a2(67,is,iwf)*u5+4*a2(68,is,iwf)*u3*rr2
     &  +2*a2(69,is,iwf)*u*t2*(6*s2+t2))*si3

        bus=a2(7,is,iwf)+2*(a2(13,is,iwf)*u+a2(14,is,iwf)*s)
     &  + 3*(a2(23,is,iwf)*u2+a2(24,is,iwf)*s2)+a2(31,is,iwf)*t2
     &  +4*a2(32,is,iwf)*u*s
     &  + s*termb2
     &  +a2(40,is,iwf)*d2suphi21+a2(44,is,iwf)*d2suphi20
     &  +a2(45,is,iwf)*d2suphi31
     &  +d2us_termsb31*rrlog + du_termsb31*ds_rrlog
     &  -(3*a2(53,is,iwf)*u2+a2(54,is,iwf)*t2)*si2
     &  -(8*a2(55,is,iwf)*u3+4*a2(56,is,iwf)*u*t2)*si3
     &  -(15*a2(58,is,iwf)*u4+12*a2(59,is,iwf)*u2*t2
     &  +3*a2(60,is,iwf)*t4)*si4
     &  -(4*a2(61,is,iwf)*u3+2*a2(62,is,iwf)*u*t2)*si2
     &  -(10*a2(64,is,iwf)*u4+6*a2(65,is,iwf)*u2*t2
     &  +2*a2(66,is,iwf)*t4)*si3
     &  -(5*a2(67,is,iwf)*u4+3*a2(68,is,iwf)*u2*t2
     &  +a2(69,is,iwf)*t4)*si2

        but=t*(2*(a2(16,is,iwf) + a2(31,is,iwf)*s+2*a2(33,is,iwf)*u)
     &  + termb2
     &  +d2ut_termsb31*rrlog + du_termsb31*dt_rrlog
     &  +2*a2(54,is,iwf)*si
     &  +4*a2(56,is,iwf)*u*si2
     &  +(6*a2(59,is,iwf)*u2+4*a2(60,is,iwf)*t2)*si3
     &  +4*a2(62,is,iwf)*u*si
     &  +(6*a2(65,is,iwf)*u2+4*a2(66,is,iwf)*t2)*si2
     &  +(6*a2(68,is,iwf)*u2+4*a2(69,is,iwf)*t2)*si)
     &  +a2(40,is,iwf)*d2utphi21+a2(44,is,iwf)*d2utphi20
     &  +a2(45,is,iwf)*d2utphi31

        bst=t*(2*((a2(18,is,iwf) + a2(31,is,iwf)*u+2*a2(34,is,iwf)*s))
     &  +a2(35,is,iwf)*st35+a2(36,is,iwf)*st36+a2(37,is,iwf)*st37
     &  +a2(38,is,iwf)*st38
     &  +a2(39,is,iwf)*st39
     &  +d2st_termsb31*rrlog + termsb31*d2st_rrlog
     &  +ds_termsb31*dt_rrlog +dt_termsb31*ds_rrlog
     &  -2*a2(54,is,iwf)*u*si2
     &  -(4*a2(56,is,iwf)*u2+8*a2(57,is,iwf)*t2)*si3
     &  -(6*a2(59,is,iwf)*u3+12*a2(60,is,iwf)*u*t2)*si4
     &  -(2*a2(62,is,iwf)*u2+4*a2(63,is,iwf)*t2)*si2
     &  -(4*a2(65,is,iwf)*u3+8*a2(66,is,iwf)*u*t2)*si3
     &  -(2*a2(68,is,iwf)*u3+4*a2(69,is,iwf)*u*t2)*si2)
     &  +a2(40,is,iwf)*d2stphi21+a2(44,is,iwf)*d2stphi20
     &  +a2(45,is,iwf)*d2stphi31

        bbb= one + a2(1,is,iwf)*u + a2(2,is,iwf)*s
     &  + a2(4,is,iwf)*u2
     &  + a2(5,is,iwf)*s2  + a2(6,is,iwf)*t2 + a2(7,is,iwf)*u*s
     &                + a2(10,is,iwf)*u3 + a2(11,is,iwf)*s3
     &  + a2(13,is,iwf)*u2*s + a2(14,is,iwf)*u*s2
     &  + a2(16,is,iwf)*u*t2
     &                  + a2(18,is,iwf)*s*t2
     &  + a2(20,is,iwf)*u4 + a2(21,is,iwf)*s4 + a2(22,is,iwf)*t4
     &  + a2(23,is,iwf)*u3*s
     &  + a2(24,is,iwf)*u*s3
     &  + a2(31,is,iwf)*u*s*t2 + a2(32,is,iwf)*u2*s2
     &  +t2*(a2(33,is,iwf)*u2+a2(34,is,iwf)*s2)
     &  + rr*(a2(35,is,iwf)*u + a2(36,is,iwf)*s) + a2(37,is,iwf)*u3/rr
     &  + term38*(a2(38,is,iwf) + a2(39,is,iwf)*u)
     &  +a2(40,is,iwf)*phi21+a2(44,is,iwf)*phi20+a2(45,is,iwf)*phi31
     &  + termsb31*rrlog
     &  + (a2(53,is,iwf)*u3+a2(54,is,iwf)*u*t2)*si
     &  + (a2(55,is,iwf)*u4+a2(56,is,iwf)*u2*t2+a2(57,is,iwf)*t4)*si2
     &  + (a2(58,is,iwf)*u5+a2(59,is,iwf)*u3*t2+a2(60,is,iwf)*u*t4)*si3
     &  + (a2(61,is,iwf)*u4+a2(62,is,iwf)*u2*t2+a2(63,is,iwf)*t4)*si
     &  + (a2(64,is,iwf)*u5+a2(65,is,iwf)*u3*t2+a2(66,is,iwf)*u*t4)*si2
     &  + (a2(67,is,iwf)*u5+a2(68,is,iwf)*u3*t2+a2(69,is,iwf)*u*t4)*si

        bbb2=bbb*bbb
        bbb3=bbb*bbb2

        fu=(bbb*au-aaa*bu)/bbb2
        fs=(bbb*as-aaa*bs)/bbb2
        ft=(bbb*at-aaa*bt)/bbb2

        fuu=(bbb*(bbb*auu-aaa*buu-2*au*bu)+2*aaa*bu*bu)/bbb3
        fssftt=(bbb*(bbb*assatt-aaa*bssbtt-2*(as*bs+at*bt))
     &        +2*aaa*(bs*bs+bt*bt))/bbb3
        fus=(bbb*(bbb*aus-aaa*bus-au*bs-bu*as)+2*aaa*bu*bs)/bbb3
        fut=(bbb*(bbb*aut-aaa*but-au*bt-bu*at)+2*aaa*bu*bt)/bbb3
        fst=(bbb*(bbb*ast-aaa*bst-as*bt-bs*at)+2*aaa*bs*bt)/bbb3

        fuu=fuu*dd1*dd1+fu*dd2
        fssftt=fssftt*(dd3**2+dd4**2)+4*fst*dd3*dd4+2*(fs*dd5+ft*dd6)
        gus=fus
        fus=(fus*dd3+fut*dd4)*dd1
        fut=(gus*dd4+fut*dd3)*dd1
        fu=fu*dd1
        gs=fs
        fs=fs*dd3+ft*dd4
        ft=gs*dd4+ft*dd3

        fsfti=(fs+ft)/r_en(i,ic)
        fsftj=(fs-ft)/r_en(j,ic)

        fubu=fu/rij
        wtjastrow = sspin*(aaa/bbb)
        fsum=fsum + wtjastrow

        v(1,i)=v(1,i) + sspin*(fsfti*rvec_en(1,i,ic)+fubu*rvec_ee(1,ij))
        v(2,i)=v(2,i) + sspin*(fsfti*rvec_en(2,i,ic)+fubu*rvec_ee(2,ij))
        v(3,i)=v(3,i) + sspin*(fsfti*rvec_en(3,i,ic)+fubu*rvec_ee(3,ij))
        v(1,j)=v(1,j) + sspin*(fsftj*rvec_en(1,j,ic)-fubu*rvec_ee(1,ij))
        v(2,j)=v(2,j) + sspin*(fsftj*rvec_en(2,j,ic)-fubu*rvec_ee(2,ij))
        v(3,j)=v(3,j) + sspin*(fsftj*rvec_en(3,j,ic)-fubu*rvec_ee(3,ij))

c Warning: The generalization to arbitrary dimension is to be checked but is most likely OK.
c We never use this for anything except 3d anyway.
c       d2ij= sspin*2*(fuu + fssftt    + 2*(fu/rij +
c    &  (ss*u2mt2*fus+tt*s2mu2*fut)/(rij*s2mt2) + 2*(ss*fs-tt*ft)/
c    &  s2mt2))
        d2ij= sspin*2*((ndim-1)*(fu/rij + (ss*fs-tt*ft)/s2mt2)
     &  + fuu + fssftt + (ss*u2mt2*fus+tt*s2mu2*fut)/(rij*s2mt2))
        div_vj(i)=div_vj(i)+d2ij/2
        div_vj(j)=div_vj(j)+d2ij/2
        d2=d2+d2ij

   40   continue
   60 continue

      value=fsum

      return
      end
c-----------------------------------------------------------------------

      subroutine psigm2(u,s,t,phi21,phi20,phi31)

      implicit real*8(a-h,o-z)

      r1=.5d0*(s+t)
      r2=.5d0*(s-t)
      call psigm(u,r1,r2,phi21,phi20,phi31,1)

      return
      end

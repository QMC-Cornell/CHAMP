      subroutine fock(u,s,t,rri,rrj,it)
! Written by Cyrus Umrigar
! Uses Fock expansion as described in:
! Fock's Expansion, Kato's Cusp Conditions and the Exponential Ansatz,
! C.R. Myers, C.J. Umrigar, J.P. Sethna and J.D. Morgan, PRA, 44, 5537 (1991).
      use constants_mod
      use contr2_mod
      use wfsec_mod
      use jaspar3_mod
      use pars_mod
      implicit real*8(a-h,o-z)

      common /gamder/ dsphi21,dtphi21,duphi21,d2sphi21,d2tphi21,
     &d2uphi21,d2stphi21,d2suphi21,d2utphi21,
     &                dsphi31,dtphi31,duphi31,d2sphi31,d2tphi31,
     &d2uphi31,d2stphi31,d2suphi31,d2utphi31,
     &                dsphi20,dtphi20,duphi20,d2sphi20,d2tphi20,
     &d2uphi20,d2stphi20,d2suphi20,d2utphi20,
     &                dsy21,dty21,duy21,d2sy21,d2ty21,
     &d2uy21,d2sty21,d2suy21,d2uty21
      common /focktmp/ fc,fcu,fcuu,fcs,fcss,fct,fctt,fcst,fcus,fcut

      dimension u(*),s(*),t(*)

      r2=(s(2)+t(2))/two
      r=dsqrt(r2)
!     r4=r2*r2
      ri=one/r
      ri2=ri*ri
      ri3=ri2*ri
      ri4=ri3*ri
      ri5=ri4*ri

      call psigm(u(1),rri,rrj,phi21,phi20,phi31,it)
      y21=r2-u(2)

! phi20-like terms

      t0=fck(6,it,iwf)*u(3)+fck(7,it,iwf)*s(3)
     &  +fck(8,it,iwf)*u(2)*s(1)+fck(9,it,iwf)*u(1)*s(2)
      t0u=three*fck(6,it,iwf)*u(2)+two*fck(8,it,iwf)*u(1)*s(1)
     &   +fck(9,it,iwf)*s(2)
      t0s=three*fck(7,it,iwf)*s(2)+fck(8,it,iwf)*u(2)
     &   +two*fck(9,it,iwf)*u(1)*s(1)
      t0us=two*fck(8,it,iwf)*u(1)+two*fck(9,it,iwf)*s(1)
      t0uu=six*fck(6,it,iwf)*u(1)+two*fck(8,it,iwf)*s(1)
      t0ss=six*fck(7,it,iwf)*s(1)+two*fck(9,it,iwf)*u(1)
      t1=(fck(4,it,iwf)*u(1)*ri-t0*ri3)/two
      t2=(fck(4,it,iwf)*ri-t0u*ri3)/two
      t3=-(fck(4,it,iwf)*u(1)*ri3-three*t0*ri5)/four

      fc=fc+(fck(4,it,iwf)*u(1)+fck(5,it,iwf)*s(1))*r+t0*ri
      fcu=fcu+fck(4,it,iwf)*r+t0u*ri
      fcuu=fcuu+t0uu*ri
      fcs=fcs+fck(5,it,iwf)*(s(2)+t(2)/two)*ri+s(1)*t1+t0s*ri
      fcss=fcss+fck(5,it,iwf)*s(1)*(s(2)+three*t(2)/two)*ri3/two
     &    +t1+s(2)*t3-s(1)*t0s*ri3+t0ss*ri
      fct=fct+fck(5,it,iwf)*s(1)*t(1)*ri/two+t(1)*t1
      fctt=fctt+fck(5,it,iwf)*s(3)*ri3/four+t1+t(2)*t3
      fcus=fcus+s(1)*t2+t0us*ri
      fcut=fcut+t(1)*t2
      fcst=fcst+fck(5,it,iwf)*t(3)*ri3/four+t(1)*(s(1)*t3-t0s*ri3/two)

! phi21
      fc=fc+fck(1,it,iwf)*phi21
      fcu=fcu+fck(1,it,iwf)*duphi21
      fcuu=fcuu+fck(1,it,iwf)*d2uphi21
      fcs=fcs+fck(1,it,iwf)*dsphi21
      fcss=fcss+fck(1,it,iwf)*d2sphi21
      fct=fct+fck(1,it,iwf)*dtphi21
      fctt=fctt+fck(1,it,iwf)*d2tphi21
      fcus=fcus+fck(1,it,iwf)*d2suphi21
      fcut=fcut+fck(1,it,iwf)*d2utphi21
      fcst=fcst+fck(1,it,iwf)*d2stphi21

      if(ifock.ge.2) then
! phi31-like terms

        t31=fck(10,it,iwf)*u(3)+fck(11,it,iwf)*s(3)
     &     +fck(12,it,iwf)*u(2)*s(1)+fck(13,it,iwf)*u(1)*s(2)
     &     +(fck(14,it,iwf)*u(1)+fck(15,it,iwf)*s(1))*t(2)
        du_t31=3*fck(10,it,iwf)*u(2)+2*fck(12,it,iwf)*u(1)*s(1)
     &        +fck(13,it,iwf)*s(2)+fck(14,it,iwf)*t(2)
        ds_t31=3*fck(11,it,iwf)*s(2)+fck(12,it,iwf)*u(2)
     &        +2*fck(13,it,iwf)*u(1)*s(1)+fck(15,it,iwf)*t(2)
        dt_t31=2*(fck(14,it,iwf)*u(1)+fck(15,it,iwf)*s(1))*t(1)
        d2u_t31=6*fck(10,it,iwf)*u(1)+2*fck(12,it,iwf)*s(1)
        d2s_t31=6*fck(11,it,iwf)*s(1)+2*fck(13,it,iwf)*u(1)
        d2t_t31=2*(fck(14,it,iwf)*u(1)+fck(15,it,iwf)*s(1))
        d2us_t31=2*(fck(12,it,iwf)*u(1)+fck(13,it,iwf)*s(1))
        d2ut_t31=2*fck(14,it,iwf)*t(1)
        d2st_t31=2*fck(15,it,iwf)*t(1)

        rp1=r+one
        rlog=two*dlog(r/rp1)
        ds_rlog=s(1)*ri2/rp1
        dt_rlog=t(1)*ri2/rp1
        d2s_rlog=ri2/rp1-s(2)*ri4*(two+three*r)/(two*rp1*rp1)
        d2t_rlog=ri2/rp1-t(2)*ri4*(two+three*r)/(two*rp1*rp1)
        d2st_rlog= -s(1)*t(1)*ri4*(two+three*r)/(two*rp1*rp1)

        fc=fc+t31*rlog
        fcu=fcu+du_t31*rlog
        fcuu=fcuu+d2u_t31*rlog
        fcs=fcs+ds_t31*rlog+t31*ds_rlog
        fcss=fcss+d2s_t31*rlog+2*ds_t31*ds_rlog+t31*d2s_rlog
        fct=fct+dt_t31*rlog+t31*dt_rlog
        fctt=fctt+d2t_t31*rlog+2*dt_t31*dt_rlog+t31*d2t_rlog
        fcus=fcus+d2us_t31*rlog+du_t31*ds_rlog
        fcut=fcut+d2ut_t31*rlog+du_t31*dt_rlog
        fcst=fcst+d2st_t31*rlog+t31*d2st_rlog+ds_t31*dt_rlog
     &      +dt_t31*ds_rlog

        if(ifock.ge.3) then
! phi31
          fc=fc+fck(3,it,iwf)*phi31
          fcu=fcu+fck(3,it,iwf)*duphi31
          fcuu=fcuu+fck(3,it,iwf)*d2uphi31
          fcs=fcs+fck(3,it,iwf)*dsphi31
          fcss=fcss+fck(3,it,iwf)*d2sphi31
          fct=fct+fck(3,it,iwf)*dtphi31
          fctt=fctt+fck(3,it,iwf)*d2tphi31
          fcus=fcus+fck(3,it,iwf)*d2suphi31
          fcut=fcut+fck(3,it,iwf)*d2utphi31
          fcst=fcst+fck(3,it,iwf)*d2stphi31

          if(ifock.eq.4) then
! phi20
            fc=fc+fck(2,it,iwf)*phi20
            fcu=fcu+fck(2,it,iwf)*duphi20
            fcuu=fcuu+fck(2,it,iwf)*d2uphi20
            fcs=fcs+fck(2,it,iwf)*dsphi20
            fcss=fcss+fck(2,it,iwf)*d2sphi20
            fct=fct+fck(2,it,iwf)*dtphi20
            fctt=fctt+fck(2,it,iwf)*d2tphi20
            fcus=fcus+fck(2,it,iwf)*d2suphi20
            fcut=fcut+fck(2,it,iwf)*d2utphi20
            fcst=fcst+fck(2,it,iwf)*d2stphi20

            fc=fc+a21*y21
            fcu=fcu+a21*duy21
            fcuu=fcuu+a21*d2uy21
            fcs=fcs+a21*dsy21
            fcss=fcss+a21*d2sy21
            fct=fct+a21*dty21
            fctt=fctt+a21*d2ty21
            fcus=fcus+a21*d2suy21
            fcut=fcut+a21*d2uty21
            fcst=fcst+a21*d2sty21
          endif
        endif
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine psigm(r12,r1,r2,phi21,phi20,phi31,it)
! Written by Cyrus Umrigar

      use constants_mod
      use atom_mod
      use pars_mod
      implicit real*8(a-h,o-z)

      common /gamder/ dsphi21,dtphi21,duphi21,d2sphi21,d2tphi21,
     &d2uphi21,d2stphi21,d2suphi21,d2utphi21,
     &                dsphi31,dtphi31,duphi31,d2sphi31,d2tphi31,
     &d2uphi31,d2stphi31,d2suphi31,d2utphi31,
     &                dsphi20,dtphi20,duphi20,d2sphi20,d2tphi20,
     &d2uphi20,d2stphi20,d2suphi20,d2utphi20,
     &                dsy21,dty21,duy21,d2sy21,d2ty21,
     &d2uy21,d2sty21,d2suy21,d2uty21


      Zfock=znuc(it)

!     the GAM wave functions are region-dependent, hence,
!     the definition of rg (greater) and rl (lesser)

      if(r1.ge.r2) then
	rg = r1
        rl = r2
        sign_t=one
       else
        rg = r2
        rl = r1
        sign_t=-one
      endif

      call GM (r12,rg,rl,phi21,phi20,phi31)

      if(sign_t.ne.one) then
        dtphi21=-dtphi21
        d2stphi21=-d2stphi21
        d2utphi21=-d2utphi21
        dtphi31=-dtphi31
        d2stphi31=-d2stphi31
        d2utphi31=-d2utphi31
        dtphi20=-dtphi20
        d2stphi20=-d2stphi20
        d2utphi20=-d2utphi20
        dty21=-dty21
        d2sty21=-d2sty21
        d2uty21=-d2uty21
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine GM(r12,r1,r2,phi21,phi20,phi31)

!     Written by C. Myers, 6/88, updated by C.M., J. Sethna and C. Umrigar
!     on 7/7/88, to remove artificial singularity in the s2 term associated
!     with collinearity of the three particles
!     Uses Fock expansion as described in:
!     Fock's Expansion, Kato's Cusp Conditions and the Exponential Ansatz,
!     C.R. Myers, C.J. Umrigar, J.P. Sethna and J.D. Morgan, PRA, 44, 5537 (1991).

!     GM calculates the helium ground state (singlet S) wave function
!     derived by Gottschalk and Maslen, J. Phys. A. 20, 2781 (1987);
!     parameters for this function are calculated in ~wfs/wfpars.f,
!     which should be called in the outermost driving program since it
!     need only be called once;
!     this function actually calculates not only the result in J. Phys. A.,
!     but also variants of that function generated by us,
!     distinguished by fflag (in common/fflags/) currently fflag = 6: GM,
!     fflag = 7: exp(GM - sub)

!     **Warning** The wavefucntion is evaluated correctly when the 2
!     electrons are collinear but not when u=0.
!     The derivatives are not evaluated correctly in the 2nd part of
!     the program even when the 2 electrons are collinear.

      use constants_mod
      use const_mod
      use contr2_mod
      use pars_mod
      implicit real*8(a-h,o-z)
      real*8 xlob,dlnrr2,x
      real*8 r1,r2,r12,s,t,rr2,term
      real*8 r,Y20,Y21,root,y,omega,alph,beta

      real*8 psi1,psi21,psi20,psi2a,psi2b,psi2c
     &,psi2d,psi31,p1p21,phi21,phi20,phi31,ss1,ss2
      real*8 osix,otwe

      common /gamder/ dsphi21,dtphi21,duphi21,d2sphi21,d2tphi21,
     &d2uphi21,d2stphi21,d2suphi21,d2utphi21,
     &                dsphi31,dtphi31,duphi31,d2sphi31,d2tphi31,
     &d2uphi31,d2stphi31,d2suphi31,d2utphi31,
     &                dsphi20,dtphi20,duphi20,d2sphi20,d2tphi20,
     &d2uphi20,d2stphi20,d2suphi20,d2utphi20,
     &                dsy21,dty21,duy21,d2sy21,d2ty21,
     &d2uy21,d2sty21,d2suy21,d2uty21

      parameter (osix = 1.d0/6.d0, otwe = 1.d0/12.d0)
      external xlob
      dlob(x)=-dlog(dcos(x))
      d2lob(x)=dtan(x)

!     various configuration-dependent variables:
!     r = hyperradius (or its sqrt, depending on notation)
!     Y00, Y20 and Y21 are solid harmonics (IC equivalents of HH)
!     y = sin(alph), where alph = 2*atan(r2/r1)
!     omega = cos(theta_{12})

      xm1 = -Zfock
      xm2 = -Zfock
      xm12 = one
      xms = -Zfock
      xma = zero

      s = r1 + r2
      t = r1 - r2
      rr2 = r1*r1 + r2*r2
      r = dsqrt(rr2)
!     Y00 = one
      Y20 = two*(r1*r1 - r2*r2)
      Y21 = rr2 - r12*r12

      psi1 = xm1*r1 + xm2*r2 + half*xm12*r12

      dlnrr2 = dlog(rr2)

      term = osix*xm12*(xma*Y20-two*xms*Y21)

      rp1=r+one
      dlnrr2 = two*dlog(r/rp1)

      psi21 = term*(one/pi-half) * dlnrr2

      phi21=psi21

      if(ifock.ge.3) then
        psi31 = (otwe*Zfock*(pi-two)/pi)
     &  *(Zfock*s*Y21-r12*(rr2-five*r12*r12/six))*dlnrr2

        p1p21 = psi1*psi21

        phi31=psi31-p1p21
      endif

      if(ianalyt_lap.eq.0) goto 10
! Analytic Fock derivatives

      s2=s*s
      s3=s2*s
      s4=s3*s
      s5=s4*s

      t2=t*t
      t3=t2*t
      t4=t3*t
      t5=t4*t
      t6=t5*t

      u=r12
      u2=u*u
      u3=u2*u
      u4=u3*u
      u5=u4*u
      u6=u5*u

      rr=r
      rr3=rr2*rr
      rr4=rr3*rr
      rr5=rr4*rr
      rr6=rr5*rr
      rr8=rr6*rr2
      rr10=rr8*rr2

      dsy21 = s

      dty21 = t

      duy21 = -2*u

      d2sy21 = 1

      d2ty21 = 1

      d2uy21 = -2

      d2sty21 = 0

      d2suy21 = 0

      d2uty21 = 0

      rlog=dlog(rr/rp1)

      dslogr = s/(2*rr2*rp1)

      dtlogr = t/(2*rr2*rp1)

      dulogr = 0

      d2slogr = one/(2*rr2*rp1) - s2*(two+three*rr)/(4*rr4*rp1*rp1)

      d2tlogr = one/(2*rr2*rp1) - t2*(two+three*rr)/(4*rr4*rp1*rp1)

      d2ulogr = 0

      d2stlogr = -s*t*(two+three*rr)/(4*rr4*rp1*rp1)

      d2sulogr = 0

      d2utlogr = 0

      term=-(pi-two)*Zfock/(three*pi)

! Analytic derivatives of phi21
      dsphi21 = term*(rlog*dsy21+dslogr*y21)

      dtphi21 = term*(rlog*dty21+dtlogr*y21)

      duphi21 = term*(rlog*duy21+dulogr*y21)

      d2sphi21 = term*(rlog*d2sy21+d2slogr*y21+2*dslogr*dsy21)

      d2tphi21 = term*(rlog*d2ty21+d2tlogr*y21+2*dtlogr*dty21)

      d2uphi21 = term*(rlog*d2uy21+d2ulogr*y21+2*dulogr*duy21)

      d2stphi21 = term*(rlog*d2sty21+d2stlogr*y21+
     & dslogr*dty21+dtlogr*dsy21)

      d2suphi21 = term*(rlog*d2suy21+d2sulogr*y21+
     & dslogr*duy21+dulogr*dsy21)

      d2utphi21 = term*(rlog*d2uty21+d2utlogr*y21+
     & dulogr*dty21+dtlogr*duy21)

! Analytic derivatives of phi31
      if(ifock.ge.3) then

        phi311 = term*(3*Zfock*s3 + 3*Zfock*s*t2 - 6*Zfock*s*u2 + u3)/12

        dsphi311 = term*(9*Zfock*s2 + 3*Zfock*t2 - 6*Zfock*u2)/12

        dsphi31 = (rlog*dsphi311+dslogr*phi311)

        dtphi311 = Zfock*s*t*term/2

        dtphi31 = (rlog*dtphi311+dtlogr*phi311)

        duphi311 = term*(-12*Zfock*s*u + 3*u2)/12

        duphi31 = (rlog*duphi311+dulogr*phi311)

        d2sphi311 = 3*Zfock*s*term/2

        d2sphi31 = (rlog*d2sphi311+d2slogr*phi311+
     &   2*dslogr*dsphi311)

        d2tphi311 = Zfock*s*term/2

        d2tphi31 = (rlog*d2tphi311+d2tlogr*phi311+
     &   2*dtlogr*dtphi311)

        d2uphi311 = term*(-12*Zfock*s + 6*u)/12

        d2uphi31 = (rlog*d2uphi311+d2ulogr*phi311+
     &   2*dulogr*duphi311)

        d2stphi311 = Zfock*t*term/2

        d2stphi31 = (rlog*d2stphi311+d2stlogr*phi311+
     &   dslogr*dtphi311+dtlogr*dsphi311)

        d2suphi311 = -(Zfock*term*u)

        d2suphi31 = (rlog*d2suphi311+d2sulogr*phi311+
     &   dslogr*duphi311+dulogr*dsphi311)

        d2utphi311 = 0

        d2utphi31 = (rlog*d2utphi311+d2utlogr*phi311+
     &   dulogr*dtphi311+dtlogr*duphi311)

      endif

  10  continue
      if(ifock.eq.4) then

        root2 = two*rr2 - r12*r12
        if(root2.lt.0.d0) then
          root = 0.d0
        else
          root = dsqrt(root2)
        endif

        root3=root2*root
        root4=root3*root
        root5=root4*root
        root6=root5*root

        omega = Y21/(two*r1*r2)
        if(dabs(omega).gt.1.d0) then
          omega = omega/dabs(omega)
        endif

        y = two*r1*r2/rr2
        if(y.gt.1.d0) then
           y = 1.d0
        endif

        alph = dasin(y)
        beta = dasin(y*omega)
        arg1 = half*(alph-beta)
        arg2 = half*(alph+beta)
        arg3 = half*(pi-alph+beta)
        arg4 = half*(pi-alph-beta)
        if(arg1.lt.zero .or. arg1.gt.half*pi) stop 'arg1'
        if(arg2.lt.zero .or. arg2.gt.half*pi) stop 'arg2'
        if(arg3.lt.zero .or. arg3.gt.half*pi) stop 'arg3'
        if(arg4.lt.zero .or. arg4.gt.half*pi) stop 'arg4'

        a = dlog((s + u)/rr)

        b1 = dlog((t + u)**2/(s*t + u*root))

        b22 = dlog((s2 - u2)/(-t2 + u2))

        b32 = dlog((s*t + u*root)/ (-(s*t) + u*root))

        c1 = u*root

        b4 = two*(xlob(arg1)-xlob(arg2)+xlob(arg3)-xlob(arg4))

        b = b1+(alph*b22-beta*b32+b4)/Pi

        ss2 = b4

        term =r12*root + half*Y20
        if(term.le.0.d0) then
           ss1 = 0.d0
         else
           ss1 = dlog(term)
           ss2 = ss2 - two*beta*ss1
        endif

        term=r12*r12 - t*t
        if(term.gt.0.d0) ss2 = ss2 + (beta-alph)*dlog(term)

        term = s*s - r12*r12
        if(term.gt.0.d0) ss2 = ss2 + (beta+alph)*dlog(term)

        term = osix*xm12*(xma*Y20-two*xms*Y21)
        psi2a =-term*dlog((s+r12)/r)

        psi2b =  otwe*xm12*(two*xma*Y21-xms*Y20)
     &     * (two*dlog(t+r12) - (ss1 - ss2/pi))

        psi2c = -osix*xm12*xms*r12*root*(one+(two/pi)*beta)

        psi2d = otwe * ( -two*etrial*rr2
     &     + r12*r12 * (two*(xm1*xm1+xm2*xm2) + xm12*xm12)
     &     + eight*xm12*r12*(xm1*r1 + xm2*r2)
     &     + four*xm2*r1*r2*(three*xm1-xm12) )

        psi20 = psi2a + psi2b + psi2c + psi2d

        phi20=psi20-half*psi1*psi1

        if(ianalyt_lap.eq.0) return
! Analytic Fock derivatives

        dsy20 = 2*t

        dsalph = t/RR2

        dsbeta = s*u2/(RR4*dsqrt(2*u2/RR2 - u4/RR4))

        dsa = (t2 - s*u)/(2*RR2*(s + u))

        dsb1 = -((Root*t + s*u)/(Root*(s*t + Root*u)))

        dsb22 = 2*s/(s2 - u2)

        dsb32 = (-2*t3*u + 2*t*u3)/(Root*s2*t2 - Root3*u2)

        dsb2 = alph*dsb22+dsalph*b22

        dsb3 =-beta*dsb32-dsbeta*b32

        dsarg1 = dsalph-dsbeta

        dsarg2 = dsalph+dsbeta

        dsb4 = (dlob(arg1)-dlob(arg3))*dsarg1 +
     &      (dlob(arg4)-dlob(arg2))*dsarg2

        dsb = dsb1+(1/Pi)*(dsb2+dsb3+dsb4)

        dsc1 = s*u/Root

        dspd = -2*etrial*s + 2*Zfock*(1 + 3*Zfock)*s - 8*Zfock*u

        dspe = Zfock*(-(Zfock*s) + u/2)

        dsphi20 = otwe*(-4*Zfock*(y21*dsa+dsy21*a) +Zfock*(y20*dsb+dsy20*b)
     &   +(4*Zfock/Pi)*(c1*dsbeta+dsc1*(beta+half*Pi)) +dspd) +dspe

        dty20 = 2*s

        dtalph = -(s/RR2)

        dtbeta = t*u2/(RR4*dsqrt(2*u2/RR2 - u4/RR4))

        dta = -t/(2*RR2)

        dtb1 =
     &    (Root*s*t - Root*s*u + 2*s2*u + t2*u - t*u2 - 2*u3)/
     &     (Root*(t + u)*(s*t + Root*u))

        dtb22 = 2*t/(-t2 + u2)

        dtb32 = (-2*s3*u + 2*s*u3)/(Root*s2*t2 - Root3*u2)

        dtb2 = alph*dtb22+dtalph*b22

        dtb3 =-beta*dtb32-dtbeta*b32

        dtarg1 = dtalph-dtbeta

        dtarg2 = dtalph+dtbeta

        dtb4 = (dlob(arg1)-dlob(arg3))*dtarg1 +
     &      (dlob(arg4)-dlob(arg2))*dtarg2

        dtb = dtb1+(1/Pi)*(dtb2+dtb3+dtb4)

        dtc1 = t*u/Root

        dtpd = -2*etrial*t - 2*Zfock*(1 + 3*Zfock)*t

        dtpe = 0

        dtphi20 = otwe*(-4*Zfock*(y21*dta+dty21*a) +Zfock*(y20*dtb+dty20*b)
     &   +(4*Zfock/Pi)*(c1*dtbeta+dtc1*(beta+half*Pi)) +dtpd) +dtpe

        duy20 = 0

        dualph = 0

        dubeta = -2*u/(RR2*dsqrt(2*u2/RR2 - u4/RR4))

        dua = 1/(s + u)

        dub1 =
     &    (2*Root*s*t - s2*t - t3 + s2*u + t2*u + 2*t*u2)/
     &   (Root*(t + u)*(s*t + Root*u))

        dub22 = 2*(s - t)*(s + t)*u/((-s + u)*(s + u)*(-t + u)*(t + u))

        dub32 =
     &    -2*s*t*(s2 + t2 - 2*u2)/
     &     (Root*(-(s*t) + Root*u)*(s*t + Root*u))

        dub2 = alph*dub22+dualph*b22

        dub3 =-beta*dub32-dubeta*b32

        duarg1 = dualph-dubeta

        duarg2 = dualph+dubeta

        dub4 = (dlob(arg1)-dlob(arg3))*duarg1 +
     &      (dlob(arg4)-dlob(arg2))*duarg2

        dub = dub1+(1/Pi)*(dub2+dub3+dub4)

        duc1 = 2*(RR - u)*(RR + u)/Root

        dupd = -8*Zfock*s + 2*(1 + 4*Zfock**2)*u

        dupe = -(-(Zfock*s) + u/2)/2

        duphi20 = otwe*(-4*Zfock*(y21*dua+duy21*a) +Zfock*(y20*dub+duy20*b)
     &   +(4*Zfock/Pi)*(c1*dubeta+duc1*(beta+half*Pi)) +dupd) +dupe

        d2sy20 = 0

        d2salph = -(s*t/RR4)

        d2sbeta =
     &    (-6*RR2*s2*u4 + s4*u4 + 2*RR2*t2*u4 +
     &       s2*t2*u4 + s2*u6 - t2*u6)/
     &     (2*RR10*(2*u2/RR2 - u4/RR4)**1.5d0)

        d2sa =
     &    -(3*s2*t2 + t4 - 2*s3*u + 2*s*t2*u - s2*u2 + t2*u2)/
     &     (4*RR4*(s + u)**2)

        d2sb1 =
     &    (Root*t + s*u)**2/(Root2*(s*t + Root*u)**2) -
     &     (t2*u - u3)/(Root3*(s*t + Root*u))

        d2sb22 = -2*(s2 + u2)/((s - u)**2*(s + u)**2)

        d2sb32 =
     &    2*s*t*(t - u)**2*u*(t + u)**2*(3*s2 + 2*t2 - 3*u2)/
     &     (Root3*(-(s*t) + Root*u)**2*(s*t + Root*u)**2)

        d2sb2 = alph*d2sb22+d2salph*b22 + 2*dsalph*dsb22

        d2sb3 =-beta*d2sb32-d2sbeta*b32 - 2*dsbeta*dsb32

        d2sarg1 = d2salph-d2sbeta

        d2sarg2 = d2salph+d2sbeta

        d2sb4 = (dlob(arg1)-dlob(arg3))*d2sarg1 +
     &       (dlob(arg4)-dlob(arg2))*d2sarg2 +
     &   half*((d2lob(arg1)+d2lob(arg3))*dsarg1**2 +
     &        (-d2lob(arg4)-d2lob(arg2))*dsarg2**2)

        d2sb = d2sb1+(1/Pi)*(d2sb2+d2sb3+d2sb4)

        d2sc1 = -(u*(-t + u)*(t + u)/Root3)

        d2spd = -2*etrial + 2*Zfock*(1 + 3*Zfock)

        d2spe = -(Zfock**2)

        d2sphi20 = otwe*(-4*Zfock*(y21*d2sa+d2sy21*a+2*dsy21*dsa) +
     &   Zfock*(y20*d2sb+d2sy20*b+2*dsy20*dsb) +
     &   (4*Zfock/Pi)*(c1*d2sbeta+d2sc1*(beta+half*Pi)+2*dsc1*dsbeta)
     &   +d2spd)+d2spe

        d2ty20 = 0

        d2talph = s*t/RR4

        d2tbeta =
     &    (2*RR2*s2*u4 - 6*RR2*t2*u4 + s2*t2*u4 +
     &     t4*u4 - s2*u6 + t2*u6)/
     &     (2*RR10*(2*u2/RR2 - u4/RR4)**1.5d0)

        d2ta = -((s - t)*(s + t))/(4*RR4)

        d2tb1 =
     &    -((Root3*s2*t2 + 2*Root4*s*t*u - 2*Root3*s2*t*u +
     &         2*Root2*s3*t*u + s3*t3*u - 2*Root5*u2 -
     &         2*Root4*s*u2 + 3*Root3*s2*u2 -
     &         2*Root2*s3*u2 + 4*s5*u2 + 7*Root3*t2*u2 -
     &         2*Root2*s*t2*u2 - 2*Root*s2*t2*u2 +
     &         6*s3*t2*u2 - 4*Root*t4*u2 + 4*Root3*t*u3 -
     &         2*Root*s2*t*u3 - 3*s3*t*u3 - 6*Root*t3*u3 -
     &         5*s*t3*u3 - 3*Root3*u4 + 2*Root2*s*u4 -
     &         8*s3*u4 - 6*s*t2*u4 + 2*Root*t*u5 + 3*s*t*u5 +
     &         4*s*u6)/(Root3*(t + u)**2*(s*t + Root*u)**2))

        d2tb22 = 2*(t2 + u2)/((t - u)**2*(t + u)**2)

        d2tb32 =
     &    2*s*t*(s - u)**2*u*(s + u)**2*(2*s2 + 3*t2 - 3*u2)/
     &     (Root3*(-(s*t) + Root*u)**2*(s*t + Root*u)**2)

        d2tb2 = alph*d2tb22+d2talph*b22 + 2*dtalph*dtb22

        d2tb3 =-beta*d2tb32-d2tbeta*b32 - 2*dtbeta*dtb32

        d2targ1 = d2talph-d2tbeta

        d2targ2 = d2talph+d2tbeta

        d2tb4 = (dlob(arg1)-dlob(arg3))*d2targ1 +
     &       (dlob(arg4)-dlob(arg2))*d2targ2 +
     &   half*((d2lob(arg1)+d2lob(arg3))*dtarg1**2 +
     &        (-d2lob(arg4)-d2lob(arg2))*dtarg2**2)

        d2tb = d2tb1+(1/Pi)*(d2tb2+d2tb3+d2tb4)

        d2tc1 = -(u*(-s + u)*(s + u)/Root3)

        d2tpd = -2*etrial - 2*Zfock*(1 + 3*Zfock)

        d2tpe = 0

        d2tphi20 = otwe*(-4*Zfock*(y21*d2ta+d2ty21*a+2*dty21*dta) +
     &   Zfock*(y20*d2tb+d2ty20*b+2*dty20*dtb) +
     &   (4*Zfock/Pi)*(c1*d2tbeta+d2tc1*(beta+half*Pi)+2*dtc1*dtbeta)
     &   +d2tpd)+d2tpe

        d2uy20 = 0

        d2ualph = 0

        d2ubeta = -2*u4/(RR6*(2*u2/RR2 - u4/RR4)**1.5d0)

        d2ua = -one/(s + u)**2

        d2ub1 =
     &    -2*(2*Root*s*t - s2*t - t3 + s2*u + t2*u + 2*t*u2)/
     &      (Root*(t + u)**2*(s*t + Root*u)) +
     &     (2*RR2 - 2*u2)*(2*Root*s*t - s2*t - t3 + s2*u + t2*u +
     &         2*t*u2)/(Root2*(t + u)*(s*t + Root*u)**2) +
     &     (2*Root5*t2 + 2*Root3*s2*t2 - 4*s5*t2 -
     &        8*s3*t4 - 4*s*t6 + 3*s3*t3*u + 3*s*t5*u -
     &        Root3*t2*u2 + 18*s3*t2*u2 + 18*s*t4*u2 +
     &        2*Root3*t*u3 + 7*s3*t*u3 + 5*s*t3*u3 +
     &        3*Root3*u4 + 3*Root*t2*u4 - 12*s*t2*u4 +
     &        6*Root*t*u5 - 6*s*t*u5 + 3*Root*u6)/
     &      (Root3*(t + u)**2*(s*t + Root*u)**2)

        d2ub22 =
     &    2*(s - t)*(s + t)*(s2*t2 + s2*u2 + t2*u2 - 3*u4)/
     &     ((-s + u)**2*(s + u)**2*(-t + u)**2*(t + u)**2)

        d2ub32 =
     &    2*s*t*(2*Root5*s*t - 4*RR2*Root*s3*t -
     &        4*RR2*Root*s*t3 + 2*Root6*u - 3*s4*t2*u -
     &        3*s2*t4*u + 8*RR2*Root*s*t*u2 - Root3*s*t*u2 +
     &        Root*s3*t*u2 + Root*s*t3*u2 - Root4*u3 +
     &        2*s2*t2*u3 - 3*Root*s*t*u4 + 3*Root2*u5)/
     &     (Root3*(-(s*t) + Root*u)**2*(s*t + Root*u)**2)

        d2ub2 = alph*d2ub22+d2ualph*b22 + 2*dualph*dub22

        d2ub3 =-beta*d2ub32-d2ubeta*b32 - 2*dubeta*dub32

        d2uarg1 = d2ualph-d2ubeta

        d2uarg2 = d2ualph+d2ubeta

        d2ub4 = (dlob(arg1)-dlob(arg3))*d2uarg1 +
     &       (dlob(arg4)-dlob(arg2))*d2uarg2 +
     &   half*((d2lob(arg1)+d2lob(arg3))*duarg1**2 +
     &        (-d2lob(arg4)-d2lob(arg2))*duarg2**2)

        d2ub = d2ub1+(1/Pi)*(d2ub2+d2ub3+d2ub4)

        d2uc1 = u*(-3*s2 - 3*t2 + 2*u2)/Root3

        d2upd = 2*(1 + 4*Zfock**2)

        d2upe = -half*half

        d2uphi20 = otwe*(-4*Zfock*(y21*d2ua+d2uy21*a+2*duy21*dua) +
     &   Zfock*(y20*d2ub+d2uy20*b+2*duy20*dub) +
     &   (4*Zfock/Pi)*(c1*d2ubeta+d2uc1*(beta+half*Pi)+2*duc1*dubeta)
     &   +d2upd)+d2upe

        d2sty20 = 2

        d2stalph = (s - t)*(s + t)/(2*RR4)

        d2stbeta =
     &    (-8*RR2*s*t*u4 + s3*t*u4 + s*t3*u4 + 2*s*t*u6)/
     &     (2*RR10*(2*u2/RR2 - u4/RR4)**1.5d0)

        d2sta = s*t/(2*RR4)

        d2stb1 =
     &    (Root3*s*t - Root*s3*t - Root*s*t3 + s2*t2*u +
     &       3*Root*s*t*u2 + s2*u3 + t2*u3 - u5)/
     &     (Root3*(s*t + Root*u)**2)

        d2stb22 = 0

        d2stb32 =
     &    (-6*t2*u + 2*u3)/(Root*s2*t2 - Root3*u2) -
     &     (-2*t3*u + 2*t*u3)*
     &       (2*s4*t + 3*s2*t3 - 5*s2*t*u2 - 3*t3*u2 + 3*t*u4)/
     &      (Root*(Root*s2*t2 - Root3*u2)**2)

        d2stb2 = alph*d2stb22+d2stalph*b22 + dsalph*dtb22 + dtalph*dsb22

        d2stb3 =-beta*d2stb32-d2stbeta*b32 - dsbeta*dtb32 - dtbeta*dsb32

        d2starg1 = d2stalph-d2stbeta

        d2starg2 = d2stalph+d2stbeta

        d2stb4 = (dlob(arg1)-dlob(arg3))*d2starg1 +
     &        (dlob(arg4)-dlob(arg2))*d2starg2 +
     &   half*((d2lob(arg1)+d2lob(arg3))*dsarg1*dtarg1 +
     &        (-d2lob(arg4)-d2lob(arg2))*dsarg2*dtarg2)

        d2stb = d2stb1+(1/Pi)*(d2stb2+d2stb3+d2stb4)

        d2stc1 = -(s*t*u/Root3)

        d2stpd = 0

        d2stpe = 0

        d2stphi20 = otwe*(-4*Zfock*(y21*d2sta+d2sty21*a+dsy21*dta+dty21*dsa)
     &   +Zfock*(y20*d2stb+d2sty20*b+dsy20*dtb+dty20*dsb) +
     &   (4*Zfock/Pi)*
     &   (c1*d2stbeta+d2stc1*(beta+half*Pi)+dsc1*dtbeta+dtc1*dsbeta)
     &   +d2stpd) +d2stpe

        d2suy20 = 0

        d2sualph = 0

        d2subeta =
     &    (4*RR2*s*u3 - s3*u3 - s*t2*u3)/
     &     (RR8*(2*u2/RR2 - u4/RR4)**1.5d0)

        d2sua = -one/(s + u)**2

        d2sub1 =
     &    (s2*t3 + t5 - Root3*s*u + Root*s3*u + Root*s*t2*u -
     &       3*s2*t*u2 - 3*t3*u2 - 3*Root*s*u3 + 2*t*u4)/
     &     (Root3*(s*t + Root*u)**2)

        d2sub22 = 4*s*u/(s2 - u2)**2

        d2sub32 =
     &    (-2*t3 + 6*t*u2)/(Root*s2*t2 - Root3*u2) -
     &     (-2*t3*u + 2*t*u3)*
     &       (-2*s4*u - 5*s2*t2*u - 2*t4*u + 7*s2*u3 +
     &         7*t2*u3 - 5*u5)/(Root*(Root*s2*t2 - Root3*u2)**2)

        d2sub2 = alph*d2sub22+d2sualph*b22 + dsalph*dub22 + dualph*dsb22

        d2sub3 =-beta*d2sub32-d2subeta*b32 - dsbeta*dub32 - dubeta*dsb32

        d2suarg1 = d2sualph-d2subeta

        d2suarg2 = d2sualph+d2subeta

        d2sub4 = (dlob(arg1)-dlob(arg3))*d2suarg1 +
     &        (dlob(arg4)-dlob(arg2))*d2suarg2 +
     &   half*((d2lob(arg1)+d2lob(arg3))*dsarg1*duarg1 +
     &        (-d2lob(arg4)-d2lob(arg2))*dsarg2*duarg2)

        d2sub = d2sub1+(1/Pi)*(d2sub2+d2sub3+d2sub4)

        d2suc1 = s*(s2 + t2)/Root3

        d2supd = -8*Zfock

        d2supe = Zfock/2

        d2suphi20 = otwe*(-4*Zfock*(y21*d2sua+d2suy21*a+dsy21*dua+duy21*dsa)
     &   +Zfock*(y20*d2sub+d2suy20*b+dsy20*dub+duy20*dsb) +
     &   (4*Zfock/Pi)*
     &   (c1*d2subeta+d2suc1*(beta+half*Pi)+dsc1*dubeta+duc1*dsbeta)
     &   +d2supd) +d2supe

        d2uty20 = 0

        d2utalph = 0

        d2utbeta =
     &    (4*RR2*t*u3 - s2*t*u3 - t3*u3)/
     &     (RR8*(2*u2/RR2 - u4/RR4)**1.5d0)

        d2uta = 0

        d2utb1 =
     &    (-(Root*s2) + 2*s3 - 3*Root*t2 + 4*s*t2 + 2*Root*t*u +
     &        2*Root*u2 - 2*s*u2)/(Root2*(t + u)*(s*t + Root*u)) -
     &     (2*Root*s*t - s2*t - t3 + s2*u + t2*u + 2*t*u2)/
     &      (Root*(t + u)**2*(s*t + Root*u)) -
     &     t*(2*Root*s*t - s2*t - t3 + s2*u + t2*u + 2*t*u2)/
     &      (Root3*(t + u)*(s*t + Root*u)) -
     &     (Root*s + t*u)*(2*Root*s*t - s2*t - t3 + s2*u + t2*u +
     &         2*t*u2)/(Root2*(t + u)*(s*t + Root*u)**2)

        d2utb22 = -4*t*u/((t - u)**2*(t + u)**2)

        d2utb32 =
     &    2*s*(-s + u)**2*(s + u)**2*
     &      (-(s2*t2) - t4 - s2*u2 - t2*u2 + 2*u4)/
     &     (Root3*(-(s*t) + Root*u)**2*(s*t + Root*u)**2)

        d2utb2 = alph*d2utb22+d2utalph*b22 + dualph*dtb22 + dtalph*dub22

        d2utb3 =-beta*d2utb32-d2utbeta*b32 - dubeta*dtb32 - dtbeta*dub32

        d2utarg1 = d2utalph-d2utbeta

        d2utarg2 = d2utalph+d2utbeta

        d2utb4 = (dlob(arg1)-dlob(arg3))*d2utarg1 +
     &        (dlob(arg4)-dlob(arg2))*d2utarg2 +
     &   half*((d2lob(arg1)+d2lob(arg3))*duarg1*dtarg1 +
     &        (-d2lob(arg4)-d2lob(arg2))*duarg2*dtarg2)

        d2utb = d2utb1+(1/Pi)*(d2utb2+d2utb3+d2utb4)

        d2utc1 = t*(s2 + t2)/Root3

        d2utpd = 0

        d2utpe = 0

        d2utphi20 = otwe*(-4*Zfock*(y21*d2uta+d2uty21*a+duy21*dta+dty21*dua)
     &   +Zfock*(y20*d2utb+d2uty20*b+duy20*dtb+dty20*dub) +
     &   (4*Zfock/Pi)*
     &   (c1*d2utbeta+d2utc1*(beta+half*Pi)+duc1*dtbeta+dtc1*dubeta)
     &   +d2utpd) +d2utpe

      endif

      return
      end
!-----------------------------------------------------------------------

      function xlob(xeval)
! Written by Chris Myers and Cyrus Umrigar
! Uses Fock expansion as described in:
! Fock's Expansion, Kato's Cusp Conditions and the Exponential Ansatz,
! C.R. Myers, C.J. Umrigar, J.P. Sethna and J.D. Morgan, PRA, 44, 5537 (1991).

      use rlobxy_mod
!     implicit none
      implicit real*8(a-h,o-z)
      real*8 xeval,yeval,xpi,xlobpi,eps2,eps22,ysing7,xlob
      real*8 xsp,xpow,xpii,xdum,A,x
      real*8 one,two,half,pi,pib2,ln2,o18,o900,o19845
      integer nlob
      parameter (one = 1.d0, two = 2.d0, half = .5d0)
      parameter (pi = 3.14159265358979323846d0, pib2 = half*pi)
      parameter (ln2 = 0.69314718055994530941d0)
      parameter (xsp = 1.29d0, xpow = 1.32d0)
      parameter (o18 = .05555555555555555d0, o900 = .1111111111111111d-2
     & ,o19845 = 1.d0/19845.d0)

      A(x) = 1.d0 + x*x*(-3.d0 + x*(2.d0))

      xlobpi = pib2*ln2

      if(xeval.lt.0.d0) nlob = int((xeval-one)/pi)
      if(xeval.ge.0.d0) nlob = int(xeval/pi)
      xpi = xeval - pi*nlob
      xpii = xpi
      if(xpi.gt.pib2) xpii = pi - xpi

      if(xpii.le.xpow)call splint(rlobx,rloby,rloby2,nsplin,xpii,yeval)

      if(xpi.le.xsp) then
	 xlob = yeval
	 goto 500
      endif

      if(xpi.gt.pi-xsp) then
         xlob = two*xlobpi - yeval
         goto 500
      endif

      eps2 = pib2 - xpii
      eps22 = eps2*eps2
!     ysing7 = xlobpi + eps2*(dlog(dabs(eps2)) - one
!    &   + eps2*eps2*(-(one/18.d0) + eps2*eps2*(-(one/900.d0)
!    &   + eps2*eps2*(-(one/19845.d0)))))
      ysing7 = xlobpi + eps2 * (dlog(dabs(eps2)) - one
     &   - eps22*(o18 + eps22*(o900 + eps22*o19845)) )
      if(eps2.eq.0.d0) ysing7 = xlobpi

      if(xpi.gt.xpow.and.xpi.le.pib2) then
         xlob = ysing7
	 goto 500
      endif

      if(xpi.gt.pib2.and.xpi.le.pi-xpow) then
         xlob = two*xlobpi - ysing7
         goto 500
      endif

      if(xpi.gt.xsp.and.xpi.le.xpow) then
         xdum = (xpii - xsp)/(xpow - xsp)
         xlob = A(xdum)*yeval + (one-A(xdum))*ysing7
         goto 500
      endif

         xdum = (xpii - xsp)/(xpow - xsp)
         xlob = two*xlobpi - (A(xdum)*yeval + (one-A(xdum))*ysing7)

500   xlob = xlob + nlob*two*xlobpi

      return
      end
!-----------------------------------------------------------------------

      subroutine scale3(iwf,it)
! Written by Cyrus Umrigar
! Uses Fock expansion as described in:
! Fock's Expansion, Kato's Cusp Conditions and the Exponential Ansatz,
! C.R. Myers, C.J. Umrigar, J.P. Sethna and J.D. Morgan, PRA, 44, 5537 (1991).

      use atom_mod
      use contr2_mod
      use jaspar3_mod
      implicit real*8(a-h,o-z)


      common /focsav/ c4sav,c5sav,c7sav,c9sav

!JT      parameter(half=0.5d0,zero=0.d0,two=2.d0,four=4.d0)
      parameter(pi=3.141592653589793d0,const2=-(pi-2.d0)/(6.d0*pi))

! Focks's terms are evaluated using scaled variables. Scaled variables in
! phi21 originate terms in the real variables that contribute to phi31.
! The coeffs fck(10,1)-fck(15,1) are defined to cancel this contribution
! from phi31.
! If wf is a Pade', for small hyper-rad wf=P(1)(1-P(2)+..) and the 2nd
! order contribution from P(1)P(2) etc. has also to be subtracted out.

! WARNING*** we think that ijas2 is wrong in the phi31-like terms
! (scalel instead of two*scalel)

      if(isc.eq.2) then
        scalel=half*scalek(iwf)
       elseif(isc.eq.4) then
        scalel=scalek(iwf)
      endif
      const3=const2*znuc(it)

      fck(10,it,iwf)=-const3*two*scalel*fck(1,it,iwf)
      fck(11,it,iwf)= const3*two*scalel*fck(1,it,iwf)/four
      fck(12,it,iwf)= zero
      fck(13,it,iwf)= zero
      fck(14,it,iwf)= zero
      fck(15,it,iwf)= fck(11,it,iwf)+const3*two*scalel*fck(1,it,iwf)/two

      return

      entry scale20(iwf,it)

! Scaled variables in phi10 originate a 2nd order that contributes to phi20.
! The coeffs c(4,7,9) need to be redifined.

      c(4,it,iwf)=c4sav+fck(2,it,iwf)*(
     & scalel*(c(1,it,iwf)+b(1,1,iwf))+b(1,1,iwf)*b(2,1,iwf)-c4sav)
      c(7,it,iwf)=c7sav+fck(2,it,iwf)*(half*
     &(scalel*(c(2,it,iwf)+a(1,iwf))+a(1,iwf)*a(2,iwf))-c7sav)
      c(9,it,iwf)=c9sav+fck(2,it,iwf)*(half*
     &(scalel*(c(2,it,iwf)+a(1,iwf))+a(1,iwf)*a(2,iwf))-c9sav)

      return
      end

      subroutine cuspcheck2(scalek,a,b,diff,isp,nspin1,ncuspc,iprin)
c Written by Cyrus Umrigar
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c impose cusp-conditions via a penalty.
c note that for the log terms from the e-e cusp, the coefficients are
c those of (s**n/2)*log(s**2/2), rather than (s**n)*log(s**2)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use atom_mod
      implicit real*8(a-h,o-z)

!JT      include '../vmc/vmc.h'

!JT      include 'fit.h'

      common /dim/ ndim
      common /pars/ a00,a20,a21,eps_fock,c0000,c1110,c2000,
     &   xm1,xm2,xm12,xms,xma,Zfock

!JT      include '../vmc/force.h'

!JT      parameter (zero=0.d0,two=2.d0,four=4.d0
!JT     &,half=0.5d0)
      parameter(d1b4=0.25d0,d1b6=1.d0/6.d0,d1b12=1.d0/12.d0
     &,d1b24=1.d0/24.d0, rt22=2.82842712474619d0,rt2i=rt22/4.d0
     &,dln2=0.6931471805599453d0, pi=3.141592653589793d0
     &,const1=(1.d0-dln2)/12.d0,const2=-(pi-2.d0)/(6.d0*pi))

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /focsav/ a4sav,a5sav,a6sav,a7sav

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata

      dimension a(*),b(*),diff(*)

c     the shifting is done in the calling routine, so do not shift here
c     for e-e cusp, but do shift for e-n cusp and for making coefs of
c     s and r pos and neg respectively.
      ishft=0

c coordinates: r,s,t

c Warning: temporarily assume only one nucleus type

c f2e    o(s^2) from phi20(r12=0)
c f2elog o(s^2 log(s)) from phi21(r12=0)
c f2n    o(r^2) from phi20(r1=0)
      f2e=-d1b12*eguess
     &+ (const1-d1b4*znuc(iwctype(1)))*znuc(iwctype(1))
      f2elog=const2*znuc(iwctype(1))
      f2n=-d1b6*(eguess+znuc(iwctype(1))*(znuc(iwctype(1))-dln2))
     &- d1b24

c Focks's terms are evaluated using scaled variables. Scaled variables
c in phi10 originate a 2nd order in the real variables that has to be
c subtracted from phi20. This is done by redifining the coeffs a(4-7).
c If wf is a Pade', for small hyper-rad wf=P(a)(1-P(b)+..) and the 2nd
c order contribution from P(a)P(b) etc. has also to be subtracted out.
      if(ifock.gt.0) then
        if(isc.eq.2) then
          scalel=half*scalek
         elseif(isc.eq.4) then
          scalel=scalek
        endif
        const3=const2*znuc(iwctype(1))
        a(4)=a4sav+a(44)*(a(1)*(b(1)+scalel)-a4sav)
        a(5)=a5sav+a(44)*(a(2)*(b(2)+half*scalel)-a5sav)
        a(6)=a6sav+a(44)*(a(2)*half*scalel-a6sav)
        a(7)=a7sav+a(44)*(a(1)*b(2)+a(2)*b(1)-a7sav)

c Same as before for terms generated from phi21 contributing to phi31.
c The coeffs a(47-52) are defined to cancel this contribution from phi31.
        a(47)=const3*(-(a(1)*b(40)+a(40)*b(1))-scalel*a(40))
        a(48)=const3*(half*(a(2)*b(40)+a(40)*b(2))+scalel*a(40)/four)
        a(49)=const3*(-(a(2)*b(40)+a(40)*b(2)))
        a(50)=const3*(half*(a(1)*b(40)+a(40)*b(1)))
        a(51)=a(50)
        a(52)=a(48)+const3*scalel*a(40)/two
      endif

      ap1=a(1)
      ap2=a(7)+rt2i*a(35)
      ap3=a(14)
      ap4=a(24)
      a1=a(2)
      a2=a(5)+rt2i*a(36)+f2e*a(44)+half*a21
      a3=a(11)
      a4=a(21)
      al2=a(38)+f2elog*a(40)
      al3=2*a(48)+half*f2elog*znuc(iwctype(1))*a(45)
      alp3=a(39)+2*a(50)
      bp1=b(1)
      bp2=b(7)+rt2i*b(35)
      bp3=b(14)
      bp4=b(24)
      b1=b(2)
      b2=b(5)+rt2i*b(36)+f2e*b(44)
      b3=b(11)
      b4=b(21)
      bl2=b(38)+f2elog*b(40)
      bl3=2*b(48)+half*f2elog*znuc(iwctype(1))*b(45)
      blp3=b(39)+2*b(50)
      if(ncuspc.eq.13) then
        ap4=a(20)
        a4=zero
        bp4=b(20)
        b4=zero
      endif

      if(ifock.ge.1) then
        diff(ishft+1) = zero
       else
        if(isp.eq.1) then
          diff(ishft+1) = a(1)-half
         else
          diff(ishft+1) = a(1)-d1b4
        endif
      endif

      diff(ishft+2) =  ap1*b1-a1*bp1+ap2- 2*ap1*b1
      diff(ishft+3) =  ap1*b2-a2*bp1
     &+ ap2*b1-a1*bp2 + ap3- ap1*(2*b2+b1**2)
      diff(ishft+4) =  ap1*b3-a3*bp1
     &+ ap2*b2-a2*bp2 + ap3*b1-a1*bp3 + ap4-
     &2*ap1*(b1*b2 + b3)
      diff(ishft+5) =  ap1*b4-a4*bp1 +ap2*b3-a3*bp2
     &+ ap3*b2-a2*bp3 + ap4*b1-a1*bp4-
     &ap1*(2*(b4+b1*b3)+b2**2)
      diff(ishft+6) =  ap2*b4-a4*bp2 +ap3*b3-a3*bp3
     &+ ap4*b2-a2*bp4- 2*ap1*(b2*b3+b1*b4)
      diff(ishft+7) =  ap3*b4-a4*bp3 +ap4*b3-a3*bp4
     &- ap1*(2*b2*b4+b3**2)
      diff(ishft+8) =  ap4*b4-a4*bp4 -2*ap1*b3*b4
      diff(ishft+9) =  ap1*b4**2
      diff(ishft+10) = alp3 - ap1*bl2 - al2*bp1
      diff(ishft+11) = alp3*b1 + ap2*bl2 - 2*ap1*b1*bl2 - ap1*bl3
     &- a1*blp3 - al3*bp1 - al2*bp2
      diff(ishft+12) = alp3*b2 + ap3*bl2 - 2*ap1*b2*bl2 + ap2*bl3
     &- 2*ap1*b1*bl3 - a2*blp3 - al3*bp2 - al2*bp3
      diff(ishft+13) = alp3*b3 + ap4*bl2 - 2*ap1*b3*bl2 + ap3*bl3
     &- 2*ap1*b2*bl3 - a3*blp3 - al3*bp3 - al2*bp4
      diff(ishft+14) = alp3*b4 - 2*ap1*b4*bl2 + ap4*bl3 - 2*ap1*b3*bl3
     &- a4*blp3 - al3*bp4
      diff(ishft+15) =-2*ap1*b4*bl3
      diff(ishft+16) = alp3*bl2 - ap1*bl2**2 - al2*blp3
      diff(ishft+17) = alp3*bl3 - 2*ap1*bl2*bl3 - al3*blp3
      diff(ishft+18) =-ap1*bl3**2

      if(iprin.ge.1) then

        write(6,'(''a1,a2,a3,a4, ap1,ap2,ap3,ap4, al2,al3,alp3'',
     &  19f10.5)')  a1,a2,a3,a4, ap1,ap2,ap3,ap4, al2,al3,alp3
        write(6,'(''b1,b2,b3,b4, bp1,bp2,bp3,bp4, bl2,bl3,blp3'',
     &  19f10.5)')  b1,b2,b3,b4, bp1,bp2,bp3,bp4, bl2,bl3,blp3
        write(6,'(''ap2-ap1*b1,ap3-ap1*b2,ap4-ap1*b3,-ap1*b4'',
     &  19f10.5)')  ap2-ap1*b1,ap3-ap1*b2,ap4-ap1*b3,-ap1*b4

        write(6,'(''coefs of s   '',9f12.6)') ap1*b1-a1*bp1+ap2,
     &  2*ap1*b1
        write(6,'(''coefs of s**2'',9f12.6)') ap1*b2-a2*bp1
     &  + ap2*b1-a1*bp2 + ap3, ap1*(2*b2+b1**2)
        write(6,'(''coefs of s**3'',9f12.6)') ap1*b3-a3*bp1
     &  + ap2*b2-a2*bp2 + ap3*b1-a1*bp3 + ap4,
     &  2*ap1*(b1*b2 + b3)
        write(6,'(''coefs of s**4'',9f12.6)') ap1*b4-a4*bp1
     &  +ap2*b3-a3*bp2 +ap3*b2-a2*bp3 +
     &  ap4*b1-a1*bp4, ap1*(2*(b4+b1*b3)+b2**2)
        write(6,'(''coefs of s**5'',9f12.6)') ap2*b4-a4*bp2
     &  + ap3*b3-a3*bp3
     &  + ap4*b2-a2*bp4, 2*ap1*(b2*b3+b1*b4)
        write(6,'(''coefs of s**6'',9f12.6)') ap3*b4-a4*bp3
     &  +ap4*b3-a3*bp4, ap1*(2*b2*b4+b3**2)
        write(6,'(''coefs of s**7'',9f12.6)') ap4*b4-a4*bp4,
     &  2*ap1*b3*b4
        write(6,'(''coefs of s**8'',12x,9f12.6)') ap1*b4**2
        write(6,'(''coefs of s2lg'',9f12.6)')
     &  alp3 - ap1*bl2 - al2*bp1
        write(6,'(''coefs of s3lg'',9f12.6)')
     &  alp3*b1 + ap2*bl2 - 2*ap1*b1*bl2 - ap1*bl3
     &  - a1*blp3 - al3*bp1 - al2*bp2
        write(6,'(''coefs of s4lg'',9f12.6)')
     &  alp3*b2 + ap3*bl2 - 2*ap1*b2*bl2 + ap2*bl3
     &  - 2*ap1*b1*bl3 - a2*blp3 - al3*bp2 - al2*bp3
        write(6,'(''coefs of s5lg'',9f12.6)')
     &  alp3*b3 + ap4*bl2 - 2*ap1*b3*bl2 + ap3*bl3
     &  - 2*ap1*b2*bl3 - a3*blp3 - al3*bp3 - al2*bp4
        write(6,'(''coefs of s6lg'',9f12.6)')
     &  alp3*b4 - 2*ap1*b4*bl2 + ap4*bl3 - 2*ap1*b3*bl3
     &  - a4*blp3 - al3*bp4
        write(6,'(''coefs of s7lg'',9f12.6)') -2*ap1*b4*bl3
        write(6,'(''coefs o s4lg2'',9f12.6)')
     &  alp3*bl2 - ap1*bl2**2 - al2*blp3
        write(6,'(''coefs o s5lg2'',9f12.6)')
     &  alp3*bl3 - 2*ap1*bl2*bl3 - al3*blp3
        write(6,'(''coefs o s6lg2'',9f12.6)') -ap1*bl3**2
      endif

c     e-n cusp
c     note that ishft is one less that the number of cusps already
c     imposed because, the next cusp-condit has argument (ishft+2)
      ishft=17

      a1=0
      a2=f2n*a(44)
      a3=0
      a4=0
      al2=0
      al3=d1b6*f2elog*a(45)
      b1=0
      b2=f2n*b(44)
      b3=0
      b4=0
      bl2=0
      bl3=d1b6*f2elog*b(45)
      do 32 i=1,ndim
        a1=a1+a(i)
   32   b1=b1+b(i)
      do 33 i=4,9
        a2=a2+a(i)
   33   b2=b2+b(i)
      do 34 i=10,19
        a3=a3+a(i)
   34   b3=b3+b(i)
      do 35 i=20,34
        a4=a4+a(i)
   35   b4=b4+b(i)
      do 36 i=35,37
        a2=a2+a(i)
   36   b2=b2+b(i)
      do 37 i=53,60
        a2=a2+a(i)
   37   b2=b2+b(i)
      do 38 i=61,66
        a3=a3+a(i)
   38   b3=b3+b(i)
      do 39 i=67,69
        a4=a4+a(i)
   39   b4=b4+b(i)
      do 40 i=47,52
        al3=al3+2*a(i)
   40   bl3=bl3+2*b(i)
      ap1=a(2)-a(3)
      ap2=2*(a(5)-a(6))+a(7)-a(9) + a(36) - a(53)-3*a(54)-2*a(55)
     &-4*a(56)-6*a(57) -3*a(58)-5*a(59)-7*a(60) -a(61)-3*a(62)-5*a(63)
     &-2*a(64)-4*a(65)-6*a(66) -a(67)-3*a(68)-5*a(69)
      ap3=3*(a(11)-a(12)) + a(13)-a(15) + 2*(a(14)-a(16)) +a(17)-a(18)
     &-a(61)-3*a(62)-5*a(63)-2*a(64)-4*a(65)-6*a(66)
      ap4=4*(a(21)-a(22)) +3*(a(24)-a(26))
     &+2*(a(27)-a(28) +a(32)-a(33)) +a(23)-a(25) +a(30)-a(31)
      alp3=2*(3*a(48)+a(49)+2*a(50)-2*a(51)-a(52))
      bp1=b(2)-b(3)
      bp2=2*(b(5)-b(6))+b(7)-b(9) + b(36) - b(53)-3*b(54)-2*b(55)
     &-4*b(56)-6*b(57) -3*b(58)-5*b(59)-7*b(60) -b(61)-3*b(62)-5*b(63)
     &-2*b(64)-4*b(65)-6*b(66) -b(67)-3*b(68)-5*b(69)
      bp3=3*(b(11)-b(12)) + b(13)-b(15) + 2*(b(14)-b(16)) +b(17)-b(18)
     &-b(61)-3*b(62)-5*b(63)-2*b(64)-4*b(65)-6*b(66)
      bp4=4*(b(21)-b(22)) +3*(b(24)-b(26))
     &+2*(b(27)-b(28) +b(32)-b(33)) +b(23)-b(25) +b(30)-b(31)
      blp3=2*(3*b(48)+b(49)+2*b(50)-2*b(51)-b(52))

      diff(ishft+2) =  ap1*b1-a1*bp1+ap2- 2*ap1*b1
      diff(ishft+3) =  ap1*b2-a2*bp1
     &+ ap2*b1-a1*bp2 + ap3- ap1*(2*b2+b1**2)
      diff(ishft+4) =  ap1*b3-a3*bp1
     &+ ap2*b2-a2*bp2 + ap3*b1-a1*bp3 + ap4-
     &2*ap1*(b1*b2 + b3)
      diff(ishft+5) =  ap1*b4-a4*bp1 +ap2*b3-a3*bp2
     &+ ap3*b2-a2*bp3 + ap4*b1-a1*bp4-
     &ap1*(2*(b4+b1*b3)+b2**2)
      diff(ishft+6) =  ap2*b4-a4*bp2 +ap3*b3-a3*bp3
     &+ ap4*b2-a2*bp4- 2*ap1*(b2*b3+b1*b4)
      diff(ishft+7) =  ap3*b4-a4*bp3 +ap4*b3-a3*bp4
     &  - ap1*(2*b2*b4+b3**2)
      diff(ishft+8) =  ap4*b4-a4*bp4 -2*ap1*b3*b4
      diff(ishft+9) =  ap1*b4**2
      diff(ishft+10) = alp3 - ap1*bl2 - al2*bp1
      diff(ishft+11) = alp3*b1 + ap2*bl2 - 2*ap1*b1*bl2 - ap1*bl3
     &- a1*blp3 - al3*bp1 - al2*bp2
      diff(ishft+12) = alp3*b2 + ap3*bl2 - 2*ap1*b2*bl2 + ap2*bl3
     &- 2*ap1*b1*bl3 - a2*blp3 - al3*bp2 - al2*bp3
      diff(ishft+13) = alp3*b3 + ap4*bl2 - 2*ap1*b3*bl2 + ap3*bl3
     &- 2*ap1*b2*bl3 - a3*blp3 - al3*bp3 - al2*bp4
      diff(ishft+14) = alp3*b4 - 2*ap1*b4*bl2 + ap4*bl3 - 2*ap1*b3*bl3
     &- a4*blp3 - al3*bp4
      diff(ishft+15) =-2*ap1*b4*bl3
      diff(ishft+16) = alp3*bl2 - ap1*bl2**2 - al2*blp3
      diff(ishft+17) = alp3*bl3 - 2*ap1*bl2*bl3 - al3*blp3
      diff(ishft+18) =-ap1*bl3**2
      if(ifock.ge.1) then
        diff(ishft+19) =a(40)-1
        diff(ishft+20) =a(44)-1
        diff(ishft+21) =a(45)-1
      endif

      if(iprin.ge.1) then
        write(6,'(''a1,a2,a3,a4, ap1,ap2,ap3,ap4, al2,al3,alp3'',
     &  19f10.5)')  a1,a2,a3,a4, ap1,ap2,ap3,ap4, al2,al3,alp3
        write(6,'(''b1,b2,b3,b4, bp1,bp2,bp3,bp4, bl2,bl3,blp3'',
     &  19f10.5)')  b1,b2,b3,b4, bp1,bp2,bp3,bp4, bl2,bl3,blp3
        write(6,'(''ap2-ap1*b1,ap3-ap1*b2,ap4-ap1*b3,-ap1*b4'',
     &  19f10.5)')  ap2-ap1*b1,ap3-ap1*b2,ap4-ap1*b3,-ap1*b4

        write(6,'(''coefs of r   '',9f12.6)') ap1*b1-a1*bp1+ap2,
     &  2*ap1*b1
        write(6,'(''coefs of r**2'',9f12.6)') ap1*b2-a2*bp1
     &  + ap2*b1-a1*bp2 + ap3, ap1*(2*b2+b1**2)
        write(6,'(''coefs of r**3'',9f12.6)') ap1*b3-a3*bp1
     &  + ap2*b2-a2*bp2 + ap3*b1-a1*bp3 + ap4,
     &  2*ap1*(b1*b2 + b3)
        write(6,'(''coefs of r**4'',9f12.6)') ap1*b4-a4*bp1
     &  +ap2*b3-a3*bp2 +ap3*b2-a2*bp3 +
     &  ap4*b1-a1*bp4, ap1*(2*(b4+b1*b3)+b2**2)
        write(6,'(''coefs of r**5'',9f12.6)') ap2*b4-a4*bp2
     &  + ap3*b3-a3*bp3
     &  + ap4*b2-a2*bp4, 2*ap1*(b2*b3+b1*b4)
        write(6,'(''coefs of r**6'',9f12.6)') ap3*b4-a4*bp3
     &  +ap4*b3-a3*bp4, ap1*(2*b2*b4+b3**2)
        write(6,'(''coefs of r**7'',9f12.6)') ap4*b4-a4*bp4,
     &  2*ap1*b3*b4
        write(6,'(''coefs of r**8'',12x,9f12.6)') ap1*b4**2
        write(6,'(''coefs of r2lg'',9f12.6)')
     &  alp3 - ap1*bl2 - al2*bp1
        write(6,'(''coefs of r3lg'',9f12.6)')
     &  alp3*b1 + ap2*bl2 - 2*ap1*b1*bl2 - ap1*bl3
     &  - a1*blp3 - al3*bp1 - al2*bp2
        write(6,'(''coefs of r4lg'',9f12.6)')
     &  alp3*b2 + ap3*bl2 - 2*ap1*b2*bl2 + ap2*bl3
     &  - 2*ap1*b1*bl3 - a2*blp3 - al3*bp2 - al2*bp3
        write(6,'(''coefs of r5lg'',9f12.6)')
     &  alp3*b3 + ap4*bl2 - 2*ap1*b3*bl2 + ap3*bl3
     &  - 2*ap1*b2*bl3 - a3*blp3 - al3*bp3 - al2*bp4
        write(6,'(''coefs of r6lg'',9f12.6)')
     &  alp3*b4 - 2*ap1*b4*bl2 + ap4*bl3 - 2*ap1*b3*bl3
     &  - a4*blp3 - al3*bp4
        write(6,'(''coefs of r7lg'',9f12.6)') -2*ap1*b4*bl3
        write(6,'(''coefs o r4lg2'',9f12.6)')
     &  alp3*bl2 - ap1*bl2**2 - al2*blp3
        write(6,'(''coefs o r5lg2'',9f12.6)')
     &  alp3*bl3 - 2*ap1*bl2*bl3 - al3*blp3
        write(6,'(''coefs o r6lg2'',9f12.6)') -ap1*bl3**2
        if(ifock.ge.1) then
          write(6,'(''coefs o a40-1'',9f12.6)') a(40)-1
          write(6,'(''coefs o a44-1'',9f12.6)') a(44)-1
          write(6,'(''coefs o a45-1'',9f12.6)') a(45)-1
        endif
      endif

      return
      end

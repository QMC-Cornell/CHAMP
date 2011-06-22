      function psi(rij,ri,rj,it)
c Written by Cyrus Umrigar, modified by Claudia Filippi
c **Warning** This routine needs to be upgraded to check rshifts
c if we add in the capability to use numerical Laplacian for
c periodic systems.
c  NOTE:  this should work for at least a 1D periodic system now
c    I have added in the scaling for een jastrow (ACM, June 2011)
      use constants_mod 
      use dets_mod
      use contr2_mod
      use wfsec_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      use pars_mod
      use jaspar1_mod
      use jaspar2_mod
      use contrl_per_mod 
      implicit real*8(a-h,o-z)

      common /chck/ bot

      dimension uu(0:max(nord,nordb,nordc)),ss(0:max(nord,norda,nordc)),tt(0:max(nord,norda,nordc))

      dlogs4(x) = 2*dlog((one-dexp(-a1(41,is,iwf)*x))/a1(41,is,iwf))

      psi=0

      u=rij
      s=ri+rj
      t=dabs(ri-rj)

      call scale_dist(rij,u,4)
      call scale_dist(ri,rri,3)
      call scale_dist(rj,rrj,3)

      if(ijas.eq.1) then
        psi=sspinn*cjas1(iwf)*rij/(one+cjas2(iwf)*rij)

       elseif(ijas.eq.2) then
        s2=s*s
        s3=s2*s
        s4=s2*s2
        si=one/s
        si2=si*si
        si3=si2*si
        t2=t*t
        t4=t2*t2
        rr2=(s2+t2)/2

        u2=u*u
        u3=u2*u
        u4=u2*u2
        u5=u3*u2
        rr=dsqrt(rr2)

        if(ifock.gt.0) call psigm(u,rri,rrj,phi21,phi20,phi31,1)

        y21=rr2-u2

        if(a1(41,is,iwf).gt.0.01d0) then
          dlnrr2=dlogs4(rr)
         else
          dlnrr2=dlog(rr2)
        endif

        top=       a1(1,is,iwf)*u + a1(2,is,iwf)*s + a1(4,is,iwf)*u2
     &  + a1(5,is,iwf)*s2  + a1(6,is,iwf)*t2 + a1(7,is,iwf)*u*s
     &                + a1(10,is,iwf)*u3 + a1(11,is,iwf)*s3
     &  + a1(13,is,iwf)*u2*s + a1(14,is,iwf)*u*s2 + a1(16,is,iwf)*u*t2
     &                  + a1(18,is,iwf)*s*t2
     &  + a1(20,is,iwf)*u4 + a1(21,is,iwf)*s4 + a1(22,is,iwf)*t4
     &                  + a1(23,is,iwf)*u3*s
     &  + a1(24,is,iwf)*u*s3
     &  + a1(31,is,iwf)*u*s*t2 + a1(32,is,iwf)*u2*s2
     &                  + t2*(a1(33,is,iwf)*u2+a1(34,is,iwf)*s2)
     &  + (a1(35,is,iwf)*u + a1(36,is,iwf)*s)*rr + a1(37,is,iwf)*u3/rr
     &  + a1(40,is,iwf)*phi21+a1(44,is,iwf)*phi20+a1(45,is,iwf)*phi31
     &                   +a21*y21
     &  + ((a1(38,is,iwf)+a1(39,is,iwf)*u) * (rr2-u2)
     &  + a1(47,is,iwf)*u3+a1(48,is,iwf)*s3+a1(49,is,iwf)*u2*s
     &                   +a1(50,is,iwf)*u*s2
     &  + (a1(51,is,iwf)*u+a1(52,is,iwf)*s)*t2) * dlnrr2
     &  + (a1(53,is,iwf)*u3+a1(54,is,iwf)*u*t2)*si
     &  + (a1(55,is,iwf)*u4+a1(56,is,iwf)*u2*t2+a1(57,is,iwf)*t4)*si2
     &  + (a1(58,is,iwf)*u5+a1(59,is,iwf)*u3*t2+a1(60,is,iwf)*u*t4)*si3
     &  + (a1(61,is,iwf)*u4+a1(62,is,iwf)*u2*t2+a1(63,is,iwf)*t4)*si
     &  + (a1(64,is,iwf)*u5+a1(65,is,iwf)*u3*t2+a1(66,is,iwf)*u*t4)*si2
     &  + (a1(67,is,iwf)*u5+a1(68,is,iwf)*u3*t2+a1(69,is,iwf)*u*t4)*si

        bot = one +  a2(1,is,iwf)*u + a2(2,is,iwf)*s + a2(4,is,iwf)*u2
     &  + a2(5,is,iwf)*s2  + a2(6,is,iwf)*t2 + a2(7,is,iwf)*u*s
     &                + a2(10,is,iwf)*u3 + a2(11,is,iwf)*s3
     &  + a2(13,is,iwf)*u2*s + a2(14,is,iwf)*u*s2 + a2(16,is,iwf)*u*t2
     &                  + a2(18,is,iwf)*s*t2
     &  + a2(20,is,iwf)*u4 + a2(21,is,iwf)*s4 + a2(22,is,iwf)*t4
     &                    + a2(23,is,iwf)*u3*s
     &  + a2(24,is,iwf)*u*s3
     &  + a2(31,is,iwf)*u*s*t2 + a2(32,is,iwf)*u2*s2
     &                    +t2*(a2(33,is,iwf)*u2+a2(34,is,iwf)*s2)
     &  + (a2(35,is,iwf)*u + a2(36,is,iwf)*s)*rr + a2(37,is,iwf)*u3/rr
     &  + a2(40,is,iwf)*phi21+a2(44,is,iwf)*phi20+a2(45,is,iwf)*phi31
     &  + ((a2(38,is,iwf)+a2(39,is,iwf)*u) * (rr2-u2)
     &  + a2(47,is,iwf)*u3+a2(48,is,iwf)*s3+a2(49,is,iwf)*u2*s
     &                    +a2(50,is,iwf)*u*s2
     &  + (a2(51,is,iwf)*u+a2(52,is,iwf)*s)*t2) * dlnrr2
     &  + (a2(53,is,iwf)*u3+a2(54,is,iwf)*u*t2)*si
     &  + (a2(55,is,iwf)*u4+a2(56,is,iwf)*u2*t2+a2(57,is,iwf)*t4)*si2
     &  + (a2(58,is,iwf)*u5+a2(59,is,iwf)*u3*t2+a2(60,is,iwf)*u*t4)*si3
     &  + (a2(61,is,iwf)*u4+a2(62,is,iwf)*u2*t2+a2(63,is,iwf)*t4)*si
     &  + (a2(64,is,iwf)*u5+a2(65,is,iwf)*u3*t2+a2(66,is,iwf)*u*t4)*si2
     &  + (a2(67,is,iwf)*u5+a2(68,is,iwf)*u3*t2+a2(69,is,iwf)*u*t4)*si

        psi=sspin*top/bot

       elseif(ijas.eq.3) then

        if(nord.eq.0) return

        uu(0)=one
        ss(0)=one
        tt(0)=one

        uu(1)=u
        ss(1)=rri+rrj
        tt(1)=rri-rrj
        do 32 jp=2,nord
          uu(jp)=uu(jp-1)*uu(1)
          ss(jp)=ss(jp-1)*ss(1)
   32     tt(jp)=tt(jp-1)*tt(1)

        ll=0
        do 34 jp=1,nord
          do 34 ju=jp,0,-1
            pc=uu(ju)
            jsx=jp-ju
            do 34 js=jsx,0,-1
              ll=ll+1
              jt=jsx-js

              if(mod(jt,2).ne.0.and.nup.eq.ndn) then
                c(ll,it,iwf)=zero
              endif

              p=pc*ss(js)*tt(jt)
              psi=psi+c(ll,it,iwf)*p

   34         continue

       if(ifock.gt.0) then
         r2=(ss(2)+tt(2))/two
         r =dsqrt(r2)
         rin=one/r
         rlog=two*dlog(r/(one+r))
         y21=r2-uu(2)
         call psigm(uu(1),rri,rrj,phi21,phi20,phi31,it)

         psi=psi+fck(2,it,iwf)*phi20
     &   +(fck(4,it,iwf)*uu(1)+fck(5,it,iwf)*ss(1))*r
     &   +(fck(6,it,iwf)*uu(3)+fck(7,it,iwf)*ss(3)
     &   +fck(8,it,iwf)*uu(2)*ss(1)
     &   +fck(9,it,iwf)*uu(1)*ss(2))*rin+a21*y21

         psi=psi+fck(1,it,iwf)*phi21+fck(3,it,iwf)*phi31
     &   +(fck(10,it,iwf)*uu(3)+fck(11,it,iwf)*ss(3)
     &   +fck(12,it,iwf)*uu(2)*ss(1)
     &   +fck(13,it,iwf)*uu(1)*ss(2)+(fck(14,it,iwf)*uu(1)
     &   +fck(15,it,iwf)*ss(1))*tt(2))*rlog
       endif

       elseif(ijas.ge.4.and.ijas.le.6) then

c       een Jastrow:
        if(nordc.le.1) return

        if(ri.gt.cutjas_en .or. rj.gt.cutjas_en) return
c       do 37 k=1,ndim
c  37     if(abs(rshift(k,i,ic)-rshift(k,j,ic)).gt.eps) return

        if(ijas.eq.4.or.ijas.eq.5) then
          call switch_scale(u,4)
          call switch_scale(rri,3)
          call switch_scale(rrj,3)
        endif
c     write(6,'(''rij,u in een'',2f12.9)') rij,u
c     write(6,'(''ri,rri in een'',2f12.9)') ri,rri

c       Extra scaling function for een jastrow in periodic case (ACM)
        psic = 0.d0
        fscale = 1.0d0
        if(iperiodic.ne.0) then
          call f_een_cuts_nd(cutjas_en, ri, rj, fscale)
        endif
        
        uu(0)=one
        ss(0)=2
        tt(0)=one
        do 40 jp=1,nordc
          uu(jp)=u**jp
          ss(jp)=rri**jp+rrj**jp
   40     tt(jp)=(rri*rrj)**jp

        ll=0
        do 50 n=2,nordc
          do 50 k=n-1,0,-1
            if(k.eq.0) then
              l_hi=n-k-2
             else
              l_hi=n-k
            endif
            do 50 l=l_hi,0,-1
              m=(n-k-l)/2
              if(2*m.eq.n-k-l) then
                ll=ll+1
                psic=psic+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
              endif
c     write(6,'(''rij,ri,rj'',9f10.5)') rij,ri,rj,u,rri,rrj
   50   continue
      psi = psi + fscale*psic


      endif

      return
      end

c-----------------------------------------------------------------------
      function psia(ri,it)

      use contr2_mod
      use wfsec_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      implicit real*8(a-h,o-z)



!JT     parameter(zero=0.d0,one=1.d0)



      psia=zero
      if(ijas.lt.3.or.ijas.gt.6) return

      if(ri.gt.cutjas_en) return

      call scale_dist(ri,rri,1)
c     write(6,'(''ri,rri in en'',2f9.5)') ri,rri

      if(ijas.eq.3) then
        psia=a(1,iwf)*rri/(one+a(2,iwf)*rri)
       elseif(ijas.eq.4.or.ijas.eq.5) then
        psia=a4(1,it,iwf)*rri/(one+a4(2,it,iwf)*rri)-asymp_jasa(it,iwf)
        do 10 i=2,norda
   10     psia=psia+a4(i+1,it,iwf)*rri**i
       elseif(ijas.eq.6) then
        psia=a4(1,it,iwf)*rri/(one+a4(2,it,iwf)*(1-rri))
        do 20 i=2,norda
   20     psia=psia+a4(i+1,it,iwf)*rri**i
      endif

      return
      end

c-----------------------------------------------------------------------

      function psib(rij,isb,ipar)

      use contr2_mod
      use wfsec_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      implicit real*8(a-h,o-z)

!JT      parameter(zero=0.d0,one=1.d0)



      psib=zero
      if(ijas.lt.3.or.ijas.gt.6) return

      if(rij.gt.cutjas_ee) return

      u=rij
      call scale_dist(rij,u,2)
c     write(6,'(''rij,u in ee'',2f9.5)') rij,u

      if(ijas.eq.4) then
        psib=sspinn*b(1,isb,iwf)*u/(one+b(2,isb,iwf)*u)-asymp_jasb(ipar+1,iwf)
        do 10 i=2,nordb
   10     psib=psib+b(i+1,isb,iwf)*u**i
       elseif(ijas.eq.5) then
        psib=b(1,isb,iwf)*u/(one+b(2,isb,iwf)*u)-asymp_jasb(ipar+1,iwf)
        do 20 i=2,nordb
   20     psib=psib+b(i+1,isb,iwf)*u**i
        psib=sspinn*psib
       elseif(ijas.eq.6) then
        psib=b(1,isb,iwf)*u/(one+b(2,isb,iwf)*(1-u))
        do 30 i=2,nordb
   30     psib=psib+b(i+1,isb,iwf)*u**i
        psib=sspinn*psib
      endif

      return
      end

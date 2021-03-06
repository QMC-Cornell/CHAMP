      function deriv_psinl(u,dkij,rshifti,rshiftj,ri,rj,rri,rrj,dki,dkj,gn,gns,it)
      use control_mod
! Written by Claudia Filippi, modified by Cyrus Umrigar
! minor modification by A.D.Guclu to add analytical scalek opt.
      use constants_mod
      use dets_mod
      use optim_mod
      use dim_mod
      use contr2_mod
      use wfsec_mod
      use contrl_per_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      use pars_mod
      use jaspar1_mod
      use jaspar2_mod
      use vardep_mod
      use cuspmat_mod
      use cuspmat4_mod
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-12)
!!!   added WAS
      common /jas_c_cut/ cutjasc,icutjasc
!!!
      dimension rshifti(3),rshiftj(3),gn(*)
      dimension uu(0:max(nord,nordb,nordc)),ss(0:max(nord,norda,nordc)),tt(0:max(nord,norda,nordc)),rrri(0:max(nord,norda,nordc)) &
     &,rrrj(0:max(nord,norda,nordc))

      dlogs4(x) = 2*dlog((one-dexp(-a1(41,is,iwf)*x))/a1(41,is,iwf))

      if(ijas.eq.1) then
!       deriv_psinl=sspinn*cjas1(iwf)*rij/(one+cjas2(iwf)*rij)
        deriv_psinl=sspinn*cjas1(iwf)*u/(one+cjas2(iwf)*u)

       elseif(ijas.eq.2) then
        s=rri+rrj
        t=dabs(rri-rrj)

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

        top=       a1(1,is,iwf)*u + a1(2,is,iwf)*s + a1(4,is,iwf)*u2 &
     &  + a1(5,is,iwf)*s2  + a1(6,is,iwf)*t2 + a1(7,is,iwf)*u*s &
     &                + a1(10,is,iwf)*u3 + a1(11,is,iwf)*s3 &
     &  + a1(13,is,iwf)*u2*s + a1(14,is,iwf)*u*s2 + a1(16,is,iwf)*u*t2 &
     &                  + a1(18,is,iwf)*s*t2 &
     &  + a1(20,is,iwf)*u4 + a1(21,is,iwf)*s4 + a1(22,is,iwf)*t4 &
     &                  + a1(23,is,iwf)*u3*s &
     &  + a1(24,is,iwf)*u*s3 &
     &  + a1(31,is,iwf)*u*s*t2 + a1(32,is,iwf)*u2*s2 &
     &                  + t2*(a1(33,is,iwf)*u2+a1(34,is,iwf)*s2) &
     &  + (a1(35,is,iwf)*u + a1(36,is,iwf)*s)*rr + a1(37,is,iwf)*u3/rr &
     &  + a1(40,is,iwf)*phi21+a1(44,is,iwf)*phi20+a1(45,is,iwf)*phi31 &
     &                   +a21*y21 &
     &  + ((a1(38,is,iwf)+a1(39,is,iwf)*u) * (rr2-u2) &
     &  + a1(47,is,iwf)*u3+a1(48,is,iwf)*s3+a1(49,is,iwf)*u2*s &
     &                   +a1(50,is,iwf)*u*s2 &
     &  + (a1(51,is,iwf)*u+a1(52,is,iwf)*s)*t2) * dlnrr2 &
     &  + (a1(53,is,iwf)*u3+a1(54,is,iwf)*u*t2)*si &
     &  + (a1(55,is,iwf)*u4+a1(56,is,iwf)*u2*t2+a1(57,is,iwf)*t4)*si2 &
     &  + (a1(58,is,iwf)*u5+a1(59,is,iwf)*u3*t2+a1(60,is,iwf)*u*t4)*si3 &
     &  + (a1(61,is,iwf)*u4+a1(62,is,iwf)*u2*t2+a1(63,is,iwf)*t4)*si &
     &  + (a1(64,is,iwf)*u5+a1(65,is,iwf)*u3*t2+a1(66,is,iwf)*u*t4)*si2 &
     &  + (a1(67,is,iwf)*u5+a1(68,is,iwf)*u3*t2+a1(69,is,iwf)*u*t4)*si

        bot = one +  a2(1,is,iwf)*u + a2(2,is,iwf)*s + a2(4,is,iwf)*u2 &
     &  + a2(5,is,iwf)*s2  + a2(6,is,iwf)*t2 + a2(7,is,iwf)*u*s &
     &                + a2(10,is,iwf)*u3 + a2(11,is,iwf)*s3 &
     &  + a2(13,is,iwf)*u2*s + a2(14,is,iwf)*u*s2 + a2(16,is,iwf)*u*t2 &
     &                  + a2(18,is,iwf)*s*t2 &
     &  + a2(20,is,iwf)*u4 + a2(21,is,iwf)*s4 + a2(22,is,iwf)*t4 &
     &                    + a2(23,is,iwf)*u3*s &
     &  + a2(24,is,iwf)*u*s3 &
     &  + a2(31,is,iwf)*u*s*t2 + a2(32,is,iwf)*u2*s2 &
     &                    +t2*(a2(33,is,iwf)*u2+a2(34,is,iwf)*s2) &
     &  + (a2(35,is,iwf)*u + a2(36,is,iwf)*s)*rr + a2(37,is,iwf)*u3/rr &
     &  + a2(40,is,iwf)*phi21+a2(44,is,iwf)*phi20+a2(45,is,iwf)*phi31 &
     &  + ((a2(38,is,iwf)+a2(39,is,iwf)*u) * (rr2-u2) &
     &  + a2(47,is,iwf)*u3+a2(48,is,iwf)*s3+a2(49,is,iwf)*u2*s &
     &                    +a2(50,is,iwf)*u*s2 &
     &  + (a2(51,is,iwf)*u+a2(52,is,iwf)*s)*t2) * dlnrr2 &
     &  + (a2(53,is,iwf)*u3+a2(54,is,iwf)*u*t2)*si &
     &  + (a2(55,is,iwf)*u4+a2(56,is,iwf)*u2*t2+a2(57,is,iwf)*t4)*si2 &
     &  + (a2(58,is,iwf)*u5+a2(59,is,iwf)*u3*t2+a2(60,is,iwf)*u*t4)*si3 &
     &  + (a2(61,is,iwf)*u4+a2(62,is,iwf)*u2*t2+a2(63,is,iwf)*t4)*si &
     &  + (a2(64,is,iwf)*u5+a2(65,is,iwf)*u3*t2+a2(66,is,iwf)*u*t4)*si2 &
     &  + (a2(67,is,iwf)*u5+a2(68,is,iwf)*u3*t2+a2(69,is,iwf)*u*t4)*si

        deriv_psinl=sspin*top/bot

       elseif(ijas.eq.3) then

        deriv_psinl=zero
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
        jj=1
        jparm=1
        do 34 jp=1,nord
          do 34 ju=jp,0,-1
            pc=uu(ju)
            jsx=jp-ju
            do 34 js=jsx,0,-1
              ll=ll+1
              jt=jsx-js

              if(mod(jt,2).ne.0.and.nup.eq.ndn) then
                c(ll,it,iwf)=zero
               else
                p=pc*ss(js)*tt(jt)
                deriv_psinl=deriv_psinl+c(ll,it,iwf)*p
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

                if(ideriv.eq.1) then
                  do 33 id=1,nvdepend(jj,it)
                    iparm=iwdepend(jj,id,it)
   33               gn(iparm)=gn(iparm)+cdep(jj,id,it)*p
                  jj=jj+1
                 elseif(ideriv.eq.2) then
                  gn(jparm)=gn(jparm)+p
                  jparm=jparm+1
                endif
              endif
   34   continue

        if(ifock.gt.0) then
          r2=half*(ss(2)+tt(2))
          r =dsqrt(r2)
          rin=one/r
          rlog=two*dlog(r/(one+r))
          y21=r2-uu(2)
          call psigm(uu(1),rri,rrj,phi21,phi20,phi31,it)

          deriv_psinl=deriv_psinl+fck(2,it,iwf)*phi20 &
     &    +(fck(4,it,iwf)*uu(1)+fck(5,it,iwf)*ss(1))*r &
     &    +(fck(6,it,iwf)*uu(3)+fck(7,it,iwf)*ss(3) &
     &    +fck(8,it,iwf)*uu(2)*ss(1) &
     &    +fck(9,it,iwf)*uu(1)*ss(2))*rin+a21*y21

          deriv_psinl=deriv_psinl+fck(1,it,iwf)*phi21 &
     &    +fck(3,it,iwf)*phi31 &
     &    +(fck(10,it,iwf)*uu(3)+fck(11,it,iwf)*ss(3) &
     &    +fck(12,it,iwf)*uu(2)*ss(1) &
     &    +fck(13,it,iwf)*uu(1)*ss(2)+(fck(14,it,iwf)*uu(1) &
     &    +fck(15,it,iwf)*ss(1))*tt(2))*rlog
        endif

      elseif(ijas.ge.4.and.ijas.le.6) then

        deriv_psinl=0
        fu=0
        fi=0
        fj=0
        if(nordc.eq.0) return

        if(rri.eq.asymp_r_en(iwf) .or. rrj.eq.asymp_r_en(iwf)) return
        do 37 k=1,ndim
   37     if(abs(rshifti(k)-rshiftj(k)).gt.eps) return

        uu(1)=u
        rrri(1)=rri
        rrrj(1)=rrj
        call switch_scale(uu(1),4)
        call switch_scale(rrri(1),3)
        call switch_scale(rrrj(1),3)

!!!!  WAS
        if(icutjasc.gt.0 .or. iperiodic.ne.0) then
           call f_een_cuts_nd (cutjas_en, ri, rj, fcut)
        endif
!!!

        uu(0)=1
        ss(0)=2
        tt(0)=1
        rrri(0)=1
        rrrj(0)=1
        do 40 jp=1,nordc
          rrri(jp)=rrri(jp-1)*rrri(1)
          rrrj(jp)=rrrj(jp-1)*rrrj(1)
          uu(jp)=uu(jp-1)*uu(1)
          ss(jp)=rrri(jp)+rrrj(jp)
   40     tt(jp)=rrri(jp)*rrrj(jp)

        ll=0
        jj=1
        jparm=1
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
                p=uu(k)*ss(l)*tt(m)

!!WAS
                if(icutjasc.gt.0 .or. iperiodic.ne.0) then
                   p = p*fcut
                endif
!!!
                deriv_psinl=deriv_psinl+c(ll,it,iwf)*p

                if(nparms.eq.1) then
                  pu=k*uu(k-1)*ss(l)*tt(m)
                  pi=uu(k)*((l+m)*rrri(l+m-1)*rrrj(m)+m*rrri(m-1)*rrrj(l+m))
                  pj=uu(k)*((l+m)*rrrj(l+m-1)*rrri(m)+m*rrrj(m-1)*rrri(l+m))
                  fu=fu+c(ll,it,iwf)*pu
                  fi=fi+c(ll,it,iwf)*pi
                  fj=fj+c(ll,it,iwf)*pj
                endif

!               ideriv=0
!               if(ll.eq.iwc4(jj)) then
!                 if(nvdepend(jj,it).gt.0) then
!                   ideriv=1
!                  else
!                   jj=jj+1
!                 endif
!                elseif(ll.eq.iwjasc(jparm,it)) then
!                 ideriv=2
!               endif

                if(nparmc(it).gt.0) then
                  ideriv=0
                  if(ll.eq.iwjasc(jparm,it)) then
                    ideriv=2
                   else
                    do 31 id=1,2*(nordc-1)
                      if(ll.eq.iwc4(id)) then
                        jj=id
                        if(nvdepend(jj,it).gt.0) ideriv=1
                      endif
   31               continue
                  endif

                  if(ideriv.eq.1) then
                    do 43 id=1,nvdepend(jj,it)
                      iparm=iwdepend(jj,id,it)
   43                 gn(iparm)=gn(iparm)+cdep(jj,id,it)*p
!                   jj=jj+1
                   elseif(ideriv.eq.2) then
                    gn(jparm)=gn(jparm)+p
                    jparm=jparm+1
                  endif
                endif
              endif
   50   continue

        if(nparms.eq.1) then
          dkki=dki
          dkkj=dkj
          dkkij=dkij
! dd1,dd2,dk2,dr,dr2 are unused variables with the 0 option
          call switch_dscale(uu(1),dd1,dd2,dkkij,dk2,dr,dr2,4,0)
          call switch_dscale(rrri(1),dd1,dd2,dkki,dk2,dr,dr2,3,0)
          call switch_dscale(rrrj(1),dd1,dd2,dkkj,dk2,dr,dr2,3,0)
          gns=gns+fu*dkkij + fi*dkki + fj*dkkj
        endif

      endif

      return
      end

!-----------------------------------------------------------------------
      function deriv_psianl(rri,dk,gn,gns,it)
! Returns the e-n Jastrow exponent as the function value and
! the derivs. wrt the parameters in gn(*).
      use optim_mod
      use contr2_mod
      use wfsec_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      implicit real*8(a-h,o-z)

!JT   parameter(one=1.d0)

      dimension gn(*)

! Note: This routine is only called with iwf=1, but parts of it are
! written for general iwf, whereas others (asymp_r) assume iwf=1.

      if(rri.eq.asymp_r_en(iwf)) then
        deriv_psianl=0
        return
      endif

      if(ijas.eq.3) then
        deriv_psianl=a(1,iwf)*rri/(one+a(2,iwf)*rri)

        do 10 jparm=1,nparma(1)
          if(iwjasa(jparm,1).eq.1) then
            top=rri
            bot=one+a(2,iwf)*rri
           elseif(iwjasa(jparm,1).eq.2) then
            top=-a(1,iwf)*rri*rri
            bot=one+a(2,iwf)*rri
            bot=bot*bot
          endif
          gn(jparm)=gn(jparm)+top/bot
   10   continue

       elseif(ijas.ge.4.and.ijas.le.6) then

        bot=one+a4(2,it,iwf)*rri
        deriv_psianl=a4(1,it,iwf)*rri/bot-asymp_jasa(it,iwf)
        do 20 i=2,norda
   20     deriv_psianl=deriv_psianl+a4(i+1,it,iwf)*rri**i

        if(nparms.eq.1) then
          fenu=a4(1,it,iwf)/(bot*bot)
          do 25 i=2,nordb
   25       fenu=fenu+a4(i+1,it,iwf)*i*rri**(i-1)
          gns=gns+fenu*dk-dasymp_jasa(it,iwf)*dasymp_r_en(iwf)
        endif

        do 30 jparm=1,nparma(it)
            if(iwjasa(jparm,it).eq.1) then
              top=rri
              bot=one+a4(2,it,iwf)*rri
              gen=top/bot-asymp_r_en(iwf)/(1+a4(2,it,iwf)*asymp_r_en(iwf))
             elseif(iwjasa(jparm,it).eq.2) then
              top=-a4(1,it,iwf)*rri*rri
              bot=one+a4(2,it,iwf)*rri
              bot=bot*bot
              gen=top/bot+a4(1,it,iwf)*asymp_r_en(iwf)**2/(1+a4(2,it,iwf)*asymp_r_en(iwf))**2
             else
              iord=iwjasa(jparm,it)-1
              gen=rri**iord-asymp_r_en(iwf)**iord
            endif
            gn(jparm)=gn(jparm)+gen
   30    continue
      endif

      return
      end

!-----------------------------------------------------------------------
      function deriv_psibnl(u,dk,gn,gns,isb,ipar)
! Returns the e-e Jastrow exponent as the function value and
! the derivs. wrt the parameters in gn(*).

      use optim_mod
      use contr2_mod
      use wfsec_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      implicit real*8(a-h,o-z)

!JT   parameter(one=1.d0)

      dimension gn(*)

! Note: This routine is only called with iwf=1, but parts of it are
! written for general iwf, whereas others (asymp_r) assume iwf=1.

! Not updated for ijas=5,6 because we will probably stay with ijas=4
! If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) stop 'ijas >= 5 not implemented in psibnl'

      if(u.eq.asymp_r_ee(iwf)) then
        deriv_psibnl=0
        return
      endif

      top=sspinn*b(1,isb,iwf)
      bot=one+b(2,isb,iwf)*u
      fee=top*u/bot
      deriv_psibnl=fee-asymp_jasb(ipar+1,iwf)
      if(ijas.ge.4.and.ijas.le.6) then
        do 10 i=2,nordb
   10     deriv_psibnl=deriv_psibnl+b(i+1,isb,iwf)*u**i
      endif
      if(nparms.eq.1) then
        feeu=top/(bot*bot)
        do 15 i=2,nordb
   15     feeu=feeu+b(i+1,isb,iwf)*i*u**(i-1)
        gns=feeu*dk-dasymp_jasb(ipar+1,iwf)*dasymp_r_ee(iwf)
      endif
      do 20 jparm=1,nparmb(isb)
        if(iwjasb(jparm,isb).eq.1) then
          top=u
          bot=one+b(2,isb,iwf)*u
          gee=sspinn*(top/bot-asymp_r_ee(iwf)/(1+b(2,isb,iwf)*asymp_r_ee(iwf)))
         elseif(iwjasb(jparm,isb).eq.2) then
          top=-b(1,isb,iwf)*u*u
          bot=one+b(2,isb,iwf)*u
          bot=bot*bot
          gee=sspinn*(top/bot+b(1,isb,iwf)*asymp_r_ee(iwf)**2/(1+b(2,isb,iwf)*asymp_r_ee(iwf))**2)
         else
          iord=iwjasb(jparm,isb)-1
          gee=u**iord-asymp_r_ee(iwf)**iord
        endif
        gn(jparm)=gn(jparm)+gee
  20  continue
      return
      end

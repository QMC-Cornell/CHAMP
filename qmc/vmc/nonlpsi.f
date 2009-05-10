c     function psinl(u,rshifti,rshiftj,rri,rrj,it)
!WAS
      function psinl(u,rshifti,rshiftj,ri, rj, rri,rrj,it)
      use control_mod
!
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use dets_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'force.h'

!JT      parameter (one=1.d0,two=2.d0,half=0.5d0,eps=1.d-12)
      parameter (eps=1.d-12)

!!!   added WAS
      common /jas_c_cut/ cutjasc,icutjasc
      common /contrl_per/ iperiodic,ibasis
!!!
      common /dim/ ndim
      common /pars/ a00,a20,a21,eps_fock,c0000,c1110,c2000,
     &   xm1,xm2,xm12,xms,xma,Zfock
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar1/ cjas1(MWF),cjas2(MWF)
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /chck/ bot

      dimension uu(0:MORDJ),ss(0:MORDJ),tt(0:MORDJ),rshifti(3),rshiftj(3)

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) stop 'ijas >= 5 not implemented in psinl'

      if(ijas.eq.1) then

        psinl=sspinn*cjas1(iwf)*u/(1+cjas2(iwf)*u)

       elseif(ijas.eq.3) then

        psinl=0.d0
        if(nord.eq.0) return

        uu(0)=1.d0
        ss(0)=1.d0
        tt(0)=1.d0

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

              if(mod(jt,2).eq.0.or.nup.ne.ndn) then
                p=pc*ss(js)*tt(jt)
                psinl=psinl+c(ll,it,iwf)*p
              endif

   34   continue

        if(ifock.gt.0) then
          r2=half*(ss(2)+tt(2))
          r =dsqrt(r2)
          rin=one/r
          rlog=two*dlog(r/(one+r))
          y21=r2-uu(2)
          call psigm(uu(1),rri,rrj,phi21,phi20,phi31,it)

          psinl=psinl+fck(2,it,iwf)*phi20
     &    +(fck(4,it,iwf)*uu(1)+fck(5,it,iwf)*ss(1))*r
     &    +(fck(6,it,iwf)*uu(3)+fck(7,it,iwf)*ss(3)
     &    +fck(8,it,iwf)*uu(2)*ss(1)
     &    +fck(9,it,iwf)*uu(1)*ss(2))*rin+a21*y21

          psinl=psinl+fck(1,it,iwf)*phi21+fck(3,it,iwf)*phi31
     &    +(fck(10,it,iwf)*uu(3)+fck(11,it,iwf)*ss(3)
     &    +fck(12,it,iwf)*uu(2)*ss(1)
     &    +fck(13,it,iwf)*uu(1)*ss(2)+(fck(14,it,iwf)*uu(1)
     &    +fck(15,it,iwf)*ss(1))*tt(2))*rlog
        endif

       elseif(ijas.ge.4.and.ijas.le.6) then

        psinl=0.d0
        if(nordc.le.1) return

        if(rri.eq.asymp_r_en(iwf) .or. rrj.eq.asymp_r_en(iwf)) return
        do 37 k=1,ndim
   37     if(abs(rshifti(k)-rshiftj(k)).gt.eps) return

        uuu=u
        rrri=rri
        rrrj=rrj
        call switch_scale(uuu,4)
        call switch_scale(rrri,3)
        call switch_scale(rrrj,3)

!!!!  WAS
        if(icutjasc.gt.0 .or. iperiodic.ne.0) then
           call f_een_cuts_nd (cutjas_en, ri, rj, fcut)
        endif
!!!
        uu(0)=one
        ss(0)=two
        tt(0)=one
        do 40 jp=1,nordc
          uu(jp)=uuu**jp
          ss(jp)=rrri**jp+rrrj**jp
   40     tt(jp)=(rrri*rrrj)**jp

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
                p=uu(k)*ss(l)*tt(m)
c     write(6,'(''n,k,l,p,fcut='',3i4,9d12.4)') n,k,l,p,fcut
!!WAS
                if(icutjasc.gt.0 .or. iperiodic.ne.0) then
                   p=p*fcut
                endif
!!!
                psinl=psinl+c(ll,it,iwf)*p
              endif
   50   continue

      endif

      return
      end

c-----------------------------------------------------------------------

      function psianl(rri,it)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) stop 'ijas >= 5 not implemented in psianl'

      psianl=0.d0
      if(ijas.le.2) return
      if(rri.eq.asymp_r_en(iwf)) return

      if(ijas.eq.3) then
        psianl=a(1,iwf)*rri/(1.d0+a(2,iwf)*rri)
       elseif(ijas.ge.4.and.ijas.le.6) then
        psianl=a4(1,it,iwf)*rri/(1.d0+a4(2,it,iwf)*rri)-asymp_jasa(it,iwf)
        do 10 i=2,norda
   10     psianl=psianl+a4(i+1,it,iwf)*rri**i
      endif

      return
      end

c-----------------------------------------------------------------------
      function psibnl(u,isb,ipar)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

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
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) stop 'ijas >= 5 not implemented in psibnl'

      psibnl=0.d0
      if(ijas.le.2) return
      if(u.eq.asymp_r_ee(iwf)) return

      fee=b(1,isb,iwf)*u/(1+b(2,isb,iwf)*u)
      psibnl=sspinn*fee

      if(ijas.eq.4) then
        psibnl=sspinn*b(1,isb,iwf)*u/(1+b(2,isb,iwf)*u)-asymp_jasb(ipar+1,iwf)
        do 10 i=2,nordb
   10     psibnl=psibnl+b(i+1,isb,iwf)*u**i
       elseif(ijas.eq.5) then
        psibnl=b(1,isb,iwf)*u/(1+b(2,isb,iwf)*u)
        do 20 i=2,nordb
   20     psibnl=psibnl+b(i+1,isb,iwf)*u**i
        psibnl=sspinn*psibnl
       elseif(ijas.eq.6) then
        psibnl=b(1,isb,iwf)*u/(1+b(2,isb,iwf)*(1-u))
        do 30 i=2,nordb
   30     psibnl=psibnl+b(i+1,isb,iwf)*u**i
        psibnl=sspinn*psibnl
      endif

      return
      end

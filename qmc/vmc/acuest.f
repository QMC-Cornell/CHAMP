      subroutine acuest
c Written by Cyrus Umrigar, modified by Claudia Filippi
c routine to accumulate estimators for energy etc.
      use all_tools_mod
      use constants_mod
      use control_mod
      use montecarlo_mod
      use variables_mod
      use mpi_mod
      use atom_mod
      use config_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use forcepar_mod
!      use doefp_mod
      use pseudo_mod
      use contrl_per_mod
      use delocc_mod
      use contr3_mod
      use div_v_mod
      use distance_mod
      use qua_mod
      use forcest_mod
      use denupdn_mod
      use stepv_mod
      use pairden_mod
      use fourier_mod
      use forcewt_mod
      use est2cm_mod
      use estsig_mod
      use estcum_mod
      use estsum_mod
      implicit real*8(a-h,o-z)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot

      dimension xstrech(3,nelec)

c statement function for error calculation
c     err(x,x2)=dsqrt(dabs(x2/iblk-(x/iblk)**2)/iblk)
      err(x,x2,i)=dsqrt(abs(x2/wcum(i)-(x/wcum(i))**2)/iblk)

c xsum = sum of values of x from metrop
c xnow = average of values of x from metrop
c xcum = accumulated sums of xnow
c xcm2 = accumulated sums of xnow**2
c xave = current average value of x
c xerr = current error of x

      if(index(mode,'mpi').ne.0) call acuest_mpi

c Warning: At present we have nblk blocks of size nstep*nproc each.
c It would be a bit better to have nblk*nproc blocks of size nstep each,
c which just requires a bit of reorganizing.
      iblk=iblk+1
      it=iblk*nstep

c write out header first time
      if(iblk.eq.1) then
        if(ndim.eq.2) then
          write(6,'(t5,''enow'',t15,''eave'',t21,''(eerr )''
     &    ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr''
     &    ,t66,''tjfave'',t72,''(tjferr'',t83,''emave'',t89,''(emave)''
     &    ,t100,''fave'',t106,''(ferr)''
     &    ,t117,''accave'',t128,''iter'')')
         elseif(nloc.eq.-3) then !MS Jellium
          write(6,'(t9,''enow'',t21,''eave'',t25,''(eerr )''
     &    ,t39,''peave'',t44,''(peerr)'',t57,''tpbave'',t63,''(tpberr''
     &    ,t76,''tjfave'',t82,''(tjferr'',t93,''fave'',t99,''(ferr)''
     &    ,t110,''accave'',t122,''iter'',t131,''sigma'')')
         else
          write(6,'(t5,''enow'',t15,''eave'',t21,''(eerr )''
     &    ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr''
     &    ,t66,''tjfave'',t72,''(tjferr'',t83,''fave'',t89,''(ferr)''
     &    ,t100,''accave'',t111,''iter'',t120,''sigma'')')
        endif
      endif

c write out current values of averages

      wsum(1)=nstep*nproc
      walker_weights_sum_block  = nstep_total
      walker_weights_sum  = nstep_total*block_iterations_nb
      call object_modified_by_index (walker_weights_sum_block_index)
      call object_modified_by_index (walker_weights_sum_index)

      do 10 ifr=1,nforce
c       wnow=wsum(ifr)/nstep
        enow=esum(ifr)/wsum(ifr)
        wcum(ifr)=wcum(ifr)+wsum(ifr)
        ecum(ifr)=ecum(ifr)+esum(ifr)
        ecm2(ifr)=ecm2(ifr)+esum(ifr)*enow
        eave=ecum(ifr)/wcum(ifr)

        if(iblk.eq.1) then
          eerr=0
         else
          eerr=err(ecum(ifr),ecm2(ifr),ifr)
        endif

        if (ifr.eq.1) then                                 !JT
         eloc_bav = enow                                   !JT
         eloc_av = eave                                    !JT
         eloc_av_err = eerr                                !JT
         call object_modified_by_index (eloc_bav_index)    !JT
         call object_modified_by_index (eloc_av_index)     !JT
         call object_modified_by_index (eloc_av_err_index) !JT
        endif                                              !JT

        if(ndim.eq.2) then
          ieerr=nint(10000000*eerr)
         else
          ieerr=nint(100000*eerr)
        endif

        if(ifr.eq.1) then

          penow=pesum/wsum(ifr)
          peinow=peisum/wsum(ifr)
          tpbnow=tpbsum/wsum(ifr)
          tjfnow=tjfsum/wsum(ifr)
          r2now=r2sum/(wsum(ifr)*nelec)

          pecm2=pecm2+pesum*penow
          peicm2=peicm2+peisum*peinow
          tpbcm2=tpbcm2+tpbsum*tpbnow
          tjfcm2=tjfcm2+tjfsum*tjfnow
          r2cm2=r2cm2+r2sum*r2now/nelec

          pecum=pecum+pesum
          peicum=peicum+peisum
          tpbcum=tpbcum+tpbsum
          tjfcum=tjfcum+tjfsum
          r2cum=r2cum+r2sum/nelec
c         acccum=acccum+accsum
          if(index(mode,'mov1').eq.0) then
            acccum=acccum+accsum
           else
            acccum=acccum+accsum/nelec
          endif

!JT          if (.not. l_opt_lin .and. .not. l_opt_ptb) then  !JT
!JTc          call grad_hess_jas_cum(wsum(ifr),enow(ifr))
!JT           call grad_hess_jas_cum(wsum(ifr),enow)
!JT          endif                     !JT

          if(iblk.eq.1) then
            peerr=0
            peierr=0
            tpberr=0
            tjferr=0
           else
            peerr=err(pecum,pecm2,ifr)
            peierr=err(peicum,peicm2,ifr)
            tpberr=err(tpbcum,tpbcm2,ifr)
            tjferr=err(tjfcum,tjfcm2,ifr)
c           if(ibasis.eq.3) emerr=nelec*0.125*bext*bext*err(r2cum,r2cm2,ifr)
          endif

          peave=pecum/wcum(ifr)
          peiave=peicum/wcum(ifr)
          tpbave=tpbcum/wcum(ifr)
          tjfave=tjfcum/wcum(ifr)
          accave=acccum/wcum(ifr)

          if(ndim.eq.2) then
            ipeerr=nint(10000000*peerr)
            ipeierr=nint(10000000*peierr)
            itpber=nint(10000000*tpberr)
            itjfer=nint(10000000*tjferr)
           else
            ipeerr=nint(100000*peerr)
            ipeierr=nint(100000*peierr)
            itpber=nint(100000*tpberr)
            itjfer=nint(100000*tjferr)
          endif

          ierror_sigma=nint(100000*error_sigma)  !JT

c magnetic energy for quantum dots...
c right definition of the potential energy does not include magnetic energy.
          if(ndim.eq.2) then
c           emave=nelec*0.125*bext*bext*r2cum/wcum(ifr)+emaglz+emagsz+emagv
            temp=0.25d0*bext*bext/(we*we)
            emave=(peave-peiave-emag)*temp+emag
            emerr=(peerr+peierr)*temp
            iemerr=nint(10000000*emerr)
c           peave=peave-emave
            peave=peave*(1-temp)+(peiave+emag)*temp-emag
c           ipeerr=ipeerr+iemerr
            ipeerr=nint(10000000*(peerr*(1-temp)+temp*peierr))
          endif

          if(ndim.eq.2) then
            write(6,'(f12.7,5(f12.7,''('',i7,'')''),17x,f10.5,i10)')
     &      enow,eave,ieerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,emave,iemerr,
     &      accave,it
           elseif(nloc.eq.-3) then !MS Jellium
            write(6,'(f12.5,4(f12.5,''('',i5,'')''),17x,f10.5,i10,f10.5,''('',i5,'')'')')
     &      enow,eave,ieerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,
     &      accave,it,sigma,ierror_sigma
           else
            write(6,'(f10.5,4(f10.5,''('',i5,'')''),17x,f10.5,i10,f10.5,''('',i5,'')'')')
     &      enow,eave,ieerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,
     &      accave,it,sigma,ierror_sigma
          endif

c         if(it.ge.1000 .and. accave.lt.0.3d0)
          if(it.ge.1000 .and. accave.lt.0.1d0)
     &    write(6,'(''Warning: Low acceptance.  Are you sure you are running the mov1 version and deltas are OK?'')')
c         if(it.ge.1000 .and. accave.lt.0.1d0)
c    &    stop 'Low acceptance. Are you sure you are running the mov1 version and delta or deltar deltat are set correctly?'

         else

          fcum(ifr)=fcum(ifr)+wsum(1)*(enow-esum(1)/wsum(1))
          fcm2(ifr)=fcm2(ifr)+wsum(1)*(enow-esum(1)/wsum(1))**2
          fave=(ecum(1)/wcum(1)-ecum(ifr)/wcum(ifr))
          ferr=err(fcum(ifr),fcm2(ifr),1)
          if(deltot(ifr).ne.0.d0) then
            fave=fave/abs(deltot(ifr))
            ferr=ferr/abs(deltot(ifr))
          endif
          if(ndim.eq.2) then
            iferr=nint(10000000*ferr)
            write(6,'(f12.7,f12.7,''('',i7,'')'',51x,f12.7,''('',i7,'')'')')
     &      enow,eave,ieerr,fave,iferr
           else
            iferr=nint(100000*ferr)
            write(6,'(f10.5,f10.5,''('',i5,'')'',51x,f10.5,''('',i5,'')'')')
     &      enow,eave,ieerr,fave,iferr
          endif

        endif
   10 continue

c zero out xsum variables for metrop

      do 20 ifr=1,nforce
        esum(ifr)=0
  20    wsum(ifr)=0
      pesum=0
      peisum=0
      tpbsum=0
      tjfsum=0
      r2sum=0
      accsum=0

      call systemflush(6)

      return

      entry acues1
c statistical fluctuation (without blocking)
      ecum1=ecum1+esum1
      ecm21=ecm21+esum1**2
c     call grad_hess_jas_save
      esum1=0
      return

      entry acusig
c statistical fluctuation, sigma, (without blocking and without p,q)
      do 30 ifr=1,nforce
        ecum1s(ifr)=ecum1s(ifr)+esum1s(ifr)
        if(ifr.eq.1) then
          ecm21s(ifr)=ecm21s(ifr)+esum1s(ifr)**2
         else
          ecm21s(ifr)=ecm21s(ifr)+esum1s(ifr)**2/wsum1s(ifr)
        endif
   30   esum1s(ifr)=0
      return

      entry zerest

c entry point to zero out all averages etc.
c the initial values of energy psi etc. is also calculated here
c although that really only needs to be done before the equil. blocks.


!     allocations
      call alloc ('wsum1s', wsum1s, nforce)
      call alloc ('esum1s', esum1s, nforce)

      iblk=0

c set quadrature points

c     if(nloc.gt.0) call gesqua(nquad,xq,yq,zq,wq)
      if(nloc.gt.0) call rotqua

c zero out estimators

c     call wf_secondary

      pecum=0
      peicum=0
      tpbcum=0
      tjfcum=0
      r2cum=0
      acccum=0
      ecum1=0
c     ecum1s=0

      pecm2=0
      peicm2=0
      tpbcm2=0
      tjfcm2=0
      r2cm2=0
      ecm21=0
c     ecm21s=0

      pesum=0
      peisum=0
      tpbsum=0
      tjfsum=0
      r2sum=0
      accsum=0

      call grad_hess_jas_init

      call alloc ('fcum', fcum, nforce)
      call alloc ('fcm2', fcm2, nforce)
      call alloc ('wcum', wcum, nforce)
      call alloc ('wsum', wsum, nforce)
      call alloc ('ecm2', ecm2, nforce)
      call alloc ('ecum1s', ecum1s, nforce)
      call alloc ('ecm21s', ecm21s, nforce)
      call alloc ('ecum', ecum, nforce)
      call alloc ('esum', esum, nforce)
      do 65 ifr=1,nforce
        ecum1s(ifr)=0
        ecm21s(ifr)=0
        ecum(ifr)=0
        ecm2(ifr)=0
        wcum(ifr)=0
        fcum(ifr)=0
        fcm2(ifr)=0
        esum(ifr)=0
   65   wsum(ifr)=0

c Zero out estimators for acceptance, force-bias trun., kin. en. and density
      do 70 i=1,NRAD
        try(i)=0
        suc(i)=0
        trunfb(i)=0
        ekin(i)=0
        ekin2(i)=0
        rprobup(i)=0
        rprobdn(i)=0
   70   rprob(i)=0


c Zero out estimators for pair densities:
      if (ifixe.ne.0) then
      allocate (den2d_t(-NAX:NAX,-NAX:NAX))
      allocate (den2d_u(-NAX:NAX,-NAX:NAX))
      allocate (den2d_d(-NAX:NAX,-NAX:NAX))
      allocate (xx0probdt(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probdu(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probdd(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probut(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probuu(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probud(0:NAX,-NAX:NAX,-NAX:NAX))
      do 75 i2=-NAX,NAX
        do 75 i3=-NAX,NAX
          den2d_t(i2,i3)=0
          den2d_u(i2,i3)=0
          den2d_d(i2,i3)=0
          do 75 i1=0,NAX
            xx0probdt(i1,i2,i3)=0
            xx0probdu(i1,i2,i3)=0
            xx0probdd(i1,i2,i3)=0
            xx0probut(i1,i2,i3)=0
            xx0probuu(i1,i2,i3)=0
   75       xx0probud(i1,i2,i3)=0
      endif
      if (ifourier.ne.0) then
      allocate(fourierrk_t(0:NAX,0:NAK1))
      allocate(fourierrk_u(0:NAX,0:NAK1))
      allocate(fourierrk_d(0:NAX,0:NAK1))
      allocate(fourierkk_t(-NAK2:NAK2,-NAK2:NAK2))
      allocate(fourierkk_u(-NAK2:NAK2,-NAK2:NAK2))
      allocate(fourierkk_d(-NAK2:NAK2,-NAK2:NAK2))
      do 76 i1=0,NAX
        do 76 i2=0,NAK1
          fourierrk_t(i1,i2)=0
          fourierrk_u(i1,i2)=0
   76     fourierrk_d(i1,i2)=0
      do 77 i1=-NAK2,NAK2
        do 77 i2=-NAK2,NAK2
          fourierkk_t(i1,i2)=0
          fourierkk_u(i1,i2)=0
   77     fourierkk_d(i1,i2)=0
      endif


c get wavefunction etc. at initial point

c secondary configs
c set n- and e-coords and n-n potentials before getting wavefn. etc.
      do 80 ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido,psijo,vold,div_vo,d2,peo,peio,eold(ifr),denergy,ifr)
   80   psi2o(ifr)=2*(dlog(dabs(psido))+psijo)+dlog(ajacob)

c primary config
c set n- and e-coords and n-n potentials before getting wavefn. etc.
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido,psijo,vold,div_vo,d2,peo,peio,eold(1),denergy,1)
      psi2o(1)=2*(dlog(dabs(psido))+psijo)
      tjfo=d2
      tjfo=-tjfo*half*hb

!      if(nefp.gt.0) then
!        call startefp
!        call sample_efp(1,xold,eold(1),0.d0)
!        call efpsav
!      endif

      call grad_hess_jas_save

c get interparticle distances
c     call distances(xold,rvec_en,r_en,rvec_ee,r_ee,pe)
      call distances(xold,pe,pei)

c Find the minimum distance of each electron to any nucleus
c     do 86 i=1,nelec
c       rmino(i)=99.d9
c       do 85 j=1,ncent
c         dist=0
c         do 84 k=1,ndim
c  84       dist=dist+(xold(k,i)-cent(k,j))**2
c         if(dist.lt.rmino(i)) then
c           rmino(i)=dist
c           nearesto(i)=j
c         endif
c  85     continue
c       rmino(i)=dsqrt(rmino(i))
c       do 86  k=1,ndim
c  86     rvmino(k,i)=xold(k,i)-cent(k,nearesto(i))

      do 86 i=1,nelec
        rmino(i)=99.d9
        do 85 k=1,ncent
          if(r_en(i,k).lt.rmino(i)) then
            rmino(i)=r_en(i,k)
            nearesto(i)=k
          endif
   85     continue
        do 86  k=1,ndim
   86     rvmino(k,i)=rvec_en(k,i,nearesto(i))


      return
      end

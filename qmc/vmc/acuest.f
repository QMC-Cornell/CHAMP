      subroutine acuest
! Written by Cyrus Umrigar, modified by Claudia Filippi
! routine to accumulate estimators for energy etc.
      use all_tools_mod
      use constants_mod
      use control_mod
      use montecarlo_mod
      use optimization_mod
      use eloc_mod
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
      use zigzag_mod
      implicit real*8(a-h,o-z)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot

      dimension xstrech(3,nelec)
      dimension zznow(nzzvars)
! statement function for error calculation
!     err(x,x2)=dsqrt(dabs(x2/iblk-(x/iblk)**2)/iblk)
      err(x,x2,i)=dsqrt(abs(x2/wcum(i)-(x/wcum(i))**2)/iblk)

! xsum = sum of values of x from metrop
! xnow = average of values of x from metrop
! xcum = accumulated sums of xnow
! xcm2 = accumulated sums of xnow**2
! xave = current average value of x
! xerr = current error of x

      if(index(mode,'mpi').ne.0) call acuest_mpi

! Warning: At present we have nblk blocks of size nstep*nproc each.
! It would be a bit better to have nblk*nproc blocks of size nstep each,
! which just requires a bit of reorganizing.
      iblk=iblk+1
      it=iblk*nstep

! write out header first time
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

! write out current values of averages

      wsum(1)=dfloat(nstep)*dfloat(nproc)
!      if (l_reweight .and. l_opt .and. nforce .eq. 1) then
!         walker_weights_sum = walker_weights_sum + walker_weights_sum_block
!         l_reset_walker_weights_sum_block = .true.
!      else
!         walker_weights_sum_block = nstep_total
!         walker_weights_sum = dfloat(nstep_total)*dfloat(block_iterations_nb)
!         call object_modified_by_index (walker_weights_sum_block_index)
!      end if
!      call object_modified_by_index (walker_weights_sum_index)

      do 10 ifr=1,nforce
!       wnow=wsum(ifr)/nstep
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
!         eloc_bav = enow                                   !JT
!         eloc_av = eave                                    !JT
         eloc_av_err = eerr                                !JT
!         call object_modified_by_index (eloc_bav_index)    !JT
!         call object_modified_by_index (eloc_av_index)     !JT
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
          r1now=r1sum/(wsum(ifr)*nelec)
          r2now=r2sum/(wsum(ifr)*nelec)
          r3now=r3sum/(wsum(ifr)*nelec)
          r4now=r4sum/(wsum(ifr)*nelec)
          if(izigzag.gt.0) then
           zznow(:)=zzsum(:)/wsum(ifr)
          endif

          pecm2=pecm2+pesum*penow
          peicm2=peicm2+peisum*peinow
          tpbcm2=tpbcm2+tpbsum*tpbnow
          tjfcm2=tjfcm2+tjfsum*tjfnow
          r1cm2=r1cm2+r1sum*r1now/nelec
          r2cm2=r2cm2+r2sum*r2now/nelec
          r3cm2=r3cm2+r3sum*r3now/nelec
          r4cm2=r4cm2+r4sum*r4now/nelec
          if(izigzag.gt.0) then
            zzcm2(:)=zzcm2(:)+zzsum(:)*zznow(:)
          endif

          pecum=pecum+pesum
          peicum=peicum+peisum
          tpbcum=tpbcum+tpbsum
          tjfcum=tjfcum+tjfsum
          r1cum=r1cum+r1sum/nelec
          r2cum=r2cum+r2sum/nelec
          r3cum=r3cum+r3sum/nelec
          r4cum=r4cum+r4sum/nelec
          if(izigzag.gt.0) then
           zzcum(:)=zzcum(:)+zzsum(:)
          endif
!         acccum=acccum+accsum
          d_node_log_cum = d_node_log_cum + d_node_log_sum
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
!           if(ibasis.eq.3) emerr=nelec*0.125*bext*bext*err(r2cum,r2cm2,ifr)
          endif

          peave=pecum/wcum(ifr)
          peiave=peicum/wcum(ifr)
          tpbave=tpbcum/wcum(ifr)
          tjfave=tjfcum/wcum(ifr)
          accave=acccum/wcum(ifr)
          d_node_ave = exp(d_node_log_cum / wcum(ifr))

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

! magnetic energy for quantum dots...
! right definition of the potential energy does not include magnetic energy.
          if(ndim.eq.2) then
!           emave=nelec*0.125*bext*bext*r2cum/wcum(ifr)+emaglz+emagsz+emagv
            temp=0.25d0*bext*bext/(we*we)
            emave=(peave-peiave-emag)*temp+emag
            emerr=(peerr+peierr)*temp
            iemerr=nint(10000000*emerr)
!           peave=peave-emave
            peave=peave*(1-temp)+(peiave+emag)*temp-emag
!           ipeerr=ipeerr+iemerr
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

!         if(it.ge.1000 .and. accave.lt.0.3d0)
          if(it.ge.1000 .and. accave.lt.0.1d0)
     &    write(6,'(''Warning: Low acceptance.  Are you sure you are running the mov1 version and deltas are OK?'')')
!         if(it.ge.1000 .and. accave.lt.0.1d0)
!    &    stop 'Low acceptance. Are you sure you are running the mov1 version and delta or deltar deltat are set correctly?'

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

! zero out xsum variables for metrop

      do 20 ifr=1,nforce
        esum(ifr)=0
  20    wsum(ifr)=0
      pesum=0
      peisum=0
      tpbsum=0
      tjfsum=0
      r1sum=0
      r2sum=0
      r3sum=0
      r4sum=0
      accsum=0
      if(izigzag.gt.0) then
        zzsum(:)=0.d0
      endif
      d_node_log_sum = 0

      call systemflush(6)

      return

      entry acues1

!     during optimization calculate modified walker weights based on distance to node
      if (l_reweight .and. l_opt) then
         if (nforce .eq. 1) then
            call object_provide_by_index(vold_index)
            d_node = 1.d0 / sqrt(sum(vold**2))
            d_node_log = log(d_node)
            d_node_log_sum = d_node_log_sum + d_node_log
            current_walker_weight = d_node**reweight_power/(d_node**reweight_power+(d_node_ave/reweight_scale)**reweight_power)
!            if (l_reset_walker_weights_sum_block) then
!               walker_weights_sum_block = 0
!               l_reset_walker_weights_sum_block = .false.
!            end if
!            walker_weights_sum_block = walker_weights_sum_block + current_walker_weight
!            write(6,*) "!fp: ", esum1, d_node, current_walker_weight
         else
            current_walker_weight = 1.d0
!            walker_weights_sum_block = nstep_total
!            walker_weights_sum = dfloat(nstep_total)*dfloat(block_iterations_nb)
!            call object_modified_by_index (walker_weights_sum_index)
         end if
         call object_modified_by_index (current_walker_weight_index)
!         call object_modified_by_index (walker_weights_sum_block_index)
      end if

! statistical fluctuation (without blocking)
      ecum1=ecum1+esum1
      ecm21=ecm21+esum1**2
!     call grad_hess_jas_save
      esum1=0

      return

      entry acusig
! statistical fluctuation, sigma, (without blocking and without p,q)
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

! entry point to zero out all averages etc.
! the initial values of energy psi etc. is also calculated here
! although that really only needs to be done before the equil. blocks.


!     allocations
      call alloc ('wsum1s', wsum1s, nforce)
      call alloc ('esum1s', esum1s, nforce)

      iblk=0

! set quadrature points

!     if(nloc.gt.0) call gesqua(nquad,xq,yq,zq,wq)
      if(nloc.gt.0) call rotqua

! zero out estimators

!     call wf_secondary

      pecum=0
      peicum=0
      tpbcum=0
      tjfcum=0
      r1cum=0
      r2cum=0
      r3cum=0
      r4cum=0
      if(izigzag.gt.0) then
        zzcum(:)=0.d0
      endif
      acccum=0
      ecum1=0
!     ecum1s=0
      d_node_log_cum = 0
!      walker_weights_sum = 0
!      call object_modified_by_index (walker_weights_sum_index)

      pecm2=0
      peicm2=0
      tpbcm2=0
      tjfcm2=0
      r1cm2=0
      r2cm2=0
      r3cm2=0
      r4cm2=0
      if(izigzag.gt.0) then
        zzcm2(:)=0.d0
      endif
      ecm21=0
!     ecm21s=0

      pesum=0
      peisum=0
      tpbsum=0
      tjfsum=0
      r1sum=0
      r2sum=0
      r3sum=0
      r4sum=0
      if(izigzag.gt.0) then
        zzsum(:)=0.d0
      endif
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

! Zero out estimators for acceptance, force-bias trun., kin. en. and density
      do 70 i=1,NRAD
        try(i)=0
        suc(i)=0
        trunfb(i)=0
        ekin(i)=0
        ekin2(i)=0
        rprobup(i)=0
        rprobdn(i)=0
   70   rprob(i)=0


! Zero out estimators for pair densities:
      if (ifixe.ne.0) then
      call alloc_range ('den2d_t', den2d_t, -NAX, NAX, -NAX, NAX)
      call alloc_range ('den2d_u', den2d_u, -NAX, NAX, -NAX, NAX)
      call alloc_range ('den2d_d', den2d_d, -NAX, NAX, -NAX, NAX)
      call alloc_range ('pot_ee2d_t', pot_ee2d_t, -NAX, NAX, -NAX, NAX)
      call alloc_range ('pot_ee2d_u', pot_ee2d_u, -NAX, NAX, -NAX, NAX)
      call alloc_range ('pot_ee2d_d', pot_ee2d_d, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probdt', xx0probdt, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probdu', xx0probdu, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probdd', xx0probdd, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probut', xx0probut, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probuu', xx0probuu, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probud', xx0probud, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('dos', dos, 0, NAX)
      dos(:) = 0
      do 75 i2=-NAX,NAX
        do 75 i3=-NAX,NAX
          do 75 i1=0,NAX
            den2d_t(i2,i3)=0
            den2d_u(i2,i3)=0
            den2d_d(i2,i3)=0
            pot_ee2d_t(i2,i3)=0
            pot_ee2d_u(i2,i3)=0
            pot_ee2d_d(i2,i3)=0
            xx0probdt(i1,i2,i3)=0
            xx0probdu(i1,i2,i3)=0
            xx0probdd(i1,i2,i3)=0
            xx0probut(i1,i2,i3)=0
            xx0probuu(i1,i2,i3)=0
   75       xx0probud(i1,i2,i3)=0
      endif
      if (ifourier.ne.0) then
      call alloc_range ('fourierrk_t', fourierrk_t, -NAX, NAX, 0, NAK1)
      call alloc_range ('fourierrk_u', fourierrk_u, -NAX, NAX, 0, NAK1)
      call alloc_range ('fourierrk_d', fourierrk_d, -NAX, NAX, 0, NAK1)
      call alloc_range ('fourierkk_t', fourierkk_t, -NAK2, NAK2, -NAK2, NAK2)
      call alloc_range ('fourierkk_u', fourierkk_u, -NAK2, NAK2, -NAK2, NAK2)
      call alloc_range ('fourierkk_d', fourierkk_d, -NAK2, NAK2, -NAK2, NAK2)
      do 76 i1=-NAX,NAX
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
      if (izigzag.ne.0) then
        call alloc_range ('zzpairden_t', zzpairden_t, -NAX, NAX, -NAX, NAX)
        call alloc_range ('zzpairdenij_t', zzpairdenij_t, -NAX, NAX, 0, (nelec-1))
        call alloc_range ('zzcorr', zzcorr, 0, NAX)
        call alloc_range ('zzcorr1', zzcorr1, 0, NAX)
        call alloc_range ('zzcorr2', zzcorr2, 0, NAX)
        call alloc_range ('zzcorrij', zzcorrij, 0, (nelec-1))
        call alloc_range ('zzcorrij1', zzcorrij1, 0, (nelec-1))
        call alloc_range ('zzcorrij2', zzcorrij2, 0, (nelec-1))
        call alloc_range ('yycorr', yycorr, 0, NAX)
        call alloc_range ('yycorr1', yycorr1, 0, NAX)
        call alloc_range ('yycorr2', yycorr2, 0, NAX)
        call alloc_range ('yycorrij', yycorrij, 0, (nelec-1))
        call alloc_range ('yycorrij1', yycorrij1, 0, (nelec-1))
        call alloc_range ('yycorrij2', yycorrij2, 0, (nelec-1))
        call alloc_range ('znncorr', znncorr, 0, NAX)
        call alloc_range ('zn2ncorr', zn2ncorr, 0, NAX)
        zzpairden_t(:,:) = 0
        zzpairdenij_t(:,:) = 0
        zzcorr(:) = 0
        zzcorr1(:) = 0
        zzcorr2(:) = 0
        zzcorrij(:) = 0
        zzcorrij1(:) = 0
        zzcorrij2(:) = 0
        yycorr(:) = 0
        yycorr1(:) = 0
        yycorr2(:) = 0
        yycorrij(:) = 0
        yycorrij1(:) = 0
        yycorrij2(:) = 0
        znncorr(:) = 0
        zn2ncorr(:) = 0
      endif

! get wavefunction etc. at initial point

! secondary configs
! set n- and e-coords and n-n potentials before getting wavefn. etc.
      do 80 ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido,psijo,vold,div_vo,d2,peo,peio,eold(ifr),denergy,ifr)
   80   psi2o(ifr)=2*(dlog(dabs(psido))+psijo)+dlog(ajacob)

! primary config
! set n- and e-coords and n-n potentials before getting wavefn. etc.
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

! get interparticle distances
!     call distances(xold,rvec_en,r_en,rvec_ee,r_ee,pe)
      call distances(xold,pe,pei)

! Find the minimum distance of each electron to any nucleus
!     do 86 i=1,nelec
!       rmino(i)=99.d9
!       do 85 j=1,ncent
!         dist=0
!         do 84 k=1,ndim
!  84       dist=dist+(xold(k,i)-cent(k,j))**2
!         if(dist.lt.rmino(i)) then
!           rmino(i)=dist
!           nearesto(i)=j
!         endif
!  85     continue
!       rmino(i)=dsqrt(rmino(i))
!       do 86  k=1,ndim
!  86     rvmino(k,i)=xold(k,i)-cent(k,nearesto(i))

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

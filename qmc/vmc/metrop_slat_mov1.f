      subroutine metrop_slat_mov1(ipass)
c     subroutine metrop6(ipass)
c Written by Cyrus Umrigar
c Uses the accelerated Metropolis method described in:
c 1) Accelerated Metropolis Method, C.J. Umrigar, PRL 71, 408 (1993).
c 2) Variational Monte Carlo Basics and Applications to Atoms and Molecules,
c    C.J. Umrigar, in {\it Quantum Monte Carlo Methods in Physics and Chemistry},
c    edited by M.~P. Nightingale and C.~J. Umrigar. NATO ASI Series, Series C,
c    Mathematical and Physical Sciences, Vol. C-525,
c    (Kluwer Academic Publishers, Boston, 1999)

      use all_tools_mod
      use montecarlo_mod
      use opt_lin_mod
      use opt_ptb_mod
      use deriv_exp_mod
      use atom_mod
      use config_mod
      use dets_mod
      use const_mod
      use const2_mod
      use dim_mod
      use forcepar_mod
      use doefp_mod
      use pseudo_mod
      use delocc_mod
      implicit real*8(a-h,o-z)


c     character*16 mode

!JT      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0,half=0.5d0)
      parameter (d3b2=1.5d0,d5b2=2.5d0,d2b3=.666666666666667d0)
      parameter (eps=1.d-10)
c     parameter (g3b2=.886226925452758d0)
      parameter (g5b2=1.329340388179137d0)
c g3b2, g5b2 are gamma3/2), gamma(5/2)


c The moves are now being made in local r,theta phi coordinates.

c The foll. additions have been made:
c 1) Slater form of Tij.
c 2) Make theta_max a function of r
c 3) Generalize to molecules. This requires geometric rejections.

c The foll. still need to be tried:
c 1) Quadratic, gaussian, Morse and Exp(-zeta*r)+co*Exp(-r) forms of Tij
c    Last 2 are prob. best

!JT      common /dim/ ndim
c     common /contr3/ mode
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /const2/ deltar,deltat
!JT      common /config/ xold(3,MELEC),xnew(3,MELEC),vold(3,MELEC)
!JT     &,vnew(3,MELEC),psi2o(MFORCE),psi2n(MFORCE),eold(MFORCE),enew(MFORCE)
!JT     &,peo,pen,peio,pein,tjfn,tjfo,psido,psijo
!JT     &,rmino(MELEC),rminn(MELEC),rvmino(3,MELEC),rvminn(3,MELEC)
!JT     &,rminon(MELEC),rminno(MELEC),rvminon(3,MELEC),rvminno(3,MELEC)
!JT     &,nearesto(MELEC),nearestn(MELEC),delttn(MELEC)
!JT      common /delocc/ denergy(MPARM)
      common /estsum/ esum1,esum(MFORCE),pesum,peisum,tpbsum,tjfsum,r2sum,accsum
      common /estsig/ wsum1s(MFORCE),esum1s(MFORCE),ecum1s(MFORCE),ecm21s(MFORCE)
      common /stats/ rejmax
      common /stepv/try(NRAD),suc(NRAD),trunfb(NRAD),rprob(NRAD),
     &ekin(NRAD),ekin2(NRAD)
      common /denupdn/ rprobup(NRAD),rprobdn(NRAD)
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /kinet/ ekineo(MELEC),ekinen(MELEC)
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
!JT     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
!JT      common /doefp/ nefp
      common /div_v/ div_vo(MELEC)
!JT      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcewt/ wsum(MFORCE),wcum(MFORCE)
      common /jel_sph2/ zconst ! RM

      dimension xaxis(3),yaxis(3),zaxis(3),idist(MELEC)

      dimension xstrech(3,MELEC)

c     area(ri,r1,r2,v)=dabs((one/sqrt(ri))*
c    &(r2**d3b2*(two*(one-v*ri)/3+.4d0*v*r2)
c    &-r1**d3b2*(two*(one-v*ri)/3+.4d0*v*r1)))

      thetamx(r,z)=deltat+(two-deltat)/(one+(z*r)**2)

c     mode='vmc_mov1    '

      deltri=one/deltar

!JT: check local energy
!      call hpsi(xold,psido,psijo,vold,div_vo,d2o,peo,peio,eold(1),denergy,1)
!      psi2o(1)=2*(dlog(dabs(psido))+psijo)
!      write(6,*) 'xold=',xold
!      write(6,*) 'eold=',eold(1)
!JT

      do 300 i=1,nelec
        fxop=one
        nearo=nearesto(i)
        if(nloc.eq.0) then
          zcusp=znuc(iwctype(nearo)) ! RM (watch out here for future extension)
!MS Jellium sphere
         else if(nloc.eq.-3) then
          zcusp=znuc(iwctype(nearo))
         else
          zcusp=zero
        endif
c Although rmino is saved, recalculate it, otherwise there may be a
c build up of errors via zaxis. I don't think so.
        rmino(i)=dsqrt(rvmino(1,i)**2+rvmino(2,i)**2+rvmino(3,i)**2)
c Choose lower and upper values of r sampling
        rbot=rmino(i)*deltri
        rtop=rmino(i)*deltar
c Calculate magnitude of the velocity in the radial direction
        voldr=zero
        do 10 k=1,ndim
   10     voldr=voldr+vold(k,i)*rvmino(k,i)
        voldr=voldr/rmino(i)

c Place x-axis along direction of angular change and
c Calculate the velocity in the phi direction
        voldp=zero
        do 20 k=1,ndim
          zaxis(k)=rvmino(k,i)/rmino(i)
          xaxis(k)=vold(k,i)-voldr*zaxis(k)
   20     voldp=voldp+xaxis(k)**2
        voldp=dsqrt(voldp)
        if(voldp.lt.eps) then
          xaxis(1)=eps*(one-zaxis(1)**2)
          xaxis(2)=eps*(-zaxis(1)*zaxis(2))
          xaxis(3)=eps*(-zaxis(1)*zaxis(3))
          voldp=eps*sqrt(one-zaxis(1)**2)
        endif
        do 30 k=1,ndim
   30     xaxis(k)=xaxis(k)/voldp

c Limit radial component of velocity.
c It may be a good idea to limit it if it is positive too.
        if(zconst.ne.0)then ! Jellium RM
          voldr=max(voldr,-2*zconst)
         else
          voldr=max(voldr,-2*znuc(iwctype(nearo)))
c         voldr=min(voldr,2*znuc(iwctype(nearo)))
        endif

c y-axis is cross-product of z and x axes
        yaxis(1)=zaxis(2)*xaxis(3)-zaxis(3)*xaxis(2)
        yaxis(2)=zaxis(3)*xaxis(1)-zaxis(1)*xaxis(3)
        yaxis(3)=zaxis(1)*xaxis(2)-zaxis(2)*xaxis(1)

c Temporary test of fbias
c       voldr=voldr*fbias
        voldp=voldp*fbias

        root=(zcusp+voldr)*(zcusp+voldr-four/rmino(i))
        if(root.ge.zero) then
          root=sqrt(root)
          zeta=half*(zcusp-voldr-root)
          if(zeta.le.zero) zeta=half*(zcusp-voldr+root)
          if(zeta.le.zero) zeta=eps
         else
          if(voldr.lt.zero) then
            zeta=-voldr-eps
           else
            zeta=one
          endif
        endif
        co=(zeta+voldr)/(one-(zeta+voldr)*rmino(i))

c       write(6,'(''rmino(i),voldr,zeta,co='',9f10.5)')
c    &  rmino(i),voldr,zeta,co,(co-zeta-co*zeta*rmino(i))/
c    &  (one+co*rmino(i))

c Use Slater approx for radial fn
c Determine the maximum value of radial function for rejection sampling
        root=dsqrt((d3b2*co-zeta)**2+two*zeta*co)
! JT: add test on co for special case such as H atom for which co=0 and rmax1=infinity
! rmax1 was being limited to the interval [rbot,rtop] in the next lines but may be the compiler does not like comparing +-infinity.
        if(co.ne.0d0) then ! JT
          rmax1=(d3b2*co-zeta+root)/(two*zeta*co)
         else              ! JT
          rmax1=rtop       ! JT
        endif              ! JT

        if(rmax1.lt.rbot) rmax1=rbot
        if(rmax1.gt.rtop) rmax1=rtop
        fmax=sqrt(rmax1)*abs(one+co*rmax1)*dexp(-zeta*rmax1)
        if(co.lt.zero) then
          rmax2=(d3b2*co-zeta-root)/(two*zeta*co)
          if(rmax2.lt.rbot) rmax2=rbot
          if(rmax2.gt.rtop) rmax2=rtop
          fmax2=sqrt(rmax2)*abs(one+co*rmax2)*dexp(-zeta*rmax2)
          fmax=max(fmax,fmax2)
        endif

c   Sample sqrt(r_f)*abs(1+co*r_f)*exp(-zeta*r_f) by rejection
        bot=sqrt(rmino(i))*abs(one+co*rmino(i))*dexp(-zeta*rmino(i))
   40   rtry=((deltar-deltri)*rannyu(0)+deltri)*rmino(i)
          top=sqrt(rtry)*abs(one+co*rtry)*dexp(-zeta*rtry)
          ratio=top/fmax
          rejmax=max(rejmax,ratio)
          if(ratio.gt.rannyu(0)) goto 50
        goto 40
   50   fxop=fxop*top/bot

c   Calculate the integral of T
        rzero=-one/co
        zrbot=zeta*rbot
        zrtop=zeta*rtop
        zebot=dsqrt(zrbot**3)*dexp(-zrbot)
        zetop=dsqrt(zrtop**3)*dexp(-zrtop)
        g52bot=gammai(d5b2,zrbot,zrbot*zebot,iflagb)
        g52top=gammai(d5b2,zrtop,zrtop*zetop,iflagt)
        if(rzero.lt.rbot .or. rzero.gt.rtop) then
          if(iflagb*iflagt.eq.1) then
            g52dif=g52top-g52bot
           else
            g52dif=g52top-g52bot+g5b2
          endif
          g32dif=d2b3*(g52dif+zetop-zebot)
          areao=dabs(g32dif+co*g52dif/zeta)/(bot*dsqrt(zeta**3))
         else
          zrzer=zeta*rzero
          zezer=dsqrt(zrzer**3)*dexp(-zrzer)
          g52zer=gammai(d5b2,zrzer,zrzer*zezer,iflagz)
          if(iflagb*iflagz.eq.1) then
            g52dif1=g52zer-g52bot
           else
            g52dif1=g52zer-g52bot+g5b2
          endif
          g32dif1=d2b3*(g52dif1+zezer-zebot)
          if(iflagt*iflagz.eq.1) then
            g52dif2=g52top-g52zer
           else
            g52dif2=g52top-g52zer+g5b2
          endif
          g32dif2=d2b3*(g52dif2+zetop-zezer)
          areao=(dabs(g32dif1+co*g52dif1/zeta)
     &          +dabs(g32dif2+co*g52dif2/zeta))/
     &          (bot*dsqrt(zeta**3))
        endif

c Sample cos(theta)
        raver=half*(rmino(i)+rtry)
        deltt=thetamx(raver,znuc(iwctype(nearo)))
        if(zconst.ne.0) deltt=thetamx(raver,zconst) ! Jellium RM
        costht=one-deltt*rannyu(0)
        zprime=rtry*costht
        sintht=dsqrt(one-costht*costht)

c For molecules deltt may not be the same for forward and reverse
c moves, so it is necessary to include this in areao and arean
        areao=areao*deltt

c Truncate phi variation if it goes through zero
c Sample phi by rejection. Note it is OK to have a) theta or sin(theta)
c and b) rtry/rold(i) or raver, as long as the forward and reverse probs.
c are consistent.
c If we do not limit term to be <=1 then use commented out lines.
        term=dmin1(voldp*raver*sintht,one)
clim    term=voldp*raver*sintht
        fmax=one+term
   60   phitry=pi*rannyu(0)
          cosphi=dcos(phitry)
          top=one+term*cosphi
clim      top=abs(one+term*cosphi)
          if(top.gt.rannyu(0)*fmax) goto 70
        goto 60
   70   fxop=fxop*top

clim    if(term.gt.one) then
clim      phizer=dacos(-one/term)
clim      areao=areao*((two/pi)*(phizer+term*dsin(phizer))-one)
clim    endif

c Calculate x and y coordinates in local coordinate system
        xprime=rtry*sintht*dcos(phitry)
        yprime=dsqrt(max(zero,rtry*rtry-xprime**2-zprime**2))
        if(rannyu(0).lt.half) yprime=-yprime

c Convert back to original coordinate system
        do 80 k=1,ndim
          rvminno(k,i)=xaxis(k)*xprime+yaxis(k)*yprime+zaxis(k)
     &    *zprime
   80     xnew(k,i)=rvminno(k,i)+cent(k,nearo)
        rminno(i)=rtry

c Do geometrical rejections for molecules
        rminn(i)=99.d9
        do 85 j=1,ncent
          dist=zero
          do 84 k=1,ndim
   84       dist=dist+(xnew(k,i)-cent(k,j))**2
          if(dist.lt.rminn(i)) then
            rminn(i)=dist
            nearestn(i)=j
          endif
   85     continue
        nearn=nearestn(i)
        rminn(i)=dsqrt(rminn(i))
        rminon(i)=zero
        dot=zero
        do 86  k=1,ndim
          rvminn(k,i)=xnew(k,i)-cent(k,nearestn(i))
          rvminon(k,i)=xold(k,i)-cent(k,nearestn(i))
          rminon(i)=rminon(i)+rvminon(k,i)**2
   86     dot=dot+rvminn(k,i)*rvminon(k,i)
        rminon(i)=dsqrt(rminon(i))
        dot=dot/(rminn(i)*rminon(i))
        costht=dot
        sintht=dsqrt(one-costht*costht)
        ravern=half*(rminn(i)+rminon(i))
        delttn(i)=thetamx(ravern,znuc(iwctype(nearestn(i))))
        if(zconst.ne.0) delttn(i)=thetamx(ravern,zconst) ! Jellium RM
        igeometrical=0
        if(rminon(i).gt.rminn(i)*deltar .or. dot.lt.one-delttn(i))then
          igeometrical=1
          p=zero
          q=one
          goto 208
        endif

c rratio^2 is needed for the density of the angular moves
        rratio=rminno(i)/rminon(i)

        if(ipr.ge.1) then
          rtest=dsqrt(rvminn(1,i)**2+rvminn(2,i)**2+rvminn(3,i)**2)
          rtest2=dsqrt(xprime**2+yprime**2+zprime**2)
          write(6,'(''rtest,rtest2,rtry'',9d14.6)')rtest,rtest2,rtry,
     &    rtest-rtry,rtest2-rtry
          write(6,'(''vold='',9d12.4)') (vold(k,i),k=1,ndim)
          write(6,'(''voldr,voldp='',9d12.4)') voldr,voldp
          write(6,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
          if(co.lt.zero) write(6,'(''rmino(i),rmax1,rmax2,rzero'',9f9.4)')
     &    rmino(i),rmax1,rmax2,rzero
          write(6,'(''rtry,costht,sintht,phitry'',9f9.4)') rtry,costht,
     &    sintht,phitry
          write(6,'(''fxop'',9f12.4)') fxop
        endif

c calculate psi at new configuration
      iel=i
      call hpsie(iel,xnew,psidn,psijn,vnew)
      psi2n(1)=2*(dlog(dabs(psidn))+psijn)

c calculate probability for reverse transition
      fxnp=one
c Choose lower and upper values of r sampling
        rbot=rminn(i)*deltri
        rtop=rminn(i)*deltar
c Calculate magnitude of the velocity in the radial direction
        vnewr=zero
        do 110 k=1,ndim
  110     vnewr=vnewr+vnew(k,i)*rvminn(k,i)
        vnewr=vnewr/rminn(i)

c Place x-axis along direction of angular change and
c Calculate the velocity in the phi direction
        vnewp=zero
        do 120 k=1,ndim
          xaxis(k)=vnew(k,i)-vnewr*rvminn(k,i)/rminn(i)
  120     vnewp=vnewp+xaxis(k)**2
        vnewp=dsqrt(vnewp)
        if(vnewp.lt.eps) then
          xaxis(1)=eps*(one-rvminn(1,i)**2)
          xaxis(2)=eps*(-rvminn(1,i)*rvminn(2,i))
          xaxis(3)=eps*(-rvminn(1,i)*rvminn(3,i))
          vnewp=eps*sqrt(one-rvminn(1,i)**2*(two-rminn(i)**2))
        endif
        do 130 k=1,ndim
  130     xaxis(k)=xaxis(k)/vnewp

c Limit radial component of velocity.
c It may be a good idea to limit it if it is positive too.
        nearn=nearestn(i)
        if(nloc.eq.0) then
          zcusp=znuc(iwctype(nearn))
         else
          zcusp=zero
        endif
        if(zconst.ne.0) then ! Jellium RM
          vnewr=max(vnewr,-2*zconst)
         else
          vnewr=max(vnewr,-2*znuc(iwctype(nearn)))
c         vnewr=min(vnewr,2*znuc(iwctype(nearn)))
        endif

c Temporary test of fbias
c       vnewr=vnewr*fbias
        vnewp=vnewp*fbias

        root=(zcusp+vnewr)*(zcusp+vnewr-four/rminn(i))
        if(root.ge.zero) then
          root=sqrt(root)
          zeta=half*(zcusp-vnewr-root)
          if(zeta.le.zero) zeta=half*(zcusp-vnewr+root)
          if(zeta.le.zero) zeta=eps
         else
          if(vnewr.lt.zero) then
            zeta=-vnewr-eps
           else
            zeta=one
          endif
        endif
        co=(zeta+vnewr)/(one-(zeta+vnewr)*rminn(i))

        bot=sqrt(rminn(i))*abs(one+co*rminn(i))*dexp(-zeta*rminn(i))
        top=sqrt(rminon(i))*abs(one+co*rminon(i))*dexp(-zeta*rminon(i))
        fxnp=fxnp*top/bot

c Calculate the integral of T
        rzero=-one/co
        zrbot=zeta*rbot
        zrtop=zeta*rtop
        zebot=dsqrt(zrbot**3)*dexp(-zrbot)
        zetop=dsqrt(zrtop**3)*dexp(-zrtop)
        g52bot=gammai(d5b2,zrbot,zrbot*zebot,iflagb)
        g52top=gammai(d5b2,zrtop,zrtop*zetop,iflagt)
        if(rzero.lt.rbot .or. rzero.gt.rtop) then
          if(iflagb*iflagt.eq.1) then
            g52dif=g52top-g52bot
           else
            g52dif=g52top-g52bot+g5b2
          endif
          g32dif=d2b3*(g52dif+zetop-zebot)
          arean=dabs(g32dif+co*g52dif/zeta)/(bot*dsqrt(zeta**3))
         else
          zrzer=zeta*rzero
          zezer=dsqrt(zrzer**3)*dexp(-zrzer)
          g52zer=gammai(d5b2,zrzer,zrzer*zezer,iflagz)
          if(iflagb*iflagz.eq.1) then
            g52dif1=g52zer-g52bot
           else
            g52dif1=g52zer-g52bot+g5b2
          endif
          g32dif1=d2b3*(g52dif1+zezer-zebot)
          if(iflagt*iflagz.eq.1) then
            g52dif2=g52top-g52zer
           else
            g52dif2=g52top-g52zer+g5b2
          endif
          g32dif2=d2b3*(g52dif2+zetop-zezer)
          arean=(dabs(g32dif1+co*g52dif1/zeta)
     &          +dabs(g32dif2+co*g52dif2/zeta))/
     &          (bot*dsqrt(zeta**3))
        endif

c For molecules deltt may not be the same for forward and reverse
c moves, so it is necessary to include this in areao and arean
        arean=arean*delttn(i)

c Truncate phi variation if it goes through zero
c Note it is OK to have a) theta or sin(theta)
c and b) rtry/rold(i) or raver, as long as the forward and reverse probs.
c are consistent.
c If we do not limit term to be <=1 then use commented out lines.
        term=dmin1(vnewp*ravern*sintht,one)
clim    term=vnewp*raver*sintht
c Determine cos(phi)
        cosphi=zero
        rnorm=zero
        do 160 k=1,ndim
          term2=rvminon(k,i)/rminon(i)-costht*rvminn(k,i)/rminn(i)
          rnorm=rnorm+term2*term2
  160     cosphi=cosphi+term2*xaxis(k)
        cosphi=cosphi/dsqrt(rnorm)
        fxnp=fxnp*(one+term*cosphi)
clim    fxnp=fxnp*abs(one+term*cosphi)

clim    if(term.gt.one) then
clim      phizer=dacos(-one/term)
clim      arean=arean*((two/pi)*(phizer+term*dsin(phizer))-one)
clim    endif

c p is the probability of accepting new move
      p=rratio**2*exp(psi2n(1)-psi2o(1))*dabs((fxnp*areao)/(fxop*arean))

        if(ipr.ge.1) then
          write(6,'(''rminn,rvminn,vnew,vnewr'',9f10.4)')
     &    rminn(i),(rvminn(k,i),k=1,ndim),(vnew(k,i),k=1,ndim),vnewr
          write(6,'(''vnew='',9d12.4)') (vnew(k,i),k=1,ndim)
          write(6,'(''vnewr,vnewp='',9d12.4)') vnewr,vnewp
          write(6,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
          if(co.lt.zero) write(6,'(''rminn(i),rmax1,rmax2,rzero'',9f9.4)')
     &    rminn(i),rmax1,rmax2,rzero
          write(6,'(''rtry,costht,sintht,phitry,cosphi'',9f9.4)') rtry,
     &    costht,sintht,phitry,cosphi
          write(6,'(''1:fxop,fxnp,areao,arean,psi2n,psi2o,p'',9f9.4)')
     &                  fxop,fxnp,areao,arean,psi2n(1),psi2o(1),p
        endif

      p=dmin1(one,p)
      q=one-p

  208 continue
c Calculate as a function of the distance to the nucleus
c 1) acceptance,  2) force-bias truncation probability,
c 3) kinetic energy and it's fluctuation
c The K.E. is not quite correct, since we should use p times new
c and q times old, and keep track of which bin the old was in
      rold=dsqrt(xold(1,i)**2+xold(2,i)**2+xold(3,i)**2)
      rnew=dsqrt(xnew(1,i)**2+xnew(2,i)**2+xnew(3,i)**2)
      itryo=min(int(delri*rold)+1,NRAD)
      itryn=min(int(delri*rnew)+1,NRAD)
      try(itryo)=try(itryo)+1
      suc(itryo)=suc(itryo)+p
      if(try(itryo).lt.0.d0) write(6,'(''itryo,try'',i5,d13.5)')itryo,
     &try(itryo)
c     if(suc(itryo).lt.0.d0) write(6,'(''itryo,suc'',i5,9d13.5)')itryo,
c    &suc(itryo),rratio,psi2n(1),psi2o(1),fxnp,areao,fxop,arean
      if(voldp*raver*sintht.gt.one) trunfb(itryo)=trunfb(itryo)+1
      if(i.le.nup) then
        rprobup(itryo)=rprobup(itryo)+q
        rprobup(itryn)=rprobup(itryn)+p
       else
        rprobdn(itryo)=rprobdn(itryo)+q
        rprobdn(itryn)=rprobdn(itryn)+p
      endif
      rprob(itryo)=rprob(itryo)+q
      rprob(itryn)=rprob(itryn)+p
      do 210 k=1,ndim
  210   r2sum=r2sum+p*xnew(k,i)**2+q*xold(k,i)**2

c accept new move with probability p
c Note when one electron moves the velocity on all electrons change.
      if(rannyu(0).lt.p) then
        idist(i)=itryn
        rmino(i)=rminn(i)
        nearesto(i)=nearestn(i)
        psi2o(1)=psi2n(1)
        do 240 k=1,ndim
          xold(k,i)=xnew(k,i)
          rvmino(k,i)=rvminn(k,i)
          do 240 ii=1,nelec
  240       vold(k,ii)=vnew(k,ii)
        psido=psidn
        psijo=psijn
        accsum=accsum+one
        call jassav(i)
        call detsav(i)
       else
        idist(i)=itryo
        do 250 k=1,ndim
  250     xnew(k,i)=xold(k,i)
c       call distancese_restore(i,rvec_en,r_en,rvec_ee,r_ee)
        if(igeometrical.eq.0) call distancese_restore(i)
      endif

  300 continue

      call object_modified_by_index (xold_index)  !JT

c loop over secondary configurations
      do 350 ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido,psijo,vold,div_vo,d2,peo,peio,eold(ifr),denergy,ifr)
  350   psi2o(ifr)=2*(dlog(dabs(psido))+psijo)+dlog(ajacob)

c primary configuration
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido,psijo,vold,div_vo,d2o,peo,peio,eold(1),denergy,1)
      psi2o(1)=2*(dlog(dabs(psido))+psijo)

      eloc = eold(1)                               !JT
      call object_modified_by_index (eold_index)   !JT
      call object_modified_by_index (eloc_index)   !JT
      call object_modified_by_index (psido_index)  !JT
      call object_modified_by_index (psijo_index)  !JT
      call object_modified_by_index (denergy_index)!JT
      call object_modified_by_index (vold_index)   !JT
      call object_modified_by_index (div_vo_index) !JT

! check deloc_exp !!!!!!!!!!!!!!!!!!!!!!!!!
!      if (l_deloc_exp_num) then
!      iexp=1
!      dzex = 0.001
!      eloc1 = eold(1)
!      call move_zex (iexp, dzex)
!      call hpsi(xold,psido,psijo,vold,div_vo,d2o,peo,eold(1),denergy,1)
!      call object_modified_by_index (eold_index)  !JT
!      call object_modified_by_index (eloc_index)  !JT
!      call object_modified_by_index (psido_index) !JT
!      call object_modified_by_index (denergy_index) !JT
!      call object_modified_by_index (vold_index) !JT
!      call object_modified_by_index (div_vo_index) !JT
!      eloc2 = eold(1)
!      call move_zex (iexp, -dzex)
!      call hpsi(xold,psido,psijo,vold,div_vo,d2o,peo,eold(1),denergy,1)
!      call object_modified_by_index (eold_index)  !JT
!      call object_modified_by_index (eloc_index)  !JT
!      call object_modified_by_index (psido_index) !JT
!      call object_modified_by_index (denergy_index) !JT
!      call object_modified_by_index (vold_index) !JT
!      call object_modified_by_index (div_vo_index) !JT
!
!      write (6,'(a,f)') 'eloc1=',eloc1
!      write (6,'(a,f)') 'eloc2=',eloc2
!      write (6,'(a,f)') 'deloc_exp_num =',(eloc2-eloc1)/dzex
!      call object_provide ('deloc_exp')
!      write (6,'(a,i3,a,f)') 'deloc_exp(',iexp,')=', deloc_exp(iexp)
!      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      tjfo=-d2o*half*hb
c form expected values of e, pe, etc.
      esum1=eold(1)
      esum(1)=esum(1)+eold(1)
      pesum=pesum+peo
      peisum=peisum+peio
      tpbsum=tpbsum+(eold(1)-peo)
      tjfsum=tjfsum+tjfo
      if(nefp.gt.0) then
        call sample_efp(1,xold,eold(1),1.d0)
        call efpsav
      endif

      call grad_hess_jas_sum(1.d0,0.d0,eold(1),eold(1),1.d0,0.d0)

      call acues1

      do 380 ifr=1,nforce
        if(ifr.eq.1) then
          esum1s(ifr)=eold(ifr)
         else
          wstro=exp(psi2o(ifr)-psi2o(1))
          wsum1s(ifr)=wstro
          esum(ifr)=esum(ifr)+eold(ifr)*wstro
          wsum(ifr)=wsum(ifr)+wstro
          esum1s(ifr)=eold(ifr)*wstro
        endif
  380 continue

      call acusig

      do 390 i=1,nelec
        ekineo(i)=ekinen(i)
        ekin(idist(i))=ekin(idist(i))+ekineo(i)
  390   ekin2(idist(i))=ekin2(idist(i))+ekineo(i)**2

      return
      end

      subroutine sample_efp(inew,x,ener,prob)
c Written by Claudia Filippi

      use atom_mod
      use const_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
      include 'fitefp.h'
!JT      include 'force.h'

      parameter(NFITCX=MEFP_FIT*MCTYPE)

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /efpbasis/ dlrdesc(MEFP_FIT,MCTYPE),dlrfixc(MEFP_FIT,MCTYPE),
     &dlrfixd(MEFP_FIT,MCTYPE),alpha(MCTYPE),rc(MCTYPE),nlefp(MCTYPE),
     &nsrbase(MCTYPE,MEFP_NL),nsrbast(MCTYPE),nsrbasx(MCTYPE),
     &nlrbase(MCTYPE),nlrfix(MCTYPE),nbpoint(MCTYPE),nbftot
      common /efpave/ rhoav(NFITCX),erhoav(NFITCX),rho2av(NFITCX,NFITCX)

      dimension x(3,*),vbase(3*MEFP_FIT),rho(MEFP_FIT),rhold(MEFP_FIT),
     &oefp(MEFP_NL)

      if(inew.eq.1) then

        do 10 ib=1,nbftot
  10      rho(ib)=0.d0

        do 50 ic=1,ncent
          it=iwctype(ic)
          ibc=nbpoint(it)
          do 50 i=1,nelec
            dx=x(1,i)-cent(1,ic)
            dy=x(2,i)-cent(2,ic)
            dz=x(3,i)-cent(3,ic)
            rad=dsqrt(dx*dx+dy*dy+dz*dz)

            call fitbasis(rad,vbase,it)
            iel=i
            if(nlefp(it).gt.0) call nlocefp(iel,ic,dx,dy,dz,rad,oefp)

c local efp potential (non-fixed short and long range)
            nbloc=nsrbase(it,1)+nlrbase(it)
            do 20 ib=1,nbloc
  20          rho(ibc+ib)=rho(ibc+ib)+vbase(ib)

c non-local efp potential
            ib0=ibc+nbloc
            do 40 l=1,nlefp(it)
              do 30 ib=1,nsrbase(it,l+1)
  30            rho(ib0+ib)=rho(ib0+ib)+vbase(ib)*oefp(l)
  40            ib0=ib0+nsrbase(it,l+1)

c fixed long-range local potential
            nbvar=nsrbasx(it)+nlrbase(it)
            do 50 ib=1,nlrfix(it)
  50          rho(ib0+ib)=rho(ib0+ib)+vbase(nbvar+ib)

       else

        do 60 ib=1,nbftot
  60      rho(ib)=rhold(ib)

      endif

      do 70 ib=1,nbftot
          rhoav(ib)=rhoav(ib)+prob*rho(ib)
          erhoav(ib)=erhoav(ib)+prob*ener*rho(ib)
          do 70 jb=1,ib
  70        rho2av(ib,jb)=rho2av(ib,jb)+prob*rho(ib)*rho(jb)

      return

      entry efpsav

      do 80 ib=1,nbftot
  80    rhold(ib)=rho(ib)

      return

      entry startefp

      do 90 ib=1,nbftot
        rhoav(ib)=0.d0
        erhoav(ib)=0.d0
        do 90 jb=1,ib
          rho2av(ib,jb)=0.d0
  90      rho2av(jb,ib)=0.d0

      return

      entry restartefp

      open(11,file='tape11',status='old',form='unformatted')
      read(11) nbftot,etot,passes,(rhoav(ib),erhoav(ib),ib=1,nbftot),
     &((rho2av(ib,jb),jb=1,nbftot),ib=1,nbftot)
      close(11)

      return
      end
c-----------------------------------------------------------------------

      subroutine fitbasis(rad,vbase,it)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'fitefp.h'

      common /efpbasis/ dlrdesc(MEFP_FIT,MCTYPE),dlrfixc(MEFP_FIT,MCTYPE),
     &dlrfixd(MEFP_FIT,MCTYPE),alpha(MCTYPE),rc(MCTYPE),nlefp(MCTYPE),
     &nsrbase(MCTYPE,MEFP_NL),nsrbast(MCTYPE),nsrbasx(MCTYPE),
     &nlrbase(MCTYPE),nlrfix(MCTYPE),nbpoint(MCTYPE),nbftot

      dimension vbase(*)

      do 10 ib=1,nsrbasx(it)
 10     vbase(ib)=bset(ib-1,0.d0,rad,alpha(it),rc(it))

      ib0=nsrbasx(it)
      do 20 ib=1,nlrbase(it)
 20     vbase(ib0+ib)=bset(-1,dlrdesc(ib,it),rad,alpha(it),rc(it))

      ib0=nsrbasx(it)+nlrbase(it)
      do 30 ib=1,nlrfix(it)
 30     vbase(ib0+ib)=bset(-1,dlrfixd(ib,it),rad,alpha(it),rc(it))

      return
      end
c-----------------------------------------------------------------------

      function bset(ibase,dlrd,rad,alpha,rc)

      use const_mod
      implicit real*8(a-h,o-z)
      include 'fitefp.h'

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

c same basis set used inn Xavier's code for vxc potential
c alpha    = coefficient for generating modified pw's
c ibase    = if -1, build a long-range function
c            otherwise, build a short-range function
c            in this latter case, ibase define
c            the argument of the cosine function
c dlrd     = long-range function description parameter
c            if lrdesc negative, exponential func
c            if lrdesc positive, inverse power law
c rc      = value of the large (cut-off) radius

c long-range functions
      if(ibase.eq.-1)then

c inverse power law functions
       if(dlrd.gt.0.0d0)then
c        rad2=rad**2
c        vbase=(rc**2*exp(-rad2/rc**2)+rad2)**(-0.5*dlrd)

         rad2=rad**2
         rad4=rad2*rad2
         vbase=(rc**4*exp(-rad4/rc**4)+rad4)**(-0.25d0*dlrd)
        else
c exponential functions
         zr=-dlrd*rad
         expzr=exp(zr)
         expzrc=exp(-dlrd*rc)
         vbase=1.d0/(expzrc+expzr-zr)
       endif

      else

c short-range functions
c       xx=alpha*rad/(1.d0+alpha*rad)
c distorted cosines
        if(ibase.eq.0) then
          vbase=1.0d0
        else
          vbase=cos(ibase*pi*rad/rc)
c         vbase=cos(ibase*pi*xx)-(-1)**ibase
        endif
        vbase=vbase*exp(-(rad/rc)**4)
c       vbase=vbase*exp(-(rad/rc)**2)
      endif

      bset=vbase

      return
      end
c-----------------------------------------------------------------------

      subroutine readbas

      use atom_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
      include 'fitefp.h'

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /efpbasis/ dlrdesc(MEFP_FIT,MCTYPE),dlrfixc(MEFP_FIT,MCTYPE),
     &dlrfixd(MEFP_FIT,MCTYPE),alpha(MCTYPE),rc(MCTYPE),nlefp(MCTYPE),
     &nsrbase(MCTYPE,MEFP_NL),nsrbast(MCTYPE),nsrbasx(MCTYPE),
     &nlrbase(MCTYPE),nlrfix(MCTYPE),nbpoint(MCTYPE),nbftot
      common /efpqua/ xq0(MEFP_QUAD),yq0(MEFP_QUAD),zq0(MEFP_QUAD),wq(MEFP_QUAD),
     &nqefp

      nqefp=0

      open(32,file='basis.fit',status='old')

      nbftot=0
      nbpoint(1)=0
      do 10 it=1,nctype

c       if(it.gt.1) nbpoint(it)=nbpoint(it-1)+nbfit

        read(32,*) nlefp(it)
        read(32,*) nsrbase(it,1),nlrbase(it),nlrfix(it)

        read(32,*) rc(it),alpha(it)
        read(32,*) (dlrdesc(ib,it),ib=1,nlrbase(it))
        read(32,*) (dlrfixd(ib,it),ib=1,nlrfix(it))
        read(32,*) (dlrfixc(ib,it),ib=1,nlrfix(it))

        nsrbast(it)=nsrbase(it,1)
        nsrbasx(it)=nsrbase(it,1)

        if(nlefp(it).gt.0) then
          read(32,*) (nsrbase(it,l+1),l=1,nlefp(it)),nqefp
          if(nqefp.gt.MEFP_QUAD) stop 'MEFP_QUAD too small'

          do 5 l=1,nlefp(it)
            ll=l+1
            if(nsrbase(it,ll).gt.nsrbasx(it)) nsrbasx(it)=nsrbase(it,ll)
   5        nsrbast(it)=nsrbast(it)+nsrbase(it,ll)
        endif

        nbfit=nsrbast(it)+nlrbase(it)+nlrfix(it)
        if(nbfit.gt.MEFP_FIT) stop 'MEFP_FIT too small'
        if(it.lt.nctype) nbpoint(it+1)=nbpoint(it)+nbfit

        write(6,*) nbfit,nsrbase(it,1),nlrbase(it),nlrfix(it)
        write(6,*) rc(it),alpha(it)
        write(6,*) (dlrdesc(ib,it),ib=1,nlrbase(it))
        write(6,*) (dlrfixd(ib,it),ib=1,nlrfix(it))
        write(6,*) (dlrfixc(ib,it),ib=1,nlrfix(it))

  10    nbftot=nbftot+nbfit

      write(6,*) nbftot
      close(32)

      if(nqefp.gt.0) call gesqua(nqefp,xq0,yq0,zq0,wq)

      return
      end
c-----------------------------------------------------------------------

      subroutine writebas(passes,etot)

      use atom_mod
      use contrl_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
      include 'fitefp.h'
!JT      include 'force.h'

      parameter(NFITCX=MEFP_FIT*MCTYPE)

!JT      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_global,nconf_new,isite,idump,irstar
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /efpbasis/ dlrdesc(MEFP_FIT,MCTYPE),dlrfixc(MEFP_FIT,MCTYPE),
     &dlrfixd(MEFP_FIT,MCTYPE),alpha(MCTYPE),rc(MCTYPE),nlefp(MCTYPE),
     &nsrbase(MCTYPE,MEFP_NL),nsrbast(MCTYPE),nsrbasx(MCTYPE),
     &nlrbase(MCTYPE),nlrfix(MCTYPE),nbpoint(MCTYPE),nbftot
      common /efpave/ rhoav(NFITCX),erhoav(NFITCX),rho2av(NFITCX,NFITCX)

      dimension vbase(MEFP_FIT),y(NFITCX),ipvt(NFITCX),v(MCENT,MEFP_NL)

      if(idump.eq.1) then
        open(11,file='tape11',status='unknown',form='unformatted')
        write(11) nbftot,etot,passes,(rhoav(ib),erhoav(ib),ib=1,nbftot),
     &  ((rho2av(ib,jb),jb=1,nbftot),ib=1,nbftot)
        close(11)
      endif

      efin=etot/passes
      do 10 ib=1,nbftot
         rhoav(ib)=rhoav(ib)/passes
         erhoav(ib)=erhoav(ib)/passes-efin*rhoav(ib)
         do 10 jb=1,ib
           rho2av(ib,jb)=rho2av(ib,jb)/passes-rhoav(ib)*rhoav(jb)
  10       rho2av(jb,ib)=rho2av(ib,jb)

      ib0=0
      do 30 it=1,nctype
        ibc=nbpoint(it)
        nbvar=nsrbast(it)+nlrbase(it)
        do 20 ib=1,nbvar
          y(ib0+ib)=erhoav(ibc+ib)
          do 20 jt=1,nctype
            jbc=nbpoint(jt)+nsrbast(jt)+nlrbase(jt)
            do 20 jb=1,nlrfix(jt)
  20          y(ib0+ib)=y(ib0+ib)-dlrfixc(jb,jt)*rho2av(jbc+jb,ibc+ib)
  30    ib0=ib0+nbvar

      ib0=nsrbast(1)+nlrbase(1)
      do 50 it=2,nctype
        ibc=nbpoint(it)
        nbvari=nsrbast(it)+nlrbase(it)
        jb0=0
        do 45 jt=1,it
          jbc=nbpoint(jt)
          nbvarj=nsrbast(jt)+nlrbase(jt)
          do 40 ib=1,nbvari
            do 40 jb=1,nbvarj
              rho2av(ib0+ib,jb0+jb)=rho2av(ibc+ib,jbc+jb)
  40          rho2av(jb0+jb,ib0+ib)=rho2av(ibc+ib,jbc+jb)
  45      jb0=jb0+nbvarj
  50    ib0=ib0+nbvari

      call dgefa(rho2av,nfitcx,ib0,ipvt,info)
      if(info.ne.0) stop 'problem in writebas'
      job=0
      call dgesl(rho2av,nfitcx,ib0,ipvt,y,job)

      write(6,'(''effective potential coefficients:'')')
      write(6,'(1p10e20.10)') (y(ib),ib=1,ib0)

      do 80 irad=1,2000
        rr=0.004d0*(irad-1)
        ibc=0
        do 75 it=1,nctype

          do 55 l=1,nlefp(it)+1
  55        v(it,l)=0.d0

          call fitbasis(rr,vbase,it)
          nbloc=nsrbase(it,1)+nlrbase(it)
          do 60 ib=1,nbloc
  60        v(it,1)=v(it,1)+y(ibc+ib)*vbase(ib)
          nbvar=nsrbasx(it)+nlrbase(it)
          do 65 ib=1,nlrfix(it)
  65        v(it,1)=v(it,1)+dlrfixc(ib,it)*vbase(ib+nbvar)

          ibc=ibc+nbloc
          do 75 l=1,nlefp(it)
            do 70 ib=1,nsrbase(it,l+1)
  70           v(it,l+1)=v(it,l+1)+y(ibc+ib)*vbase(ib)
  75        ibc=ibc+nsrbase(it,l+1)

  80      write(6,'(1p10e20.10)')rr,((v(i,l),l=1,nlefp(i)+1),i=1,nctype)

      return
      end
c-----------------------------------------------------------------------

      subroutine nlocefp(i,ic,dx,dy,dz,r,oefp)

      use atom_mod
      use dets_mod
      use slater_mod
      use distance_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
      include 'fitefp.h'
!JT      include 'force.h'

!JT      parameter (zero=0.d0,one=1.d0)

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent

      common /efpbasis/ dlrdesc(MEFP_FIT,MCTYPE),dlrfixc(MEFP_FIT,MCTYPE),
     &dlrfixd(MEFP_FIT,MCTYPE),alpha(MCTYPE),rc(MCTYPE),nlefp(MCTYPE),
     &nsrbase(MCTYPE,MEFP_NL),nsrbast(MCTYPE),nsrbasx(MCTYPE),
     &nlrbase(MCTYPE),nlrfix(MCTYPE),nbpoint(MCTYPE),nbftot
      common /efpqua/ xq0(MEFP_QUAD),yq0(MEFP_QUAD),zq0(MEFP_QUAD),wq(MEFP_QUAD),
     &nqefp

!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
!JT      common /slater/ slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD)
!JT     &,fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
!JT     &,fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
!JT     &,detu(MDETUD),detd(MDETUD)
!JT     &,ddeti_deti(3,MELEC,MDETUD),d2edeti_deti(MELEC,MDETUD),deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD),d2deti_det(MPARMD),d2det_det
!JT     &,detij_det(MPARMD,MPARMD)
!JT      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      dimension oefp(*),xn(3)

      it=iwctype(ic)
      do 20 l=1,nlefp(it)
   20   oefp(l)=zero

      ri=one/r
      do 35 iq=1,nqefp
        costh=dx*xq0(iq)+dy*yq0(iq)+dz*zq0(iq)
        costh=costh*ri

        xn(1)=r*xq0(iq)+cent(1,ic)
        xn(2)=r*yq0(iq)+cent(2,ic)
        xn(3)=r*zq0(iq)+cent(3,ic)

c Warning: we probably need to evaluate rvec_en,r_en before making this call
        call nonlocd(i,xn,rvec_en,r_en,detu,detd,slmui,slmdi,deter)

        do 30 l=1,nlefp(it)
   30     oefp(l)=oefp(l)+wq(iq)*yl0(l,costh)*deter

   35   continue

      deter=zero
c     do 40 idet=1,ndet
c  40   deter=deter+cdet(idet,1)*detu(idet)*detd(idet)
      do 40 icsf=1,ncsf
        do 40 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
   40     deter=deter+csf_coef(icsf,1)*cdet_in_csf(idet_in_csf,icsf)*detu(idet)*detd(idet)

      do 50 l=1,nlefp(it)
   50   oefp(l)=oefp(l)/deter

      return
      end

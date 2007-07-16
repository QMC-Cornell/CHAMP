c test basis_fns.f
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'numbas.h'
      include 'force.h'
      parameter(eps=1.d-5)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis2/ n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),zex2(MRWF,MCTYPE,MWF),iwrwf2(MBASIS)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC)
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en_sav(3,MELEC,MCENT),dphin_num(3,MBASIS),d2phin_num(MBASIS)


      numr=0
      nelec=1
      ncent=1
      ndim=3
      betaq=1.d0
      iwctype(1)=1
      iwf=1
      read(5,*) nbasis
      if(nbasis.gt.MBASIS) stop 'nbasis>MBASIS'
      nbasis_ctype(1)=nbasis
      nrbas(1)=nbasis
      read(5,*) (n_bas(i),l_bas(i),m_bas(i),i=1,nbasis)
      read(5,*) (zex2(i,1,1),i=1,nbasis)
      do 5 i=1,nbasis
        iwrwf(i,1)=i
    5   n_bas2(i,1)=n_bas(i)
      read(5,*) (rvec_en(k,1,1),k=1,3)
      r_en(1,1)=0
      do 10 k=1,3
        rvec_en_sav(k,1,1)=rvec_en(k,1,1)
   10   r_en(1,1)=r_en(1,1)+(rvec_en(k,1,1)**2)
      r_en(1,1)=dsqrt(r_en(1,1))

      call basis_fns(1,rvec_en,r_en)

      do 20 ib=1,nbasis
        d2phin_num(ib)=-6*phin(ib,1)/eps**2
   20   write(6,'(i4,9f10.6)') ib,phin(ib,1),(dphin(k,ib,1),k=1,3),d2phin(ib,1)

c calculate derivatives
      do 90 k=1,3
c +eps in k direction
        do 30 kk=1,3
   30     rvec_en(kk,1,1)=rvec_en_sav(kk,1,1)
        rvec_en(k,1,1)=rvec_en_sav(k,1,1)+eps
        r_en(1,1)=0
        do 40 kk=1,3
   40     r_en(1,1)=r_en(1,1)+(rvec_en(kk,1,1)**2)
        r_en(1,1)=dsqrt(r_en(1,1))
        call basis_fns(1,rvec_en,r_en)
        do 50 ib=1,nbasis
          dphin_num(k,ib)=phin(ib,1)
   50   d2phin_num(ib)=d2phin_num(ib)+phin(ib,1)/(eps**2)
c -eps in k direction
        rvec_en(k,1,1)=rvec_en_sav(k,1,1)-eps
        r_en(1,1)=0
        do 60 kk=1,3
   60     r_en(1,1)=r_en(1,1)+(rvec_en(kk,1,1)**2)
        r_en(1,1)=dsqrt(r_en(1,1))
        call basis_fns(1,rvec_en,r_en)
        do 90 ib=1,nbasis
          dphin_num(k,ib)=(dphin_num(k,ib)-phin(ib,1))/(2*eps)
   90     d2phin_num(ib)=d2phin_num(ib)+phin(ib,1)/(eps**2)

      do 100 ib=1,nbasis
        write(6,'(i4,9f11.7)') ib,phin(ib,1),(dphin_num(k,ib),k=1,3),d2phin_num(ib)
        write(6,'(i4,9f11.7)') ib,phin(ib,1),(dphin(k,ib,1),k=1,3),d2phin(ib,1)
        do 95 k=1,3
   95     if(abs((dphin_num(k,ib)-dphin(k,ib,1))/dphin(k,ib,1)).gt.1.d-3)
     &    write(6,'(''***Warning: error in dphin, dimension'',i3,'' basis'',i3,f9.5)') k,ib,dphin_num(k,ib)/dphin(k,ib,1)
  100     if(abs((d2phin_num(ib)-d2phin(ib,1))/d2phin(ib,1)).gt.1.d-3)
     &    write(6,'(''***Warning: error in d2phin, basis'',i3,f9.5)') ib,d2phin_num(ib)/d2phin(ib,1)
 
      stop
      end
c-----------------------------------------------------------------------

      subroutine basis_fns(iel,rvec_en,r_en)
c Written by Cyrus Umrigar
c Calculate 3-dim localised basis functions and their derivatives
c iel = 0, for all electrons
c    != 0, for electron iel

c In input:
c n1s,n2s,...     > 0 : Slater basis
c nsa,npa,nda     = 0 : asymptotic basis
c n1s,nsa,...     < 0 : Gaussian basis
c Here:
c n_bas2(irb,ict) > 0 : Slater basis
c                 = 0 : asymptotic basis
c                 < 0 : Gaussian basis

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'numbas.h'
      include 'force.h'

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis2/ n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),zex2(MRWF,MCTYPE,MWF),iwrwf2(MBASIS)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC)
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

      dimension wfv(3,MRWF),xc(3)
     &,th(0:ML_BAS,0:ML_BAS),dth(3,0:ML_BAS,0:ML_BAS),ph(-ML_BAS:ML_BAS),dph(3,-ML_BAS:ML_BAS)

c Here we have additional normalization factors beyond those in basis_norm, viz., sqrt((2*l+1)/(4*pi))
      data cd0,cd1,cd2,cf0,cf1,cf2,cf3/0.5d0,1.73205080756888d0,0.866025403784439d0,1.d0,1.d0,1.d0,1.d0/

      ider=1

c Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
       else
        nelec1=iel
        nelec2=iel
      endif

      do 40 ie=nelec1,nelec2
        ib=0
        do 40 ic=1,ncent
          ict=iwctype(ic)

          if(r_en(ie,ic).eq.0.d0) then
            r_en(ie,ic)=1.d-99
            rvec_en(1,ie,ic)=1.d-99
          endif
          xc(1)=rvec_en(1,ie,ic)
          xc(2)=rvec_en(2,ie,ic)
          xc(3)=rvec_en(3,ie,ic)
          r=r_en(ie,ic)
          r2=r*r
          ri=1/r
          ri2=ri*ri
          ri3=ri2*ri
          ri4=ri3*ri
          ri5=ri4*ri

          if(numr.eq.0) then

            do 10 irb=1,nrbas(ict)
              n=n_bas2(irb,ict)
              rn=abs(n)
c             rnm2=max(0.d0,rn-2)
              rm3=r**(rn-3)
              rm2=rm3*r
              rm1=rm2*r
c Slater r^(n-1)*Exp(-zeta*r)
              if(n.gt.0) then
                zr=zex2(irb,ict,iwf)*r
                ex=dexp(-zr)
                wfv(3,irb)=rm3*((rn-1)*(rn-2-zr)-(rn-1-zr)*zr)*ex
               elseif(n.eq.0) then
c Warning: Asymptotic and Gaussian not yet tested.
c Asymptotic r^(rn-1)*Exp(-zeta*r), where rn=beta+1, beta=betaq/zeta-1, zeta=sqrt(-2*E_ion)?
                stop 'asymptotic not yet fully tested'
                rn=betaq/zex2(irb,ict,iwf)
                rm3=r**(rn-3)
                rm2=rm3*r
                rm1=rm2*r
                wfv(3,irb)=rm3*((rn-1)*(rn-2-zr)-(rn-1-zr)*zr)*ex
               elseif(n.lt.0) then
c Gaussian  r^(n-1)*Exp(-zeta*r^2) (to be checked)
c               rnm2=rn-2
                zr=2*zex2(irb,ict,iwf)*r2
                ex=dexp(-0.5d0*zr)
                wfv(3,irb)=rm3*((rn-1)*(rn-2)-(2*rn-1-zr)*zr)*ex
              endif
              wfv(1,irb)=rm1*ex
   10         wfv(2,irb)=rm2*((rn-1)-zr)*ex

c       write(6,'(''ic,ict,nrbas(ict),wfv='',3i5,29d12.5)') ic,ict,nrbas(ict),(wfv(1,irb),irb=1,nrbas(ict))

           else

            rk=r
c           do 20 irb=1,nrbas(ict)
c             call splfit(rk,irb,ict,iwf,wfv(1,irb),ider)
c  20         if(wfv(1,irb).eq.0.d0) wfv(1,irb)=DBLMIN

          endif

          do 40 ib2=1,nbasis_ctype(ict)

            xx=xc(1)
            yy=xc(2)
            zz=xc(3)
            xx2=xx*xx
            yy2=yy*yy
            zz2=zz*zz

            xhat=xx*ri
            yhat=yy*ri
            zhat=zz*ri

c Phi function

            ph(0)=1
            dph(1,0)=0
            dph(2,0)=0

            ph(1)=xx
            dph(1,1)=1
            dph(2,1)=0

            ph(-1)=yy
            dph(1,-1)=0
            dph(2,-1)=1

            ph(2)=xx2-yy2
            dph(1,2)= 2*xx
            dph(2,2)=-2*yy

            ph(-2)=2*xx*yy
            dph(1,-2)=2*yy
            dph(2,-2)=2*xx

c           ph(3)=(xx2-yy2)*xx-2*xx*yy2
            ph(3)=ph(2)*ph(1)-ph(-2)*ph(-1)
            dph(1,3)=3*ph(2)
            dph(2,3)=-3*ph(-2)

            ph(-3)=ph(-2)*ph(1)+ph(-1)*ph(2)
            dph(1,-3)=3*ph(-2)
            dph(2,-3)=3*ph(2)

c Theta function

            th(0,0)=1
            dth(1,0,0)=0
            dth(2,0,0)=0
            dth(3,0,0)=0

            th(1,0)=ri*zz
            dth(1,1,0)=-ri3*xx*zz
            dth(2,1,0)=-ri3*yy*zz
            dth(3,1,0)= ri3*(xx2+yy2)

            th(1,1)=ri
            dth(1,1,1)=-ri3*xx
            dth(2,1,1)=-ri3*yy
            dth(3,1,1)=-ri3*zz

c           th(2,0)=cd0*ri3*(3*zz**2-r2)
            th(2,0)=cd0*(3*zhat**2-1)
            dth(1,2,0)=-6*cd0*ri*zhat**2*xhat
            dth(2,2,0)=-6*cd0*ri*zhat**2*yhat
            dth(3,2,0)= 6*cd0*ri4*(xx2+yy2)*zz

            th(2,1)=cd1*ri2*zz
            dth(1,2,1)=-2*cd1*ri4*zz*xx
            dth(2,2,1)=-2*cd1*ri4*zz*yy
            dth(3,2,1)=   cd1*ri4*(xx2+yy2-zz2)

            th(2,2)=cd2*ri2
            dth(1,2,2)=-2*cd2*ri4*xx
            dth(2,2,2)=-2*cd2*ri4*yy
            dth(3,2,2)=-2*cd2*ri4*zz

            th(3,0)=cf0*ri3*zz*(2*zz2-3*(xx2+yy2))
            dth(1,3,0)= 3*cf0*ri5*(xx2+yy2-4*zz2)*zz*xx
            dth(2,3,0)= 3*cf0*ri5*(xx2+yy2-4*zz2)*zz*yy
c           dth(3,3,0)=-3*cf0*ri5*(4*(xx2+yy2)**2-zz2*(7*(xx2+yy2)+6*zz2))
            dth(3,3,0)=-3*cf0*ri5*((xx2+yy2)**2-4*zz2*(xx2+yy2))

            th(3,1)=cf1*ri3*(4*zz2-(xx2+yy2))
            dth(1,3,1)=cf1*ri5*(xx2+yy2-14*zz2)*xx
            dth(2,3,1)=cf1*ri5*(xx2+yy2-14*zz2)*yy
            dth(3,3,1)=cf1*ri5*(11*(xx2+yy2)-4*zz2)*zz

            th(3,2)=cf2*ri3*zz
            dth(1,3,2)=-3*cf2*ri5*zz*xx
            dth(2,3,2)=-3*cf2*ri5*zz*yy
            dth(3,3,2)=   cf2*ri5*(xx2+yy2-2*zz2)

            th(3,3)=cf3*ri3
            dth(1,3,3)=-3*cf3*ri5*xx
            dth(2,3,3)=-3*cf3*ri5*yy
            dth(3,3,3)=-3*cf3*ri5*zz

            ib=ib+1
            irb=iwrwf(ib2,ict)
            l=l_bas(ib)
            m=m_bas(ib)
            mabs=abs(m)
            ylm=th(l,mabs)*ph(m)
            phin(ib,ie)=ylm*wfv(1,irb)
            dphin(1,ib,ie)=ylm*xhat*wfv(2,irb)+(th(l,mabs)*dph(1,m)+dth(1,l,mabs)*ph(m))*wfv(1,irb)
            dphin(2,ib,ie)=ylm*yhat*wfv(2,irb)+(th(l,mabs)*dph(2,m)+dth(2,l,mabs)*ph(m))*wfv(1,irb)
            dphin(3,ib,ie)=ylm*zhat*wfv(2,irb)+(                    dth(3,l,mabs)*ph(m))*wfv(1,irb)
            d2phin(ib,ie)=ylm*(wfv(3,irb)+2*ri*wfv(2,irb)-l*(l+1)*ri2*wfv(1,irb))
c     write(6,'(''ie,ib,ib2,irb,phin(ib,ie),dphin,d2phin'',4i5,9d12.5)')
c    & ie,ib,ib2,irb,phin(ib,ie),(dphin(kk,ib,ie),kk=1,ndim),d2phin(ib,ie),wfv(1,irb)
c    &,wfv(1,irb)
   40 continue

c     do 2000 l=1,nbasis
c2000   write(6,'(''basisfns='',20d12.5)') phin(l,nelec),(dphin(ie,l,nelec),ie=1,ndim),d2phin(l,nelec)

      return
      end
c-----------------------------------------------------------------------

      subroutine basis_fns_2d(iel,rvec_en,r_en)
c Written by Cyrus Umrigar
c Calculate 2-dim localised basis functions and their derivatives for all electrons

c In input:
c n1s,n2s,...     > 0 : Slater basis
c nsa,npa,nda     = 0 : asymptotic basis
c n1s,nsa,...     < 0 : Gaussian basis
c Here:
c n_bas2(irb,ict) > 0 : Slater basis
c                 = 0 : asymptotic basis
c                 < 0 : Gaussian basis

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'numbas.h'
      include 'force.h'

c     common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis2/ n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),zex2(MRWF,MCTYPE,MWF),iwrwf2(MBASIS)
c     common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC)
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

      dimension wfv(3,MRWF),xc(3)
     &,ph(-ML_BAS:ML_BAS),dph(3,-ML_BAS:ML_BAS)

      ider=1

c Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
       else
        nelec1=iel
        nelec2=iel
      endif

      do 40 ie=nelec1,nelec2
        ib=0
        do 40 ic=1,ncent
          ict=iwctype(ic)

          if(r_en(ie,ic).eq.0.d0) then
            r_en(ie,ic)=1.d-99
            rvec_en(1,ie,ic)=1.d-99
          endif
          xc(1)=rvec_en(1,ie,ic)
          xc(2)=rvec_en(2,ie,ic)
          r=r_en(ie,ic)
          r2=r*r
          ri=1/r
          ri2=ri*ri
          ri3=ri2*ri

          if(numr.eq.0) then

            do 10 irb=1,nrbas(ict)
              n=n_bas2(irb,ict)
              rn=abs(n)
c             rnm2=max(0.d0,rn-2)
              rm3=r**(rn-3)
              rm2=rm3*r
              rm1=rm2*r
c Slater r^(n-1)*Exp(-zeta*r)
              if(n.gt.0) then
                zr=zex2(irb,ict,iwf)*r
                ex=dexp(-zr)
                wfv(3,irb)=rm3*((rn-1)*(rn-2-zr)-(rn-1-zr)*zr)*ex
               elseif(n.eq.0) then
c Warning: Asymptotic and Gaussian not yet tested.
c Asymptotic r^(rn-1)*Exp(-zeta*r), where rn=beta+1, beta=betaq/zeta-1, zeta=sqrt(-2*E_ion)?
                stop 'asymptotic not yet fully tested'
                rn=betaq/zex2(irb,ict,iwf)
                rm3=r**(rn-3)
                rm2=rm3*r
                rm1=rm2*r
                wfv(3,irb)=rm3*((rn-1)*(rn-2-zr)-(rn-1-zr)*zr)*ex
               elseif(n.lt.0) then
c Gaussian  r^(n-1)*Exp(-zeta*r^2) (to be checked)
c               rnm2=rn-2
                zr=2*zex2(irb,ict,iwf)*r2
                ex=dexp(-0.5d0*zr)
                wfv(3,irb)=rm3*((rn-1)*(rn-2)-(2*rn-1-zr)*zr)*ex
              endif
              wfv(1,irb)=rm1*ex
   10         wfv(2,irb)=rm2*((rn-1)-zr)*ex

c       write(6,'(''ic,ict,nrbas(ict),wfv='',3i5,29d12.5)') ic,ict,nrbas(ict),(wfv(1,irb),irb=1,nrbas(ict))

           else

            rk=r
c           do 20 irb=1,nrbas(ict)
c  20         call splfit(rk,irb,ict,iwf,wfv(1,irb),ider)

          endif

          do 40 ib2=1,nbasis_ctype(ict)

            xx=xc(1)
            yy=xc(2)
            xhat=xx*ri
            yhat=yy*ri

c Phi function

            ph(0)=1
            dph(1,0)=0
            dph(2,0)=0
            ph(1)=xhat
            ph(-1)=yhat
            dph(1,1)=yhat**2*ri
            dph(2,1)=-xhat*yhat*ri
            dph(1,-1)=-xhat*yhat*ri
            dph(2,-1)=xhat**2*ri
            do 30 i=2,ML_BAS
              ph(i)=ph(i-1)*ph(1)-ph(1-i)*ph(-1)
              ph(-i)=ph(1-i)*ph(1)+ph(i-1)*ph(-1)
              dph(1,i) = i*ph(-i)*ph(-1)*ri
              dph(2,i) =-i*ph(-i)*ph(1)*ri
              dph(1,-i)=-i*ph(i)*ph(-1)*ri
c             write(6,'(''xx,yy,ri,xhat,yhat,ph(i),ph(-i)='',19f9.5)')
c    &        xx,yy,ri,xhat,yhat,ph(i),ph(-i),dph(1,i),dph(2,i),dph(1,-i)
   30         dph(2,-i)= i*ph(i)*ph(1)*ri

            ib=ib+1
            irb=iwrwf(ib2,ict)
            m=m_bas(ib)
            ylm=ph(m)
            phin(ib,ie)=ylm*wfv(1,irb)
            dphin(1,ib,ie)=ylm*xhat*wfv(2,irb)+dph(1,m)*wfv(1,irb)
            dphin(2,ib,ie)=ylm*yhat*wfv(2,irb)+dph(2,m)*wfv(1,irb)
            d2phin(ib,ie)=ylm*(wfv(3,irb)+ri*wfv(2,irb)-m*m*ri2*wfv(1,irb))
c     write(6,'(''mbas'',9i5)') m
c     write(6,'(''ie,ib,ph(m),(dph(k,m),ddph(m)'',2i5,9d12.5)')
c    & ie,ib,ph(m),(dph(k,m),k=1,ndim),-m*m*ri2*wfv(1,irb),wfv(1,irb),wfv(2,irb),wfv(3,irb),ylm*(wfv(3,irb)+ri*wfv(2,irb)),ylm*(-m*m*ri2*wfv(1,irb))
c     write(6,'(''ie,ib,ib2,irb,phin(ib,ie),dphin,d2phin'',4i5,9d12.5)')
c    & ie,ib,ib2,irb,phin(ib,ie),(dphin(kk,ib,ie),kk=1,ndim),d2phin(ib,ie)
   40 continue

c     do 2000 l=1,nbasis
c2000   write(6,'(''basisfns='',20d12.5)') phin(l,nelec),(dphin(ie,l,nelec),ie=1,ndim),d2phin(l,nelec)

      return
      end

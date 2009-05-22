      subroutine cbasis_fns_num(iel,rvec_en,r_en)
c written by Amit Ghosal starting from basis_fns.f

      use atom_mod
      use basis1_mod
      use const_mod
      use numbas_mod
      use basis2_mod
      use wfsec_mod
      implicit real*8(a-h,o-z)
!JT	include 'vmc.h'
!JT	include 'pseudo.h'
!JT	include 'numbas.h'
!JT	include 'force.h'

c New temporary variables defined ************************
      complex*16 temp_p,temp_n,cph,cdph,cphin,cdphin,cd2phin,ylm
c Definition of new temporary variables ends here ********

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere
!JT      common /basis/ zex(MBASIS,MWF),betaq
!JT     &,n1s(MCTYPE)
!JT     &,n2s(MCTYPE),n2p(-1:1,MCTYPE)
!JT     &,n3s(MCTYPE),n3p(-1:1,MCTYPE),n3d(-2:2,MCTYPE)
!JT     &,n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE),n4f(-3:3,MCTYPE)
!JT     &,n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE),n5f(-3:3,MCTYPE)
!JT     &,n5g(-4:4,MCTYPE)
!JT     &,n6d(-2:2,MCTYPE),n6f(-3:3,MCTYPE),n6g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
!JT     &,n7g(-4:4,MCTYPE),n7h(-5:5,MCTYPE),n7i(-6:6,MCTYPE)
!JT     &,n8i(-6:6,MCTYPE),n8j(-7:7,MCTYPE)
!JT     &,n9k(-8:8,MCTYPE)
!JT     &,n10l(-9:9,MCTYPE)
!JT     &,n11m(-10:10,MCTYPE)
!JT     &,n12n(-11:11,MCTYPE)
!JT     &,n13o(-12:12,MCTYPE)
!JT     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
!JT      common /basis2/ zex2(MRWF,MCTYPE,MWF),n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
!JT     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
!JT     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),iwrwf2(MBASIS)
c     common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
c     &,d2phin(MBASIS,MELEC)
!JT      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
!JT     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
!JT     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

!JT      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /cphifun/ cphin(MBASIS,MELEC),cdphin(3,MBASIS,MELEC,MELEC)
     &,cd2phin(MBASIS,MELEC)

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

      dimension wfv(3,MRWF),xc(3)
     &,cph(-ML_BAS:ML_BAS),cdph(3,-ML_BAS:ML_BAS)

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
c Asymptotic r^(rn-1)*Exp(-zeta*r), where rn=beta+1,
c beta=betaq/zeta-1, zeta=sqrt(-2*E_ion)?
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

        write(6,'(''ic,ict,nrbas(ict),wfv='',3i5,29d12.5)') ic,ict,nrbas(ict),(wfv(1,irb),irb=1,nrbas(ict))

           else

            rk=r
            do 20 irb=1,nrbas(ict)
ccc   20         call splfit_bas(rk,irb,ict,iwf,wfv(1,irb),ider)
              call splfit_bas(rk,irb,ict,iwf,wfv(1,irb),ider)
   20         if(wfv(1,irb).eq.0.d0) wfv(1,irb)=DBLMIN

          endif

          do 40 ib2=1,nbasis_ctype(ict)

            xx=xc(1)
            yy=xc(2)
            xhat=xx*ri
            yhat=yy*ri

c Phi function
c ****************** The following is the modified part ********************
            theta=datan2(yy,xx)
            cph(0)=dcmplx(1.d0,0.d0)
            cdph(1,0)=dcmplx(0.d0,0.d0)
            cdph(2,0)=dcmplx(0.d0,0.d0)
            do 30 i=1,ML_BAS
              th_i=i*theta
              cosith=dcos(th_i)
              sinith=dsin(th_i)
              cph(i)=dcmplx(cosith,sinith)
              cph(-i)=dcmplx(cosith,-sinith)
              temp_p=dcmplx(-sinith,cosith)
              temp_n=dcmplx(-sinith,-cosith)
              cdph(1,i)=-i*yhat*ri*temp_p
              cdph(2,i)=i*xhat*ri*temp_p
              cdph(1,-i)=-i*yhat*ri*temp_n
   30         cdph(2,-i)=i*xhat*ri*temp_n
c ****************** The modified part ends here ***************************

            ib=ib+1
            irb=iwrwf(ib2,ict)
            m=m_bas(ib)
            ylm=cph(m)
c	    if(abs(m).gt.ML_BAS) then
c	      write(*,*) 'm=',m
c	      stop
c	    endif
	    cphin(ib,ie)=ylm*wfv(1,irb)
            cdphin(1,ib,ie,1)=ylm*xhat*wfv(2,irb)+cdph(1,m)*wfv(1,irb)
            cdphin(2,ib,ie,1)=ylm*yhat*wfv(2,irb)+cdph(2,m)*wfv(1,irb)
            cd2phin(ib,ie)=ylm*(wfv(3,irb)+ri*wfv(2,irb)-m*m*ri2*wfv(1,irb))
   40 continue

      return
      end

c------------------------------------------------------------------------

      subroutine cbasis_fns_fd(iel,rvec_en,r_en)

c Written by A.D.Guclu, Feb 2004
c to replace basis_fns for quantum dots.
c Calculate Fock-Darwin basis functions and their derivatives for all
c or 1 electron

c arguments: iel=0 -> all electron
c               >0 -> only electron iel
c            rvec_en=vector electron-nucleus
c                    (or electron-dot center in this context)

c output: complex functions cphin,dcphin, and d2cphin are calculated

c Fock-Darwin wavefunctions are:
c phi = norm * r^abs(m) L_n^abs(m)(we r^2) exp(-0.5 we r^2) exp(-imtheta)
c (norm is taken care before in basis_norm_dot)

c For the first derivative, use chain rule.
c Also remember: d/dx L_n^k(x) =-L_(n-1)^(k+1)(x)

c Laplacian can be directly obtained from Shrodinger equation:
c -1/2 d2cphin + 1/2 we^2 r^2 phin = (2n+1+abs(m)) we phin
c (part due to angular momentum excess cancels out )

c In this version we replace we by zex*we for basis set optimization.

      use coefs_mod
      use basis1_mod
      use const_mod
      use wfsec_mod
      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'
!JT      include 'force.h'

      complex*16 ctemp1,ctemp2
      complex*16 cphin,cdphin,cd2phin

c     common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere
!JT       common /basis/ zex(MBASIS,MWF),betaq
!JT     &,n1s(MCTYPE)
!JT     &,n2s(MCTYPE),n2p(-1:1,MCTYPE)
!JT     &,n3s(MCTYPE),n3p(-1:1,MCTYPE),n3d(-2:2,MCTYPE)
!JT     &,n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE),n4f(-3:3,MCTYPE)
!JT     &,n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE),n5f(-3:3,MCTYPE)
!JT     &,n5g(-4:4,MCTYPE)
!JT     &,n6d(-2:2,MCTYPE),n6f(-3:3,MCTYPE),n6g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
!JT     &,n7g(-4:4,MCTYPE),n7h(-5:5,MCTYPE),n7i(-6:6,MCTYPE)
!JT     &,n8i(-6:6,MCTYPE),n8j(-7:7,MCTYPE)
!JT     &,n9k(-8:8,MCTYPE)
!JT     &,n10l(-9:9,MCTYPE)
!JT     &,n11m(-10:10,MCTYPE)
!JT     &,n12n(-11:11,MCTYPE)
!JT     &,n13o(-12:12,MCTYPE)
!JT     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis3/ n_fd(MBASIS),m_fd(MBASIS),n_cf(MBASIS),ncfmax
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /cphifun/ cphin(MBASIS,MELEC),cdphin(3,MBASIS,MELEC,MELEC)
     &,cd2phin(MBASIS,MELEC)
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
!JT      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

c Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
      else
        nelec1=iel
        nelec2=iel
      endif

      ic=1

      do 100 ie=nelec1,nelec2
        if(r_en(ie,ic).eq.0.d0) then
           r_en(ie,ic)=1.d-99
           rvec_en(1,ie,ic)=1.d-99
        endif
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)
        r=r_en(ie,ic)
        r2=r*r
        theta=datan2(x2,x1)

        do 90 ib=1,nbasis
c calculate all basis dependent functions to avoid repeating:
	  wez=we*zex(ib,iwf)
	  rho=r2*wez
	  fexp=dexp(-0.5d0*rho)
          n=n_fd(ib)
          m=m_fd(ib)
          mabs=abs(m)
          thetam=m*theta
          cosm=dcos(thetam)
          sinm=dsin(thetam)
          flagr=flaguerre(n,mabs,rho)
          pow=r**mabs
          pow2=pow/r2
          fac_phi=fexp*pow*flagr
          fac_dphi=fexp*pow2

c calculate cphin:
          cphin(ib,ie)=dcmplx(cosm,-sinm)
          cphin(ib,ie)=cphin(ib,ie)*fac_phi

c calculate cdphin:
          temp=flagr*(mabs-rho)
c following "if" is necessary
          if(n.gt.0) temp=temp-2*rho*flaguerre(n-1,mabs+1,rho)
          ctemp1=dcmplx(temp*cosm,-temp*sinm)

          temp=m*flagr
          ctemp2=dcmplx(temp*sinm,temp*cosm)

          cdphin(1,ib,ie,1)=(x1*ctemp1+x2*ctemp2)*fac_dphi
          cdphin(2,ib,ie,1)=(x2*ctemp1-x1*ctemp2)*fac_dphi

c calculate cd2phin:
          cd2phin(ib,ie)=wez*(rho-2*(2*n+mabs+1))*cphin(ib,ie)

 90     enddo
 100  enddo

      return
      end

c---------------------------------------------------------------------

      function flaguerre(n,m,x)
c calculate laguerre function using recursion relation.
c n and m must be positive real

      implicit real*8(a-h,o-z)

      p2=0.d0
      p1=1.d0
      do 10 j=0,n-1
        p3=p2
        p2=p1
        p1=((2*j+m+1-x)*p2-(j+m)*p3)/(j+1)
 10   enddo

      flaguerre=p1

      return
      end

c----------------------------------------------------------------------

      subroutine cbasis_fns_cf(rvec_en)

c Written by A.D.Guclu, sep 2004
c inspired from the code by Gun Sang Jeon & J.Jain

c Calculates projected composite fermion basis functions and their derivatives for all
c electron case.

c arguments:
c            rvec_en=vector electron-nucleus
c                    (or electron-dot center in this context)

c output: complex functions cphin, and dcphin

c projected cf functions are:
c phi(j) = norm * z^(m+n_landau) * (d^n/dz_j^n \Pi_k (z_j-z_k)^p) / \Pi_k (z_j-z_k)^p
c (norm is taken care before in basis_norm_dot)

c local variables: djnkj = d^n/dz_j^m K_j
c                 didjnkj = d/dz_i d^n/dz_j^m K_j
c                  where e^K = \Pi_k (z_j-z_k)^p)
c                 djnflux = (d^n/dz_j^n \Pi_k (z_j-z_k)^p) / \Pi_k (z_j-z_k)^p
c                 didjnflux = d/dz_i djnflux
c Laplacian is zero and ignored here

c note: one complication here is that the orbitals depend on all electrons...

      use coefs_mod
      use const_mod
      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'
!JT      include 'force.h'


      complex*16 zj,zk,zjk,djnflux,didjnflux
      complex*16 djnkj(MBASIS+1)
      complex*16 didjnkj(MBASIS+1,MELEC)
      complex*16 df(0:MBASIS+1),ddf(0:MBASIS+1,MELEC)
      complex*16 ctemp1,ctemp2,zjpowmn1,zjpowmn

      complex*16 cphin,cdphin,cd2phin

c     common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /basis3/ n_fd(MBASIS),m_fd(MBASIS),n_cf(MBASIS),ncfmax
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /cphifun/ cphin(MBASIS,MELEC),cdphin(3,MBASIS,MELEC,MELEC)
     &,cd2phin(MBASIS,MELEC)
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot

      dimension rvec_en(3,MELEC,MCENT)

      ic=1
      do je=1,nelec
        zj=dcmplx(rvec_en(1,je,ic),-rvec_en(2,je,ic))
        do in=1,ncfmax+1
          djnkj(in)=dcmplx(0.d0,0.d0)
          didjnkj(in,je)=dcmplx(0.d0,0.d0)    ! undefined
        enddo

c calculate djnkj and didjnkj in one go:
        do ke=1,nelec
          if(ke.ne.je) then
            zk=dcmplx(rvec_en(1,ke,ic),-rvec_en(2,ke,ic))
            zjk=zj-zk
c   ctemp1 should be precaculated for faster processing...
            ctemp1=dcmplx(-1.d0,0.d0)/zjk
            ctemp2=dcmplx(1.d0,0.d0)
            do in=1,ncfmax+1
              ctemp2=ctemp1*ctemp2
              djnkj(in)=djnkj(in)-ctemp2
              didjnkj(in,ke)=ctemp2*ctemp1
            enddo
          endif
        enddo
c normalize:
        djnkj(1)=djnkj(1)*nv
        tempfac=1.d0            ! integer but can be big?
        do in=1,ncfmax
          tempfac=tempfac*in
          djnkj(in+1)=nv*tempfac*djnkj(in+1)
          do ie=1,nelec
            didjnkj(in,ie)=nv*tempfac*didjnkj(in,ie)   ! we dont need for in=ncfmax+1
          enddo
        enddo

c now get derivatives of the flux term for all landau levels.
c they are calculated a part because they only depend on landau level index.
c        write(*,*) ncfmax
        do in=0,ncfmax
          df(in)=djnflux(in,djnkj)
          do ie=1,nelec
            ddf(in,ie)=didjnflux(in,ie,je,djnkj,didjnkj)
          enddo
        enddo

c loop on all basis functions
        do ib=1,nbasis
          n=n_cf(ib)
          m=m_fd(ib)
          mn=m+n
          zjpowmn1=zj**(mn-1)      ! should we care about zj near zero?
          zjpowmn=zjpowmn1*zj

c calculate orbital:
          cphin(ib,je)=zjpowmn*df(n)
c          write(*,*) 'je,n,m,zj,cphin',je,n,m,zj,cphin(ib,je)
c calculate gradient of the orbital with respect to the electron ie's coord.:
          do ie=1,nelec
            if(ie.eq.je) then
              cdphin(1,ib,je,ie)=mn*zjpowmn1*df(n)+zjpowmn*ddf(n,ie)
              cdphin(2,ib,je,ie)=cdphin(1,ib,je,ie)*dcmplx(0.d0,-1.d0)
            else
              cdphin(1,ib,je,ie)=zjpowmn*ddf(n,ie)
              cdphin(2,ib,je,ie)=cdphin(1,ib,je,ie)*dcmplx(0.d0,-1.d0)
            endif
          enddo
        enddo

      enddo

      return
      end

c-----------------------------------------------------------------------
      function djnflux(in,dk)
c written by A.D.Guclu, sep 2004
c following Gun Sang Jeon & J.Jain's code.

c calculates (d^n/dj^n e^k)/e^k for given d^n/dj^n k

c code not well optimized. the powers of dk's should be
c precalculated for faster processing.

      include 'vmc.h'

      complex*16 djnflux

c arguments:
      integer in
      complex*16 dk(MBASIS+1)

c locals:
      complex*16 df

      if(in.eq.0) then

        df=dcmplx(1.d0,0.d0)

      elseif(in.eq.1) then

        df=dk(1)

      elseif(in.eq.2) then

        df=dk(1)*dk(1) + dk(2)

      elseif(in.eq.3) then

        df=dk(1)*(dk(1)*dk(1)+3*dk(2))+dk(3)

      elseif(in.eq.4) then

        df=dk(1)**4
        df=df+6*dk(1)*dk(1)*dk(2)
        df=df+3*dk(2)*dk(2)
        df=df+4*dk(1)*dk(3)
        df=df+dk(4)

      elseif(in.eq.5) then

        df=dk(1)**5
        df=df+10*dk(1)**3*dk(2)
        df=df+10*dk(1)*dk(1)*dk(3)
        df=df+10*dk(2)*dk(3)
        df=df+5*dk(1)*(3*dk(2)*dk(2)+dk(4))+dk(5)

      elseif(in.eq.6) then

        df=dk(1)**6+15*dk(1)**4*dk(2)+15*dk(2)**3
        df=df+20*dk(1)**3*dk(3)+10*dk(3)**2
        df=df+15*dk(2)*dk(4)+15*dk(1)**2*(3*dk(2)**2+dk(4))
        df=df+6*dk(1)*(10*dk(2)*dk(3)+dk(5))+dk(6)

      elseif(in.eq.7) then

        df=dk(1)**7+21*dk(1)**5*dk(2)+35*dk(1)**4*dk(3)
        df=df+105*dk(2)**2*dk(3)+35*dk(3)*dk(4)
        df=df+35*dk(1)**3*(3*dk(2)**2+dk(4))
        df=df+21*dk(2)*dk(5)+21*dk(1)**2*(10*dk(2)*dk(3)+dk(5))
        df=df+7*dk(1)*(15*dk(2)**3+10*dk(3)**2+15*dk(2)*dk(4)+dk(6))
        df=df+dk(7)

      elseif(in.eq.8) then

        df=dk(1)**8+28*dk(1)**6*dk(2)+105*dk(2)**4+56*dk(1)**5*dk(3)
        df=df+210*dk(2)**2*dk(4)+35*dk(4)**2+70*dk(1)**4*(3*dk(2)**2+dk(4))
        df=df+56*dk(3)*dk(5)+56*dk(1)**3*(10*dk(2)*dk(3)+dk(5))
        df=df+28*dk(2)*(10*dk(3)**2+dk(6))
        df=df+28*dk(1)**2*(15*dk(2)**3+10*dk(3)**2+15*dk(2)*dk(4)+dk(6))
        df=df+8*dk(1)*(105*dk(2)**2*dk(3)+35*dk(3)*dk(4)+21*dk(2)*dk(5)+dk(7))
        df=df+dk(8)

      else

        stop 'nlandau bigger than 8 in cbasis_fns.f!'

      endif

      djnflux=df

      return
      end

c---------------------------------------------------------------------

      function didjnflux(in,ie,je,dk,ddk)
c written by A.D.Guclu, sep 2004

c calculates d/di ((d^n/dj^n e^k)/e^k) for given d/di d^n/dj^n k and d^n/dj^n k

c far from being optimized. too risky!

      include 'vmc.h'

      complex*16 didjnflux

c arguments:
      integer in,ie,je
      complex*16 dk(MBASIS+1),ddk(MBASIS+1,MELEC)

c locals:
      complex*16 ddf



      if(in.eq.0) then

        ddf=dcmplx(0.d0,0.d0)

      elseif(in.eq.1) then

        if(ie.eq.je) then
          ddf=dk(2)
        else
          ddf=ddk(1,ie)
        endif

      elseif(in.eq.2) then

        if(ie.eq.je) then
          ddf=2*dk(1)*dk(2)+dk(3)
        else
          ddf=2*dk(1)*ddk(1,ie)+ddk(2,ie)
        endif

      elseif(in.eq.3) then

        if(ie.eq.je) then
          ddf=3*(dk(1)*dk(1)*dk(2)+dk(1)*dk(3)+dk(2)*dk(2))+dk(4)
        else
          ddf=3*(ddk(1,ie)*dk(2)+dk(1)*ddk(2,ie)+dk(1)*dk(1)*ddk(1,ie))+ddk(3,ie)
        endif

      elseif(in.eq.4) then

        if(ie.eq.je) then
          ddf=4*dk(1)**3*dk(2)
          ddf=ddf+6*(2*dk(1)*dk(2)*dk(2)+dk(1)*dk(1)*dk(3))
          ddf=ddf+6*dk(2)*dk(3)
          ddf=ddf+4*(dk(2)*dk(3)+dk(1)*dk(4))
          ddf=ddf+dk(5)
        else
          ddf=4*dk(1)**3*ddk(1,ie)
          ddf=ddf+6*(2*dk(1)*ddk(1,ie)*dk(2)+dk(1)*dk(1)*ddk(2,ie))
          ddf=ddf+6*dk(2)*ddk(2,ie)
          ddf=ddf+4*(ddk(1,ie)*dk(3)+dk(1)*ddk(3,ie))
          ddf=ddf+ddk(4,ie)
        endif

      elseif(in.eq.5) then

        if(ie.eq.je) then
          ddf=5*dk(1)**4*dk(2)
          ddf=ddf+10*(3*dk(1)**2*dk(2)**2+dk(1)**3*dk(3))
          ddf=ddf+10*(2*dk(1)*dk(2)*dk(3)+dk(1)**2*dk(4))
          ddf=ddf+10*(dk(3)*dk(3)+dk(2)*dk(4))
          ddf=ddf+15*(dk(2)**3+2*dk(1)*dk(2)*dk(3))
          ddf=ddf+5*(dk(2)*dk(4)+dk(1)*dk(5))
          ddf=ddf+dk(6)
        else
          ddf=5*dk(1)**4*ddk(1,ie)
          ddf=ddf+10*(3*dk(1)**2*ddk(1,ie)*dk(2)+dk(1)**3*ddk(2,ie))
          ddf=ddf+10*(2*dk(1)*ddk(1,ie)*dk(3)+dk(1)**2*ddk(3,ie))
          ddf=ddf+10*(ddk(2,ie)*dk(3)+dk(2)*ddk(3,ie))
          ddf=ddf+15*(ddk(1,ie)*dk(2)*dk(2)+2*dk(1)*dk(2)*ddk(2,ie))
          ddf=ddf+5*(ddk(1,ie)*dk(4)+dk(1)*ddk(4,ie))
          ddf=ddf+ddk(5,ie)
        endif

      elseif(in.eq.6) then

        if(ie.eq.je) then
	  ddf=6*dk(1)**5*dk(2)+60*dk(1)**3*dk(2)**2+15*dk(1)**4*dk(3)+45*dk(2)**2*dk(3)
	  ddf=ddf+60*dk(1)**2*dk(2)*dk(3)+20*dk(1)**3*dk(4)+20*dk(3)*dk(4)
	  ddf=ddf+15*(dk(3)*dk(4)+dk(2)*dk(5))
	  ddf=ddf+90*(dk(1)*dk(2)**3+dk(1)**2*dk(2)*dk(3))
	  ddf=ddf+15*(2*dk(1)*dk(2)*dk(4)+dk(1)**2*dk(5))
	  ddf=ddf+60*(dk(2)**2*dk(3)+dk(1)*dk(3)**2+dk(1)*dk(2)*dk(4))
	  ddf=ddf+6*(dk(2)*dk(5)+dk(1)*dk(6))+dk(7)
        else
	  ddf=6*dk(1)**5*ddk(1,ie)
	  ddf=ddf+15*(4*dk(1)**3*ddk(1,ie)*dk(2)+dk(1)**4*ddk(2,ie)+3*dk(2)**2*ddk(2,ie))
	  ddf=ddf+20*(3*dk(1)**2*ddk(1,ie)*dk(3)+dk(1)**3*ddk(3,ie)+dk(3)*ddk(3,ie))
	  ddf=ddf+15*(ddk(2,ie)*dk(4)+dk(2)*ddk(4,ie))
	  ddf=ddf+90*(dk(1)*ddk(1,ie)*dk(2)**2+dk(1)**2*dk(2)*ddk(2,ie))
	  ddf=ddf+15*(2*dk(1)*ddk(1,ie)*dk(4)+dk(1)**2*ddk(4,ie))
	  ddf=ddf+60*(ddk(1,ie)*dk(2)*dk(3)+dk(1)*ddk(2,ie)*dk(3)+dk(1)*dk(2)*ddk(3,ie))
	  ddf=ddf+6*(ddk(1,ie)*dk(5)+dk(1)*ddk(5,ie))+ddk(6,ie)
        endif

      elseif(in.eq.7) then

        if(ie.eq.je) then
	  ddf=7*dk(1)**6*dk(2)
	  ddf=ddf+21*(5*dk(1)**4*dk(2)*dk(2)+dk(1)**5*dk(3))
	  ddf=ddf+35*(4*dk(1)**3*dk(2)*dk(3)+dk(1)**4*dk(4))
	  ddf=ddf+105*(2*dk(2)*dk(3)*dk(3)+dk(2)**2*dk(4))
	  ddf=ddf+35*(dk(4)*dk(4)+dk(3)*dk(5))
	  ddf=ddf+105*(3*dk(1)**2*dk(2)*dk(2)**2+2*dk(1)**3*dk(2)*dk(3))
	  ddf=ddf+35*(3*dk(1)**2*dk(2)*dk(4)+dk(1)**3*dk(5))
	  ddf=ddf+21*(dk(3)*dk(5)+dk(2)*dk(6))
	  ddf=ddf+210*(2*dk(1)*dk(2)*dk(2)*dk(3)+dk(1)**2*dk(3)*dk(3)+dk(1)**2*dk(2)*dk(4))
	  ddf=ddf+21*(2*dk(1)*dk(2)*dk(5)+dk(1)**2*dk(6))
	  ddf=ddf+105*(dk(2)*dk(2)**3+3*dk(1)*dk(2)**2*dk(3))
	  ddf=ddf+70*(dk(2)*dk(3)**2+2*dk(1)*dk(3)*dk(4))
	  ddf=ddf+105*(dk(2)*dk(2)*dk(4)+dk(1)*dk(3)*dk(4)+dk(1)*dk(2)*dk(5))
	  ddf=ddf+7*(dk(2)*dk(6)+dk(1)*dk(7))+dk(8)
        else
	  ddf=7*dk(1)**6*ddk(1,ie)
	  ddf=ddf+21*(5*dk(1)**4*ddk(1,ie)*dk(2)+dk(1)**5*ddk(2,ie))
	  ddf=ddf+35*(4*dk(1)**3*ddk(1,ie)*dk(3)+dk(1)**4*ddk(3,ie))
	  ddf=ddf+105*(2*dk(2)*ddk(2,ie)*dk(3)+dk(2)**2*ddk(3,ie))
	  ddf=ddf+35*(ddk(3,ie)*dk(4)+dk(3)*ddk(4,ie))
	  ddf=ddf+105*(3*dk(1)**2*ddk(1,ie)*dk(2)**2+2*dk(1)**3*dk(2)*ddk(2,ie))
	  ddf=ddf+35*(3*dk(1)**2*ddk(1,ie)*dk(4)+dk(1)**3*ddk(4,ie))
	  ddf=ddf+21*(ddk(2,ie)*dk(5)+dk(2)*ddk(5,ie))
	  ddf=ddf+210*(2*dk(1)*ddk(1,ie)*dk(2)*dk(3)+dk(1)**2*ddk(2,ie)*dk(3)+dk(1)**2*dk(2)*ddk(3,ie))
	  ddf=ddf+21*(2*dk(1)*ddk(1,ie)*dk(5)+dk(1)**2*ddk(5,ie))
	  ddf=ddf+105*(ddk(1,ie)*dk(2)**3+3*dk(1)*dk(2)**2*ddk(2,ie))
	  ddf=ddf+70*(ddk(1,ie)*dk(3)**2+2*dk(1)*dk(3)*ddk(3,ie))
	  ddf=ddf+105*(ddk(1,ie)*dk(2)*dk(4)+dk(1)*ddk(2,ie)*dk(4)+dk(1)*dk(2)*ddk(4,ie))
	  ddf=ddf+7*(ddk(1,ie)*dk(6)+dk(1)*ddk(6,ie))+ddk(7,ie)
        endif

      elseif(in.eq.8) then

        if(ie.eq.je) then
          ddf=8*dk(1)**7*dk(2)
          ddf=ddf+28*(6*dk(1)**5*dk(2)*dk(2)+dk(1)**6*dk(3))
          ddf=ddf+105*(4*dk(2)**3*dk(3))
          ddf=ddf+56*(5*dk(1)**4*dk(2)*dk(3)+dk(1)**5*dk(4))
          ddf=ddf+210*(2*dk(2)*dk(3)*dk(4)+dk(2)**2*dk(5))
          ddf=ddf+70*dk(4)*dk(5)
          ddf=ddf+70*(4*dk(1)**3*dk(2)*(3*dk(2)**2+dk(4))+dk(1)**4*(6*dk(2)*dk(3)+dk(5)))
          ddf=ddf+56*(dk(4)*dk(5)+dk(3)*dk(6))
          ddf=ddf+56*(3*dk(1)**2*dk(2)*(10*dk(2)*dk(3)+dk(5))+dk(1)**3*(10*(dk(3)*dk(3)+dk(2)*dk(4))+dk(6)))
          ddf=ddf+28*(dk(3)*(10*dk(3)**2+dk(6))+dk(2)*(20*dk(3)*dk(4)+dk(7)))
          ddf=ddf+28*(2*dk(1)*dk(2)*(15*dk(2)**3+10*dk(3)**2+15*dk(2)*dk(4)+dk(6))+
     &                dk(1)**2*(45*dk(2)**2*dk(3)+20*dk(3)*dk(4)+15*(dk(3)*dk(4)+dk(2)*dk(5))+dk(7)))
          ddf=ddf+8*(dk(2)*(105*dk(2)**2*dk(3)+35*dk(3)*dk(4)+21*dk(2)*dk(5)+dk(7))+
     &               dk(1)*(105*(2*dk(2)*dk(3)*dk(3)+dk(2)**2*dk(4))+35*(dk(4)**2+dk(3)*dk(5))+
     &                      21*(dk(3)*dk(5)+dk(2)*dk(6))+dk(8)))
          ddf=ddf+dk(9)
        else
          ddf=8*dk(1)**7*ddk(1,ie)
          ddf=ddf+28*(6*dk(1)**5*ddk(1,ie)*dk(2)+dk(1)**6*ddk(2,ie))
          ddf=ddf+105*(4*dk(2)**3*ddk(2,ie))
          ddf=ddf+56*(5*dk(1)**4*ddk(1,ie)*dk(3)+dk(1)**5*ddk(3,ie))
          ddf=ddf+210*(2*dk(2)*ddk(2,ie)*dk(4)+dk(2)**2*ddk(4,ie))
          ddf=ddf+70*dk(4)*ddk(4,ie)
          ddf=ddf+70*(4*dk(1)**3*ddk(1,ie)*(3*dk(2)**2+dk(4))+dk(1)**4*(6*dk(2)*ddk(2,ie)+ddk(4,ie)))
          ddf=ddf+56*(ddk(3,ie)*dk(5)+dk(3)*ddk(5,ie))
          ddf=ddf+56*(3*dk(1)**2*ddk(1,ie)*(10*dk(2)*dk(3)+dk(5))+dk(1)**3*(10*(ddk(2,ie)*dk(3)+dk(2)*ddk(3,ie))+ddk(5,ie)))
          ddf=ddf+28*(ddk(2,ie)*(10*dk(3)**2+dk(6))+dk(2)*(20*dk(3)*ddk(3,ie)+ddk(6,ie)))
          ddf=ddf+28*(2*dk(1)*ddk(1,ie)*(15*dk(2)**3+10*dk(3)**2+15*dk(2)*dk(4)+dk(6))+
     &                dk(1)**2*(45*dk(2)**2*ddk(2,ie)+20*dk(3)*ddk(3,ie)+15*(ddk(2,ie)*dk(4)+dk(2)*ddk(4,ie))+ddk(6,ie)))
          ddf=ddf+8*(ddk(1,ie)*(105*dk(2)**2*dk(3)+35*dk(3)*dk(4)+21*dk(2)*dk(5)+dk(7))+
     &               dk(1)*(105*(2*dk(2)*ddk(2,ie)*dk(3)+dk(2)**2*ddk(3,ie))+35*(ddk(3,ie)*dk(4)+dk(3)*ddk(4,ie))+
     &                      21*(ddk(2,ie)*dk(5)+dk(2)*ddk(5,ie))+ddk(7,ie)))
          ddf=ddf+ddk(8,ie)
        endif

      else

        stop 'nlandau bigger than 8 not permitted in cbasis_fns.f!'

      endif

      didjnflux=ddf

      return
      end





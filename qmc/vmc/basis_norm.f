      subroutine basis_norm(iwf,iflag)
c Written by Cyrus Umrigar
c Set normalization of basis fns.
c In 3d:
c Norm of radial part: ((2*zeta)^{2n+1}/(2n)!)^{1/2} for Slaters (n>0).
c                      (2*(2*zeta)^(n+1/2)/Gamma(n+1/2))^{1/2} for gaussians
c where                Gamma(n+1/2) = integral {t^{n+1/2} e^-t dt} = (2n-1)!! sqrt(pi)/2^n for gaussians
c                      Gamma(1/2)=sqrt(pi), Gamma(a+1)=a*Gamma(a), Gamma(a)=(a-1)!
c obtained by Integrate[r^2 (r^{n-1} Exp[-zeta r])^2,{r,0,Infinity}] for Slaters
c obtained by Integrate[r^2 (r^{n-1} Exp[-zeta r^2])^2,{r,0,Infinity}] for Gaussians
c Norm of angular part: ((2*l+1)/(4*pi))^{1/2}.
c In 2d:
c Norm of radial part: ((2*zeta)^{2n}/(2n-1)!)^{1/2}.
c Norm of angular part: (min(m+1,2)/(2*pi))^{1/2}.
c If numr =0 or -1 we are using analytic basis functions and we use normalization
c                  for angular and radial parts.  The -1 is just to tell it to order
c                  basis functions by all the s's first, then all the p's etc.
c                  instead of 1s,2s,2p,3s,3p,3d,...
c         =1       we are using numerical basis functions and we use normalization
c                  for angular part only
c Whether one is using Slater or gaussian basis fns. is inputted by having
c n1s,n2s etc. be either > 0 or < 0.
c The two old versions of the code used unnormalized Gaussian and asymptotic functions,
c and, the same normal. for Gaussians as for Slaters.

      use all_tools_mod
      use control_mod
      use orbitals_mod
      use coefs_mod
      use basis1_mod
      use const_mod
      implicit real*8(a-h,o-z)
      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
c    &,iwctype(MCENT),nctype,ncent
      common /contrl_per/ iperiodic,ibasis
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
      common /basis2/ zex2(MRWF,MCTYPE,MWF),n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),iwrwf2(MBASIS)
c anorm stored for reuse in fit.  Since iwf=1 in fit, we omit iwf dependence.
      common /basisnorm/ anorm(MBASIS)
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap       !JT
     &,ifock,i3body,irewgt,iaver,istrch                                !JT
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt   !JT

      do 20 ib=1,nbasis
        n=n_bas(ib)
        if(ndim.eq.3) then
          l=l_bas(ib)
          if(numr.le.0) then
            if(n.gt.0) then
              anorm(ib)=sqrt((2*zex(ib,iwf))**(2*n+1)*(2*l+1)/(fact(2*n)*4*pi))
             else
              n1=abs(n)
              anorm(ib)=sqrt(2*(2*zex(ib,iwf))**(n1+0.5d0)*(2*l+1)/(gamma1(n1)*4*pi))
            endif
           elseif(numr.eq.1 .or. numr.eq.-1 .or. numr.eq.-2 .or. numr.eq.-3) then
            anorm(ib)=sqrt((2*l+1)/(4*pi))
           else
            stop 'numr must be between -3 and 1'
          endif
         elseif(ndim.eq.2) then
          m=m_bas(ib)
          if(numr.le.0 .and. ibasis.lt.4) then
            anorm(ib)=sqrt((2*zex(ib,iwf))**(2*n)*min(abs(m)+1,2)/(fact(2*n-1)*2*pi))
           elseif(numr.le.0 .and. (ibasis.ge.4 .and. ibasis.le.6)) then
            anorm(ib)=dsqrt(1/pi)
           else
c Warning: temporarily commented out diff norm for m=0
c           anorm(ib)=sqrt(min(abs(m)+1,2)/(2*pi))
            anorm(ib)=sqrt(1/pi)
          endif
        endif
        if(iflag.eq.1) then
          do iorb=1,norb
            coef(ib,iorb,iwf)=coef(ib,iorb,iwf)*anorm(ib)
          enddo
        endif
   20 continue
      if(ipr.ge.0) write(6,'(''anorm='',20f10.6)') (anorm(ib),ib=1,nbasis)

      call object_modified ('anorm')   !JT
!      write(*,*) 'JT icusp=',icusp
!      if (icusp.eq.1) then
!       call ie
!       call equiv_bas
!       call impose_cusp_en_orb_occ
!      endif
      call object_modified ('coef')    !JT

      return
      end
c-----------------------------------------------------------------------

      function fact(n)
      implicit real*8(a-h,o-z)

      fact=1
      do 10 i=2,n
   10   fact=fact*i
      return
      end
c-----------------------------------------------------------------------

      function gamma1(n)
c W. Al-Saidi 4/11/07
c Used for norm of 3D Guassians, gamma1(n) = gamma(n+1/2) = integral {t^{n+1/2} e^-t dt} = (2n-1)!! sqrt(pi)/2^n
      implicit real*8(a-h,o-z)

      spi=1.77245385090552d0 ! sqrt{pi}

      if(n.eq.0) then
        gamma1= spi
       elseif(n.eq.1) then
         gamma1= spi/2
       elseif(n.eq.2) then
         gamma1=3*spi/4
       elseif(n.eq.3) then
         gamma1=15*spi/8
       elseif(n.eq.4) then
         gamma1=105*spi/16
       elseif(n.eq.5) then
         gamma1=945*spi/32
       elseif(n.eq.6) then
         gamma1=10395*spi/64
       else
         write(6,'(''gamma1 is not implemented  in basis_norm.f for n='',i3)') n
         stop "gamma1 is not implemented  in basis_norm.f for this n"
      endif
      return
      end

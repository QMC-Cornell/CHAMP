      subroutine basis_norm_dot(iwf,iflag)

c Written by A.D.Guclu, Feb 2004.to replace basis_norm.f
c for complex wf quantum dot calculations
c calculates normalization of Fock-Darwin basis functions to
c avoid multiple calculations in cbasis*.f
c output -> anorm is calculated
c           coef updated if iflag=1
c iflag is useless for the moment since we don't optimize in fit
c any "exponent term" in coefficients for dots. maybe could be useful later?
c also see basis_norm.dot for real wfs

      use all_tools_mod
      use const_mod
      use coefs_mod
      use basis1_mod
      use basisnorm_mod
      implicit real*8(a-h,o-z)

!     include 'vmc.h'
!     include 'force.h'

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
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
      common /basis3/ n_fd(MBASIS),m_fd(MBASIS),n_cf(MBASIS),ncfmax
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /compferm/ emagv,nv,idot
!JT      common /basisnorm/ anorm(MBASIS)

      do 20 ib=1,nbasis
        n=n_fd(ib)
        mabs=abs(m_fd(ib))
        fac=1.d0
        if(mabs.gt.0) then
          do 5 i=1,mabs
             fac=fac*(n+i)
 5        enddo
        endif
        if(idot.eq.3) then
c we need this to be consistent with G.S.Jeon&J.Jain's convention:
          anorm(ib)=(-1.d0)**n_cf(ib)
        else
          anorm(ib)=dsqrt((zex(ib,iwf)*we)**(mabs+1)/pi/fac)
        endif
c        write(*,*) 'ib,n_fd,m_fd,anorm=',ib,n_fd(ib),m_fd(ib),anorm(ib)
        if(iflag.eq.1) then
          do 10 iorb=1,norb
            coef(ib,iorb,iwf)=coef(ib,iorb,iwf)*anorm(ib)
 10       enddo
        endif
 20   enddo

      call object_modified ('anorm')   !JT
      call object_modified ('coef')    !JT

      return
      end

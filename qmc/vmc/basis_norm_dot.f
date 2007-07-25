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

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE)
     &,n5f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis3/ n_fd(MBASIS),m_fd(MBASIS),n_cf(MBASIS),ncfmax
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /compferm/ emagv,nv,idot
      common /basisnorm/ anorm(MBASIS)

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

      return
      end

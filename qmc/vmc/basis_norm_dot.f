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
      use basis3_mod
      implicit real*8(a-h,o-z)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot

      call alloc ('anorm', anorm, nbasis)

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

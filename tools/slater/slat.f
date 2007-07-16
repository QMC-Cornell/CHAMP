      program slater
c Written by Cyrus Umrigar
c Writes out radial Slater basis functions on a grid

      implicit real*8(a-h,o-z)
      include 'numbas.h'

      common /dim/ ndim
      common nrwf,n_rbas(MRWF),l_rbas(MRWF),zex(MRWF),rwf(MRWF),anorm(MRWF)

c nrbas = number of radial basis functions for each center

c igrid = 1 linear r(i+1)=exp_h_bas+r(i), r(1)=r0_bas
c         2 exponential r(i+1)=exp_h_bas*r(i), r(1)=r0_bas
c         3 shifted exponential r(i+1)=r0_bas*(exp_h_bas**(i-1)-1)
c           r(n) is read in, r0_bas=r(n)/(exp_h_bas**(nr-1)-1)

      ndim=3

      read(5,*) nr,rn,exp_h_bas
      r0_bas=rn/(exp_h_bas**(nr-1)-1)
      read(5,*) nrwf
      if(nrwf.gt.MRWF) stop 'nrwf > MRWF'
      do 10 i=1,nrwf
   10 read(5,*) n_rbas(i),l_rbas(i),zex(i)

      call basis_norm(1,1)

      do 40 ir=1,nr
        r=r0_bas*(exp_h_bas**(ir-1)-1)+1.d-20
        call basis_fns(r)
   40   write(6,'(1p,99d20.12)') r,(anorm(i)*rwf(i),i=1,nrwf)

      stop
      end
c-----------------------------------------------------------------------
      subroutine basis_fns(r)
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
      include 'numbas.h'

      common nrwf,n_rbas(MRWF),l_rbas(MRWF),zex(MRWF),rwf(MRWF),anorm(MRWF)

      pu=4*datan(1.d0)

            do 10 irb=1,nrwf
              n=n_rbas(irb)
              rn=abs(n)
c             rnm2=max(0.d0,rn-2)
              rm3=r**(rn-3)
              rm2=rm3*r
              rm1=rm2*r
c Slater r^(n-1)*Exp(-zeta*r)
              if(n.gt.0) then
                zr=zex(irb)*r
                ex=dexp(-zr)
               elseif(n.eq.0) then
c Warning: Asymptotic and Gaussian not yet tested.
c Asymptotic r^(rn-1)*Exp(-zeta*r), where rn=beta+1, beta=betaq/zeta-1, zeta=sqrt(-2*E_ion)?
                stop 'asymptotic not yet fully tested'
                rn=betaq/zex(irb)
                rm3=r**(rn-3)
                rm2=rm3*r
                rm1=rm2*r
               elseif(n.lt.0) then
c Gaussian  r^(n-1)*Exp(-zeta*r^2) (to be checked)
c               rnm2=rn-2
                zr=2*zex(irb)*r2
                ex=dexp(-0.5d0*zr)
              endif
   10         rwf(irb)=rm1*ex

      return
      end
c-----------------------------------------------------------------------

      subroutine basis_norm(iwf,iflag)
c Written by Cyrus Umrigar
c Set normalization of basis fns.
c In 3d:
c Norm of radial part: ((2\zeta)^{2n+1}/(2n)!)^{1/2}. 
c Norm of angular part: ((2*l+1)/(4*pi))^{1/2}. 
c In 2d:
c Norm of radial part: ((2\zeta)^{2n}/(2n-1)!)^{1/2}. 
c Norm of angular part: (min(m+1,2)/(2*pi))^{1/2}. 
c At present for Gaussians and asymptotic we are using same norm as for Slaters,
c but the correct norm for Gaussians is sqrt(2^(2n+3)*zeta^(2n+1))/Gamma(n+1/2).
c Gamma(1/2)=sqrt(pi) and Gamma(a+1)=a*Gamma(a).
c The old version of the code used unnormalized Gaussian and asymptotic functions.

      implicit real*8(a-h,o-z)
      include 'numbas.h'

      common /dim/ ndim
      common nrwf,n_rbas(MRWF),l_rbas(MRWF),zex(MRWF),rwf(MRWF),anorm(MRWF)

      pi=4*datan(1.d0)

      do 20 ib=1,nrwf
        n=n_rbas(ib)
        if(ndim.eq.3) then
          l=l_rbas(ib)
          if(numr.eq.0) then
c           anorm(ib)=sqrt((2*zex(ib))**(2*n+1)*(2*l+1)/(fact(2*n)*4*pi))
            anorm(ib)=sqrt((2*zex(ib))**(2*n+1)/(fact(2*n)))
           else
            anorm(ib)=sqrt((2*l+1)/(4*pi))
          endif
c        elseif(ndim.eq.2) then
c         m=m_bas(ib)
c         if(numr.eq.0) then
c           anorm(ib)=sqrt((2*zex(ib))**(2*n)*min(abs(m)+1,2)/(fact(2*n-1)*2*pi))
c          else
ccWarning: temporarily commented out diff norm for m=0
cc          anorm(ib)=sqrt(min(abs(m)+1,2)/(2*pi))
c           anorm(ib)=sqrt(1/pi)
c         endif
        endif
   20 continue

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

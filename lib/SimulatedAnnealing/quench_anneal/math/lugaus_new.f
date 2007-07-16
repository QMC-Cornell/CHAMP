cDate of last modification: Fri May  6 16:15:57 EDT 1994
      subroutine lugaus_new(A,B,n,np,dum1,idum1,idum2,eps,dseed,newA,r,
     &  problog,pronly,T,dum2)
c lugaus: purpose: sample from multivariate gaussian/calculate probablility
c         density prop. to exp[-(r.A.r/2+B.r)/T], where A is asumed positive
c 	  definite.
c in:
c A      = positive definite symmetric matrix of dimension n
c          and physical dimension np
c B      = vector of order n
c pronly = true: calculate probablity for given r
c          false: sample r
c eps    = for compatibility with mvgaus (not used)
c dseed  = seed from random number generator
c newA   = true for the first call with matrix A
c          false for repeat call with same A, but in this case on input
c          A has to be the LU decomposition previously obtained
c r      = n-vector or which probability will be calculated, if
c          pronly = .true.
c out:
c A      = LU decomposition of original A=L.L^{Transpose}
c B      = L^{-1}B
c r      = random n-vector sampled from gaussian
c problo = log( sqrt{det[A/T]/(2 pi)^n} exp(-B.A^(-1).B/2T) exp[-(r.A.r/2+B.r)/T)
c In this normalization factor:
c (a) the factors of 1/sqrt(2 pi) are the usual normalization
c     factors of the N(0,1) distribution.
c (b) Sqrt(det[A]) is the Jacobian due to the transformation to
c     principal axes.
c (c) The exponential factor comes about by completing the square and
c     shifting the coordinates.
      implicit real*8(a-h,o-z)
      parameter (zero=0,one=1,two=2,half=one/two)
      dimension A(np,np),B(n),r(n)
      logical newA,pronly,zeroT
      save detlog
      if(T.eq.zero) then
c T0 has no meaning here it's defined to avoid zero devides
        zeroT=.true.
        T0=one
        T0log=zero
      else
        zeroT=.false.
        T0=T
        T0log=log(T0)
        T0sqrt=sqrt(T0)
      endif
      if(pronly) then
         s2=zero
         do i=1,n
            ri=r(i)
            s2=s2-B(i)*ri
            do j=1,i
               w=one
               if(i.eq.j) w=half
               s2=s2-w*ri*A(i,j)*r(j)
            enddo
         enddo
         call chlsky(A,n,np,detlog)
         call lxb(A,n,np,B)
         s=0
         do i=1,n
            s=s-B(i)**2
         enddo
         problog=detlog+(s2+half*s)/T0-half*n*T0log
         return
      else
         if(newA) then
            call chlsky(A,n,np,detlog)
            do i=1,n
               do j=i+1,n
                  A(i,j)=A(j,i)
               enddo
            enddo
            call lxb(A,n,np,B)
         endif
         s2=zero
         do i=1,n
            if(zeroT) then
              s=zero
            else
              s=T0sqrt*gaushn()
            endif
            s2=s2-s**2
            r(i)=s-B(i)
         enddo
         call uxb(A,n,np,r)
c The variables s are the shifted gaussian coordinates. The Jacobian
c part of the normalization constant is needed to get the probability
c density of the r coordinates.
         problog=detlog+half*(s2/T0-n*T0log)
         return
      endif
      end

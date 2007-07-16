      subroutine chlsky(a,n,np,detlog)
c chlsky: purpose: cholesky decomposition of a and determinant
c in: matrix a of order n stored with physical dimension np
c out: lower triangular matrix stored in lower portion of a
c      detlog= log of sqrt[det{a/(2 pi)}]
c note: lower triangular portion of original a is destroyed
      implicit real*8(a-h,o-z)
      dimension a(np,np)
      parameter (ZERO=0,ONE=1,TWO=2,SXTN=16,SXTNTH=1/SXTN,
     &  REALMIN=1d-300)
      parameter (SQO2PI=0.398942280401432678d0)
      detexp=ZERO
      d1=ONE
      do j=1,n
         if(j.gt.1) then
            jm1=j-1
            do k=j,n
               sum=0
               do ip=1,jm1
                  sum=sum+a(k,ip)*a(j,ip)
               enddo
               a(k,j)=a(k,j)-sum
            enddo
         endif

         if(a(j,j).le.ZERO) a(j,j)=REALMIN

         s=ONE/sqrt(a(j,j))
         do k=j,n
            a(k,j)=a(k,j)*s
         enddo
         d1=d1*a(j,j)*SQO2PI
         if(d1.gt.ONE) then
            d1=d1*SXTNTH
            detexp=detexp+4
         else
            if(d1.lt.SXTNTH) then
               d1=d1*SXTN
               detexp=detexp-4
            endif
         endif
      enddo
      detlog=detexp*log(TWO)+log(d1)
      return
      end

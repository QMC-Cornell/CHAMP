      function gammai(a,x,xae,iflag)
! Written by Cyrus Umrigar
! Calculate the incomplete gamma function

! To speed up the evaluation several changes have been made.
! 1) The series is used for x < a+4 rather that x < a+1
! 2) Some of the abs() have been removed. OK for a>0
! 3) xae is passed in rather than being re-evaluated here.
!    xae = x**a * exp(-x)
! 4) Warning: To speed up the do loops, the test has been removed
! and the number of times the loop is performed is fixed by formulae
! that are appropriate for a=3/2 and 5/2, but not for general a.
! The formulae are:
! a) n=21+3*(x-a) for series, b) max(20+a-x,7) for cont frac
!--------------------------
!   x     5/2        3/2
!      need form  need form
!--------------------------
!  .01   6   14     6   17
!  .5   12   15    13   18
! 1.5   17   18    18   21
! 3.5   24   24    25   27
! 4.5   27   27    28   30
! 5.5   30   30    31   33
! 6.5   32   33
!-------------------------
! 5.5              16   16
! 7.5   13   15    13   14
! 9.5   11   13    12   12
!12.5    9   10    10    9
!25.5    7    7     7    7

! We have solved the problem of round-off as follows:
! The problem occurs for very large x. If both branches of this routine
! were changed to return the complementary fn gammac, then the problem
! would occur for very small x.  The problem is that for
! x large gammai is calculated by calculating its complement gammac and
! subtracting it from gamma.  However for x very large, gammai is almost
! zero and for x > 43 all significant digits are lost when gammac is
! subtracted from gamma. Since the expression for the area in metrop
! involves the difference of gammai's, the gamma's cancel giving area=0.
! The solution is to return gammai if x < a+4 and gammac otherwise.
! This requires changing all the metrop routines that call gammai.
! Returns   gammai, iflag=1,    x <  a+4
!          -gammac, iflag=-1,   x >= a+4 even though it is called gammai

      implicit real*8(a-h,o-z)
!     parameter (itmax=100,eps=1.d-14, d3b2=1.5d0,d5b2=2.5d0,
!    & g3b2=.886226925452758d0,g5b2=1.329340388179137d0)

      if(x.lt.a+4)then
!       write(6,*)'series'
        iflag=1
        ap=a
        sum=1/a
        del=sum
        iter=int(21+3*(x-a))
!i      do 10 n=1,itmax
        do 10 n=1,iter
          ap=ap+1
          del=del*x/ap
          sum=sum+del
!i        if(del.lt.sum*eps)go to 20
   10   continue
!i      stop 'a too large, itmax too small'
!i 20   gammai=sum*dexp(-x+a*dlog(x))
        gammai=sum*xae
       else
!       write(6,*) 'contin frac'
        iflag=-1
!       gold=0
        a0=1
        a1=x
        b0=0
        b1=1
        fac=1
        iter=max(20+nint(a-x),7)
!i      do 30 n=1,itmax
        do 30 n=1,iter
          an=dfloat(n)
          ana=an-a
          a0=(a1+a0*ana)*fac
          b0=(b1+b0*ana)*fac
          anf=an*fac
          a1=x*a0+anf*a1
          b1=x*b0+anf*b1
            fac=1/a1
            g=b1*fac
!           if(abs((g-gold)/g).lt.eps)go to 40
!i          if(abs((g-gold)).lt.eps*g)go to 40
!           gold=g
   30   continue
!i      stop 'a too large, itmax too small'
!i 40   gammcf=g*dexp(-x+a*dlog(x))
        gammcf=g*xae
        gammai=-gammcf
!c      if(a.eq.d3b2) then
!c        gammai=g3b2-gammcf
!c       elseif(a.eq.d5b2) then
!c        gammai=g5b2-gammcf
!c       else
!c        gammai=gamm(a)-gammcf
!c      endif
      endif
!     write(6,*) xae,dexp(-x+a*dlog(x))
      return
      end
!-----------------------------------------------------------------------

      function gamm(a)
! Written by Cyrus Umrigar
! Evaluate the Gamma function.
! For a=3/2 the error of this formula is 5.d-11
      implicit real*8(a-h,o-z)
      dimension cof(6)
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0, &
     & -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/.5d0,1.d0,5.5d0/
      x=a-one
      tmp=x+fpf
!     tmp=(x+half)*dlog(tmp)-tmp
      tmp=tmp**(x+half)*dexp(-tmp)
      ser=one
      do 10 j=1,6
        x=x+one
   10   ser=ser+cof(j)/x
!     gammln=tmp+dlog(stp*ser)
!     gamm=dexp(gammln)
      gamm=stp*ser*tmp
      return
      end


!-----------------------------------------------------------------------

      function gammai_u0(x)
!     Written by Abhijit Mehta
!     Calculates the upper incomplete gamma function for a = 0
!       \int_x^\inf  (t^(a-1) exp(-t)) dt
!     Based on routine for Exponential integral E_1(x) in NR
!     Uses a partial fraction expansion for x > 1
!     and uses a series for 0 < x < 1.

      implicit real*8(a-h,o-z)
      parameter (itmax=100,eps=1.d-14,bignum=1.d40,euler=0.577215664902)
      if(x.le.0)then
         write(6,*) 'x must be greater than 0 in gammai(0,x)'
         stop 'x must be greater than 0 in gammai(0,x)'
      else if(x.gt.1)then
         b = x+1
         c = bignum
         d = 1./b
         h = d
         do i=1,itmax
            a = -i*i
            b = b+2
            d = 1./(a*d+b)
            c = b+a/c
            del = c*d
            h = h*del
            if(abs(del-1.).lt.eps)then
               gammai_u0 = h*dexp(-x)
               return
            endif
         enddo
         write(6,*) 'gammai_u0(x) failed in continued fraction (x>1)'
         stop 'gammai_u0(x) failed in continued fraction (x>1)'
      else
         gammai_u0 = -dlog(x) - euler
         f = 1
         do i=1,itmax
            f = -f*x/i
            gammai_u0 = gammai_u0 - f/i
            if(abs(f/i).lt.abs(gammai_u0)*eps) return
         enddo
         write(6,*) 'gammai_u0(x) failed in series (x<=1)'
         stop 'gammai_u0(x) failed in series (x<=1)'
      endif
      return
      end



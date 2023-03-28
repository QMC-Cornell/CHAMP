      FUNCTION ei(x)
      INTEGER MAXIT
      REAL*8 ei,x,EPS,EULER,FPMIN
! EULER is the Euler-Mascheroni constant, EulerGamma in Mathematica
      PARAMETER (EPS=6.d-8,EULER=.5772156649015328d0,MAXIT=100)
      PARAMETER (FPMIN=1.d-30)
      INTEGER k
      REAL fact,prev,sum,term
      if(x.le.0.) pause 'bad argument in ei'
      if(x.lt.FPMIN)then
        ei=log(x)+EULER
      else if(x.le.-log(EPS))then
        sum=0.
        fact=1.
        do 11 k=1,MAXIT
          fact=fact*x/k
          term=fact/k
          sum=sum+term
          if(term.lt.EPS*sum)goto 1
11      continue
        pause 'series failed in ei'
1       ei=sum+log(x)+EULER
      else
        sum=0.
        term=1.
        do 12 k=1,MAXIT
          prev=term
          term=term*k/x
          if(term.lt.EPS)goto 2
          if(term.lt.prev)then
            sum=sum+term
          else
            sum=sum-prev
            goto 2
          endif
12      continue
2       ei=exp(x)*(1.+sum)/x
      endif
      return
      END

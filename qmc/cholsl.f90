      SUBROUTINE cholsl(a,n,np,p,b,x)
      INTEGER n,np
      REAL*8 a(np,np),b(n),p(n),x(n)
      INTEGER i,k
      REAL sum
      do 12 i=1,n
        sum=b(i)
        do 11 k=i-1,1,-1
          sum=sum-a(i,k)*x(k)
11      continue
        x(i)=sum/p(i)
12    continue
      do 14 i=n,1,-1
        sum=x(i)
        do 13 k=i+1,n
          sum=sum-a(k,i)*x(k)
13      continue
        x(i)=sum/p(i)
14    continue
      return
      END

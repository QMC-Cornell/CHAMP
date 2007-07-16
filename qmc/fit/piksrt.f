      subroutine piksrt(n,arr)
c Written by Claudia Filippi

      implicit real*8(a-h,o-z)

      integer arr(n),a

      do 12 j=2,n
        a=arr(j)
        do 11 i=j-1,1,-1
          if(arr(i).le.a) goto 10
          arr(i+1)=arr(i)
11      continue
        i=0
10      arr(i+1)=a
12    continue

      return
      end

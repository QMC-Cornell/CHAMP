      subroutine svbksb(u,w,v,m,n,mp,np,b,x,work)

      implicit real*8(a-h,o-z)

      parameter (zero=0.d0)
      dimension b(*),u(mp,*),v(np,*),w(*),x(*),work(*)

      do 12 j=1,n
        s=zero
        if(dabs(w(j)).gt.zero) then
          s=ddot(m,b,1,u(1,j),1)
          s=s/w(j)
        endif
 12     work(j)=s

      do 14 j=1,n
        s=zero
        do 13 jj=1,n
 13       s=s+v(j,jj)*work(jj)
 14     x(j)=s

      return
      end

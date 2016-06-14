      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter,itmax,converged)
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      DOUBLE PRECISION ftol,p(mp,np),y(mp),funk
!JT      PARAMETER (NMAX=200,ITMAX=5000)
      PARAMETER (NMAX=200)
      EXTERNAL funk
!U    USES amotry,funk
      INTEGER i,ihi,ilo,inhi,j,m,n
      DOUBLE PRECISION rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      LOGICAL converged
      converged = .false. !JT
      iter=0
1     do 12 n=1,ndim
        sum=0.d0
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        converged = .true. !JT
        return
      endif
!JT      if (iter.ge.ITMAX) pause 'ITMAX exceeded in amoeba'
      if (iter.ge.ITMAX) return !JT
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0d0)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0d0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum)
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END

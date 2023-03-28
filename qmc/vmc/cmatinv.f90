      subroutine cmatinv(a,n,det)

! from Wolfgang's qdot.f
! invert a complex nxn matrix using gauss elimination
! with row pivoting. Note matrix must be dimension (n,n) or equivalently
      use const_mod
      implicit double precision (a-h,o-z)

      parameter (nmax=100)
      complex(dpc) det,adiag,adiagi,t,cone,czero,atemp
      complex(dpc) a(nelec,nelec)

      dimension atemp(nmax),ipvt(nmax)
      cone=dcmplx(1.d0,0.d0)
      czero=dcmplx(0.d0,0.d0)
      if(n.gt.nmax) then
        write (6,'(1x,'' nmax too small in cmati'')')
        stop
      endif
      do 10 i=1,n
   10   ipvt(i)=i
      det=cone

! loop through columns
      do 20 i=1,n
        adiag=a(ipvt(i),i)
        idiag=i

! find best pivot element in column and record pivot

      do 30 k=i,n
      if(abs(a(ipvt(k),i)).gt.abs(adiag)) then
         adiag=a(ipvt(k),i)
         idiag=k
         endif
   30 continue
      if(idiag.ne.i) then
         det=-det
         itemp=ipvt(i)
         ipvt(i)=ipvt(idiag)
         ipvt(idiag)=itemp
      endif
      det=adiag*det

! row reduce matrix

      a(ipvt(i),i)=cone
      adiagi=cone/adiag
      do 40 k=1,n
   40 a(ipvt(i),k)=a(ipvt(i),k)*adiagi
      do 50 j=1,n
      if(j.ne.ipvt(i)) then
         t=-a(j,i)
         a(j,i)=czero
         do 60 k=1,n
   60    a(j,k)=a(j,k)+t*a(ipvt(i),k)
         endif
   50 continue
   20 continue

! interchange elements to unpivot inverse matrix
! the following is equivalent to:
!      anew(i,ipvt(j))=aold(ipvt(i),j)

      do 70 j=1,n
      do 80 i=1,n
   80 atemp(i)=a(i,j)
      do 90 i=1,n
   90 a(i,j)=atemp(ipvt(i))
   70 continue
      do 100 i=1,n
      do 110 j=1,n
  110 atemp(j)=a(i,j)
      do 120 j=1,n
  120 a(i,ipvt(j))=atemp(j)
  100 continue
      return
      end

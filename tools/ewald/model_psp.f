c write out a model pseudopotential that is cubic or quartic at the
c origin and goes as -1/r at infinity

      implicit real*8(a-h,o-z)
  
      read(2,*) Z,a
c     read(2,*) n,rmax,arg_ps
      read(2,*) n,arg_ps
      do 10 i=1,n
c       r=max(r0_ps*(arg_ps**(i-1)-1),1.d-9)
        read(2,*) r,y
   10   write(6,'(1p9d20.12)') r, y, y+z/r

      stop
      end

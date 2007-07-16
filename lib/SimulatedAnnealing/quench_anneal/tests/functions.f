      real*8 function func(ndata, nparm, parm, diff, ioffset)

      implicit none

      include 'dimen.h'

      integer ndata, nparm, ioffset
      real*8 parm(nparm), diff(ndata)

      integer icalls
      common / fcn_calls / icalls

      real*8 x, y
      common / coord / x(MDATA), y(MDATA)

      integer i
      real*8 valley

      func=0.0d0

      if ( ioffset .ge. 0 ) then

         icalls = icalls + 1
         do i=1,ndata
            diff(i) = valley( parm, nparm, x(i+ioffset) ) - y(i+ioffset)
            func = func + diff(i)**2
         enddo

      else

         do i=1,ndata
            func = func + diff(i)**2
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine init_func(ptrue, nparm, ndata)

      implicit none

      include 'dimen.h'

      integer nparm, ndata
      real*8 ptrue(nparm)

      real*8 x, y
      common / coord / x(MDATA), y(MDATA)

      real*8 valley

      integer i, iseed, idata
      real*8 widthf, chisqtrue, f, r, err


c     Choose x values.
      do i = 1, ndata
         x(i) = dble(i)
      enddo

c     Initialize random number generator.
      iseed = 12321
      call ransi(iseed)

c     Width of noise distribution added to data points.
      widthf = 0.0d0

      chisqtrue = 0.0d0
      do idata = 1, ndata

         call ransr( r, 1 )
         f = valley( ptrue, nparm, x(idata) )
         err = f * widthf * ( r - 0.5d0 )
         y(idata) = f + err
         chisqtrue = chisqtrue + err**2

      enddo

c     write(6,'((''y   '',6(1x,g12.6)))') (y(i),i=1,ndata)
c     write(6,'((''parm'',6(1x,g12.6)))') parm

      return
      end
c-----------------------------------------------------------------------

      real*8 function valley(parm, nparm, x)

      implicit none

      include 'dimen.h'

      real*8 PI
      parameter( PI = 3.141592653589793d0 )

      integer nparm
      real*8 parm(nparm), x

c     valley = 10.0d0**4 * ( parm(2) - parm(1) - parm(1)**2 ) * x
      valley = 10.0d0**4 * ( parm(2) - sin( PI * parm(1) ) ) * x
c    & + 10.0d0**2 * ( parm(2) - parm(1) ) * x**2
c    & + 10.0d0**2 * parm(2) * x**2
     & + ( parm(1)**2 + parm(2)**2 ) * x**3
c     write(6,*) x, valley

      return
      end

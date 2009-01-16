      subroutine r8evbicub(xget,yget,x,nx,y,ny,ilinx,iliny,
     >                   f,inf2,ict,fval,ier)

C  evaluate a 2d cubic Spline interpolant on a rectilinear
C  grid -- this is C2 in both directions.

C  this subroutine calls two subroutines:
C     herm2xy  -- find cell containing (xget,yget)
C     fvbicub  -- evaluate interpolant function and (optionally) derivatives

C  input arguments:
C  ================

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inf2
!============
      integer nx,ny                     ! grid sizes
      REAL*8 xget,yget                    ! target of this interpolation
      REAL*8 x(nx)                        ! ordered x grid
      REAL*8 y(ny)                        ! ordered y grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
      integer iliny                     ! iliny=1 => assume y evenly spaced

      REAL*8 f(0:3,inf2,ny)               ! function data

C       f 2nd dimension inf2 must be .ge. nx
C       contents of f:

C  f(0,i,j) = f @ x(i),y(j)
C  f(1,i,j) = d2f/dx2 @ x(i),y(j)
C  f(2,i,j) = d2f/dy2 @ x(i),y(j)
C  f(3,i,j) = d4f/dx2dy2 @ x(i),y(j)

C      (these are spline coefficients selected for continuous 2-
C      diffentiability, see mkbicub[w].for)

      integer ict(6)                    ! code specifying output desired

C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return df/dy  (0, don't)
C  ict(4)=1 -- return d2f/dx2  (0, don't)
C  ict(5)=1 -- return d2f/dy2  (0, don't)
C  ict(6)=1 -- return d2f/dxdy (0, don't)
c                   the number of non zero values ict(1:6)
c                   determines the number of outputs...

c  new dmc December 2005 -- access to higher derivatives (even if not
c  continuous-- but can only go up to 3rd derivatives on any one coordinate.
c     if ict(1)=3 -- want 3rd derivatives
c          ict(2)=1 for d3f/dx3
c          ict(3)=1 for d3f/dx2dy
c          ict(4)=1 for d3f/dxdy2
c          ict(5)=1 for d3f/dy3
c               number of non-zero values ict(2:5) gives no. of outputs
c     if ict(1)=4 -- want 4th derivatives
c          ict(2)=1 for d4f/dx3dy
c          ict(3)=1 for d4f/dx2dy2
c          ict(4)=1 for d4f/dxdy3
c               number of non-zero values ict(2:4) gives no. of outputs
c     if ict(1)=5 -- want 5th derivatives
c          ict(2)=1 for d5f/dx3dy2
c          ict(3)=1 for d5f/dx2dy3
c               number of non-zero values ict(2:3) gives no. of outputs
c     if ict(1)=6 -- want 6th derivatives
c          d6f/dx3dy3 -- one value is returned.

C output arguments:
C =================

      REAL*8 fval(*)                      ! output data
      integer ier                       ! error code =0 ==> no error

C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C  fval(3) receives the third output (depends on ict(...) spec)
C  fval(4) receives the fourth output (depends on ict(...) spec)
C  fval(5) receives the fourth output (depends on ict(...) spec)
C  fval(6) receives the fourth output (depends on ict(...) spec)

C  examples:
C    on input ict = [1,1,1,0,0,1]
C   on output fval= [f,df/dx,df/dy,d2f/dxdy], elements 5 & 6 not referenced.

C    on input ict = [1,0,0,0,0,0]
C   on output fval= [f] ... elements 2 -- 6 never referenced.

C    on input ict = [0,0,0,1,1,0]
C   on output fval= [d2f/dx2,d2f/dy2] ... elements 3 -- 6 never referenced.

C    on input ict = [0,0,1,0,0,0]
C   on output fval= [df/dy] ... elements 2 -- 6 never referenced.

C  ier -- completion code:  0 means OK
C-------------------
C  local:

      integer i,j                       ! cell indices

C  normalized displacement from (x(i),y(j)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C    yparam=0 @y(j)  yparam=1 @y(j+1)

      REAL*8 xparam,yparam

C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))

      REAL*8 hx,hy
      REAL*8 hxi,hyi

C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1

C  ** the interface is very similar to herm2ev.for; can use herm2xy **
C---------------------------------------------------------------------

      call r8herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i,j,xparam,yparam,hx,hxi,hy,hyi,ier)
      if(ier.ne.0) return

      call r8fvbicub(ict,1,1,
     >   fval,i,j,xparam,yparam,hx,hxi,hy,hyi,f,inf2,ny)

      return
      end
C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 2d fcn
C   --vectorized-- dmc 10 Feb 1999

      subroutine r8fvbicub(ict,ivec,ivecd,
     >   fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >   fin,inf2,ny)

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ny,inf2,iadr,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 z36th,xp,xpi,xp2,xpi2,cx,cxi,hx2,yp,ypi,yp2,ypi2,cy
      REAL*8 cyi,hy2,cxd,cxdi,cyd,cydi
!============
      integer ict(6)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)

      integer ii(ivec),jj(ivec)         ! target cells (i,j)
      REAL*8 xparam(ivec),yparam(ivec)
                          ! normalized displacements from (i,j) corners

      REAL*8 hx(ivec),hy(ivec)            ! grid spacing, and
      REAL*8 hxi(ivec),hyi(ivec)          ! inverse grid spacing 1/(x(i+1)-x(i))
                                        ! & 1/(y(j+1)-y(j))

      REAL*8 fin(0:3,inf2,ny)             ! interpolant data (cf "evbicub")

      REAL*8 fval(ivecd,*)                ! output returned

C  for detailed description of fin, ict and fval see subroutine
C  evbicub comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.

C  note that the index inputs ii,jj and parameter inputs
C     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument

C  to use this routine in scalar mode, pass in ivec=ivecd=1

C---------------
C  Spline evaluation consists of a "mixing" of the interpolant
C  data using the linear functionals xparam, xpi = 1-xparam,
C  yparam, ypi = 1-yparam, and the cubic functionals
C  xparam**3-xparam, xpi**3-xpi, yparam**3-yparam, ypi**3-ypi ...
C  and their derivatives as needed.

      integer v
      REAL*8 sum

      REAL*8, parameter :: sixth = 0.166666666666666667_r8

C---------------
C   ...in x direction

      z36th=sixth*sixth
      iadr=0

      if(ict(1).le.2) then

C  get desired values:

         if(ict(1).eq.1) then

C  function value:

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cx=xp*(xp2-1.0_r8)
               cxi=xpi*(xpi2-1.0_r8)
               hx2=hx(v)*hx(v)

C   ...and in y direction

               yp=yparam(v)
               ypi=1.0_r8-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cy=yp*(yp2-1.0_r8)
               cyi=ypi*(ypi2-1.0_r8)
               hy2=hy(v)*hy(v)

               sum=xpi*(ypi*fin(0,i,j)  +yp*fin(0,i,j+1))+
     >              xp*(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1))

               sum=sum+sixth*hx2*(
     >              cxi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+
     >              cx*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))

               sum=sum+sixth*hy2*(
     >              xpi*(cyi*fin(2,i,j)  +cy*fin(2,i,j+1))+
     >              xp*(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))

               sum=sum+z36th*hx2*hy2*(
     >              cxi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+
     >              cx*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(2).eq.1) then

C  df/dx:

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0_r8*xp2-1.0_r8
               cxdi=-3.0_r8*xpi2+1.0_r8

C   ...and in y direction

               yp=yparam(v)
               ypi=1.0_r8-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cy=yp*(yp2-1.0_r8)
               cyi=ypi*(ypi2-1.0_r8)
               hy2=hy(v)*hy(v)

               sum=hxi(v)*(
     >              -(ypi*fin(0,i,j)  +yp*fin(0,i,j+1))
     >              +(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1)))

               sum=sum+sixth*hx(v)*(
     >              cxdi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+
     >              cxd*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))

               sum=sum+sixth*hxi(v)*hy2*(
     >              -(cyi*fin(2,i,j)  +cy*fin(2,i,j+1))
     >              +(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))

               sum=sum+z36th*hx(v)*hy2*(
     >              cxdi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+
     >              cxd*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(3).eq.1) then

C  df/dy:

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cx=xp*(xp2-1.0_r8)
               cxi=xpi*(xpi2-1.0_r8)
               hx2=hx(v)*hx(v)

C   ...and in y direction

               yp=yparam(v)
               ypi=1.0_r8-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0_r8*yp2-1.0_r8
               cydi=-3.0_r8*ypi2+1.0_r8

               sum=hyi(v)*(
     >              xpi*(-fin(0,i,j)  +fin(0,i,j+1))+
     >              xp*(-fin(0,i+1,j)+fin(0,i+1,j+1)))

               sum=sum+sixth*hx2*hyi(v)*(
     >              cxi*(-fin(1,i,j)  +fin(1,i,j+1))+
     >              cx*(-fin(1,i+1,j)+fin(1,i+1,j+1)))

               sum=sum+sixth*hy(v)*(
     >              xpi*(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1))+
     >              xp*(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))

               sum=sum+z36th*hx2*hy(v)*(
     >              cxi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+
     >              cx*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(4).eq.1) then

C  d2f/dx2:

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp

C   ...and in y direction

               yp=yparam(v)
               ypi=1.0_r8-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cy=yp*(yp2-1.0_r8)
               cyi=ypi*(ypi2-1.0_r8)
               hy2=hy(v)*hy(v)

               sum=(
     >              xpi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+
     >              xp*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))

               sum=sum+sixth*hy2*(
     >              xpi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+
     >              xp*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(5).eq.1) then

C  d2f/dy2:

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cx=xp*(xp2-1.0_r8)
               cxi=xpi*(xpi2-1.0_r8)
               hx2=hx(v)*hx(v)

C   ...and in y direction

               yp=yparam(v)
               ypi=1.0_r8-yp

               sum=(
     >              xpi*(ypi*fin(2,i,j)  +yp*fin(2,i,j+1))+
     >              xp*(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))

               sum=sum+sixth*hx2*(
     >              cxi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+
     >              cx*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(6).eq.1) then

C  d2f/dxdy:

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0_r8*xp2-1.0_r8
               cxdi=-3.0_r8*xpi2+1.0_r8

C   ...and in y direction

               yp=yparam(v)
               ypi=1.0_r8-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0_r8*yp2-1.0_r8
               cydi=-3.0_r8*ypi2+1.0_r8

               sum=hxi(v)*hyi(v)*(
     >              fin(0,i,j)  -fin(0,i,j+1)
     >              -fin(0,i+1,j)+fin(0,i+1,j+1))

               sum=sum+sixth*hx(v)*hyi(v)*(
     >              cxdi*(-fin(1,i,j)  +fin(1,i,j+1))+
     >              cxd*(-fin(1,i+1,j)+fin(1,i+1,j+1)))

               sum=sum+sixth*hxi(v)*hy(v)*(
     >              -(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1))
     >              +(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))

               sum=sum+z36th*hx(v)*hy(v)*(
     >              cxdi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+
     >              cxd*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

C-------------------------------------------------

      else if(ict(1).eq.3) then
         if(ict(2).eq.1) then
c  evaluate d3f/dx3 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               yp=yparam(v)
               ypi=1.0_r8-yp
               yp2=yp*yp
               ypi2=ypi*ypi
               cy=yp*(yp2-1.0_r8)
               cyi=ypi*(ypi2-1.0_r8)
               hy2=hy(v)*hy(v)
               sum=hxi(v)*(
     >              -(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))
     >              +(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))

               sum=sum+sixth*hy2*hxi(v)*(
     >              -(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))
     >              +(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(3).eq.1) then
c  evaluate d3f/dx2dy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               xp=xparam(v)
               xpi=1.0_r8-xp
               yp=yparam(v)
               ypi=1.0_r8-yp
               yp2=yp*yp
               ypi2=ypi*ypi
               cyd=3.0_r8*yp2-1.0_r8
               cydi=-3.0_r8*ypi2+1.0_r8

               sum=hyi(v)*(
     >              xpi*(-fin(1,i,j)  +fin(1,i,j+1))+
     >              xp*(-fin(1,i+1,j) +fin(1,i+1,j+1)))

               sum=sum+sixth*hy(v)*(
     >              xpi*(cydi*fin(3,i,j) +cyd*fin(3,i,j+1))+
     >              xp*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(4).eq.1) then
c  evaluate d3f/dxdy2
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi
               cxd=3.0_r8*xp2-1.0_r8
               cxdi=-3.0_r8*xpi2+1.0_r8
               yp=yparam(v)
               ypi=1.0_r8-yp

               sum=hxi(v)*(
     >              -(ypi*fin(2,i,j)  +yp*fin(2,i,j+1))
     >              +(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))

               sum=sum+sixth*hx(v)*(
     >              cxdi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+
     >              cxd*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(5).eq.1) then
c  evaluate d3f/dy3 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cx=xp*(xp2-1.0_r8)
               cxi=xpi*(xpi2-1.0_r8)
               hx2=hx(v)*hx(v)

               sum=hyi(v)*(
     >              xpi*(-fin(2,i,j)  +fin(2,i,j+1))+
     >              xp*(-fin(2,i+1,j) +fin(2,i+1,j+1)))

               sum=sum+sixth*hx2*hyi(v)*(
     >              cxi*(-fin(3,i,j)  +fin(3,i,j+1))+
     >              cx*(-fin(3,i+1,j) +fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

c-----------------------------------
c  access to 4th derivatives

      else if(ict(1).eq.4) then
         if(ict(2).eq.1) then
c  evaluate d4f/dx3dy (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               yp=yparam(v)
               ypi=1.0_r8-yp
               yp2=yp*yp
               ypi2=ypi*ypi
               cyd=3.0_r8*yp2-1.0_r8
               cydi=-3.0_r8*ypi2+1.0_r8

               sum=hxi(v)*hyi(v)*(
     >              +( fin(1,i,j)  -fin(1,i,j+1))
     >              +(-fin(1,i+1,j)+fin(1,i+1,j+1)))

               sum=sum+sixth*hy(v)*hxi(v)*(
     >              -(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))
     >              +(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(3).eq.1) then
c  evaluate d4f/dx2dy2
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

               xp=xparam(v)
               xpi=1.0_r8-xp
               yp=yparam(v)
               ypi=1.0_r8-yp

               sum=xpi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+
     >              xp*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(4).eq.1) then
c  evaluate d4f/dxdy3 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0_r8*xp2-1.0_r8
               cxdi=-3.0_r8*xpi2+1.0_r8

               sum=hyi(v)*hxi(v)*(
     >              +( fin(2,i,j)  -fin(2,i,j+1))
     >              +(-fin(2,i+1,j)+fin(2,i+1,j+1)))

               sum=sum+sixth*hx(v)*hyi(v)*(
     >              cxdi*(-fin(3,i,j)  +fin(3,i,j+1))+
     >              cxd*(-fin(3,i+1,j) +fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

c-----------------------------------
c  access to 5th derivatives

      else if(ict(1).eq.5) then
         if(ict(2).eq.1) then
c  evaluate d5f/dx3dy2 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

               yp=yparam(v)
               ypi=1.0_r8-yp

               sum=hxi(v)*(
     >              -(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))
     >              +(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(3).eq.1) then
c  evaluate d5f/dx2dy3 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)

               xp=xparam(v)
               xpi=1.0_r8-xp

               sum=hyi(v)*(
     >              xpi*(-fin(3,i,j)  +fin(3,i,j+1))+
     >              xp*(-fin(3,i+1,j)+fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

c-----------------------------------
c  access to 6th derivatives

      else if(ict(1).eq.6) then
c  evaluate d6f/dx3dy3 (not continuous)
         iadr=iadr+1
         do v=1,ivec
            i=ii(v)
            j=jj(v)
            sum=hxi(v)*hyi(v)*(
     >              +( fin(3,i,j)  -fin(3,i,j+1))
     >              +(-fin(3,i+1,j)+fin(3,i+1,j+1)))
            fval(v,iadr)=sum
         enddo
      endif

      return
      end
C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 2d fcn
C   --vectorized-- dmc 10 Feb 1999
C    --optimized for VARIATION along x axis ONLY--

      subroutine r8fvbicubx(ict,ivec,ivecd,
     >   fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >   fin,inf2,ny)

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ny,inf2,iadr,j,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 z36th,yp,ypi,yp2,ypi2,cy,cyi,hy2,xp,xpi,xp2,xpi2,cx
      REAL*8 cxi,hx2,cxd,cxdi,cyd,cydi
!============
      integer ict(6)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)

      integer ii(ivec),jj               ! target cells (i,j)
      REAL*8 xparam(ivec),yparam
                          ! normalized displacements from (i,j) corners

      REAL*8 hx(ivec),hy                  ! grid spacing, and
      REAL*8 hxi(ivec),hyi                ! inverse grid spacing 1/(x(i+1)-x(i))
                                        ! & 1/(y(j+1)-y(j))

      REAL*8 fin(0:3,inf2,ny)             ! interpolant data (cf "evbicub")

      REAL*8 fval(ivecd,*)                ! output returned

C  for detailed description of fin, ict and fval see subroutine
C  evbicub comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.

C  note that the index inputs ii,jj and parameter inputs
C     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument

C  to use this routine in scalar mode, pass in ivec=ivecd=1

C---------------
C  Spline evaluation consists of a "mixing" of the interpolant
C  data using the linear functionals xparam, xpi = 1-xparam,
C  yparam, ypi = 1-yparam, and the cubic functionals
C  xparam**3-xparam, xpi**3-xpi, yparam**3-yparam, ypi**3-ypi ...
C  and their derivatives as needed.

      integer v
      REAL*8 sum

      REAL*8, parameter :: sixth = 0.166666666666666667_r8

C---------------
C   ...in x direction

      z36th=sixth*sixth
      iadr=0

      if(ict(1).le.2) then

C  get desired values:

         if(ict(1).eq.1) then

C  function value:

            j=jj

C   ...and in y direction

            yp=yparam
            ypi=1.0_r8-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cy=yp*(yp2-1.0_r8)
            cyi=ypi*(ypi2-1.0_r8)
            hy2=hy*hy

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cx=xp*(xp2-1.0_r8)
               cxi=xpi*(xpi2-1.0_r8)
               hx2=hx(v)*hx(v)

               sum=xpi*(ypi*fin(0,i,j)  +yp*fin(0,i,j+1))+
     >              xp*(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1))

               sum=sum+sixth*hx2*(
     >              cxi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+
     >              cx*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))

               sum=sum+sixth*hy2*(
     >              xpi*(cyi*fin(2,i,j)  +cy*fin(2,i,j+1))+
     >              xp*(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))

               sum=sum+z36th*hx2*hy2*(
     >              cxi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+
     >              cx*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(2).eq.1) then

C  df/dx:

            j=jj

C   ...and in y direction

            yp=yparam
            ypi=1.0_r8-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cy=yp*(yp2-1.0_r8)
            cyi=ypi*(ypi2-1.0_r8)
            hy2=hy*hy

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0_r8*xp2-1.0_r8
               cxdi=-3.0_r8*xpi2+1.0_r8

               sum=hxi(v)*(
     >              -(ypi*fin(0,i,j)  +yp*fin(0,i,j+1))
     >              +(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1)))

               sum=sum+sixth*hx(v)*(
     >              cxdi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+
     >              cxd*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))

               sum=sum+sixth*hxi(v)*hy2*(
     >              -(cyi*fin(2,i,j)  +cy*fin(2,i,j+1))
     >              +(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))

               sum=sum+z36th*hx(v)*hy2*(
     >              cxdi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+
     >              cxd*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(3).eq.1) then

C  df/dy:

            j=jj

C   ...and in y direction

            yp=yparam
            ypi=1.0_r8-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0_r8*yp2-1.0_r8
            cydi=-3.0_r8*ypi2+1.0_r8

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cx=xp*(xp2-1.0_r8)
               cxi=xpi*(xpi2-1.0_r8)
               hx2=hx(v)*hx(v)

               sum=hyi*(
     >              xpi*(-fin(0,i,j)  +fin(0,i,j+1))+
     >              xp*(-fin(0,i+1,j)+fin(0,i+1,j+1)))

               sum=sum+sixth*hx2*hyi*(
     >              cxi*(-fin(1,i,j)  +fin(1,i,j+1))+
     >              cx*(-fin(1,i+1,j)+fin(1,i+1,j+1)))

               sum=sum+sixth*hy*(
     >              xpi*(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1))+
     >              xp*(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))

               sum=sum+z36th*hx2*hy*(
     >              cxi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+
     >              cx*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(4).eq.1) then

C  d2f/dx2:

            j=jj

C   ...and in y direction

            yp=yparam
            ypi=1.0_r8-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cy=yp*(yp2-1.0_r8)
            cyi=ypi*(ypi2-1.0_r8)
            hy2=hy*hy

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp

               sum=(
     >              xpi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+
     >              xp*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))

               sum=sum+sixth*hy2*(
     >              xpi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+
     >              xp*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(5).eq.1) then

C  d2f/dy2:

            j=jj

C   ...and in y direction

            yp=yparam
            ypi=1.0_r8-yp

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cx=xp*(xp2-1.0_r8)
               cxi=xpi*(xpi2-1.0_r8)
               hx2=hx(v)*hx(v)

               sum=(
     >              xpi*(ypi*fin(2,i,j)  +yp*fin(2,i,j+1))+
     >              xp*(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))

               sum=sum+sixth*hx2*(
     >              cxi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+
     >              cx*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(6).eq.1) then

C  d2f/dxdy:

            j=jj

C   ...and in y direction

            yp=yparam
            ypi=1.0_r8-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0_r8*yp2-1.0_r8
            cydi=-3.0_r8*ypi2+1.0_r8

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)

C  in x direction...

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0_r8*xp2-1.0_r8
               cxdi=-3.0_r8*xpi2+1.0_r8

               sum=hxi(v)*hyi*(
     >              fin(0,i,j)  -fin(0,i,j+1)
     >              -fin(0,i+1,j)+fin(0,i+1,j+1))

               sum=sum+sixth*hx(v)*hyi*(
     >              cxdi*(-fin(1,i,j)  +fin(1,i,j+1))+
     >              cxd*(-fin(1,i+1,j)+fin(1,i+1,j+1)))

               sum=sum+sixth*hxi(v)*hy*(
     >              -(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1))
     >              +(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))

               sum=sum+z36th*hx(v)*hy*(
     >              cxdi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+
     >              cxd*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

C-------------------------------------------------

      else if(ict(1).eq.3) then
         if(ict(2).eq.1) then
c  evaluate d3f/dx3 (not continuous)
            j=jj
            yp=yparam
            ypi=1.0_r8-yp
            yp2=yp*yp
            ypi2=ypi*ypi
            cy=yp*(yp2-1.0_r8)
            cyi=ypi*(ypi2-1.0_r8)
            hy2=hy*hy

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               sum=hxi(v)*(
     >              -(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))
     >              +(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))

               sum=sum+sixth*hy2*hxi(v)*(
     >              -(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))
     >              +(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(3).eq.1) then
c  evaluate d3f/dx2dy
            j=jj

            yp=yparam
            ypi=1.0_r8-yp
            yp2=yp*yp
            ypi2=ypi*ypi
            cyd=3.0_r8*yp2-1.0_r8
            cydi=-3.0_r8*ypi2+1.0_r8

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               xp=xparam(v)
               xpi=1.0_r8-xp

               sum=hyi*(
     >              xpi*(-fin(1,i,j)  +fin(1,i,j+1))+
     >              xp*(-fin(1,i+1,j) +fin(1,i+1,j+1)))

               sum=sum+sixth*hy*(
     >              xpi*(cydi*fin(3,i,j) +cyd*fin(3,i,j+1))+
     >              xp*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(4).eq.1) then
c  evaluate d3f/dxdy2
            j=jj
            yp=yparam
            ypi=1.0_r8-yp

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi
               cxd=3.0_r8*xp2-1.0_r8
               cxdi=-3.0_r8*xpi2+1.0_r8

               sum=hxi(v)*(
     >              -(ypi*fin(2,i,j)  +yp*fin(2,i,j+1))
     >              +(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))

               sum=sum+sixth*hx(v)*(
     >              cxdi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+
     >              cxd*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(5).eq.1) then
c  evaluate d3f/dy3 (not continuous)
            iadr=iadr+1
            j=jj
            do v=1,ivec
               i=ii(v)

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cx=xp*(xp2-1.0_r8)
               cxi=xpi*(xpi2-1.0_r8)
               hx2=hx(v)*hx(v)

               sum=hyi*(
     >              xpi*(-fin(2,i,j)  +fin(2,i,j+1))+
     >              xp*(-fin(2,i+1,j) +fin(2,i+1,j+1)))

               sum=sum+sixth*hx2*hyi*(
     >              cxi*(-fin(3,i,j)  +fin(3,i,j+1))+
     >              cx*(-fin(3,i+1,j) +fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

c-----------------------------------
c  access to 4th derivatives

      else if(ict(1).eq.4) then
         if(ict(2).eq.1) then
c  evaluate d4f/dx3dy (not continuous)
            iadr=iadr+1
            j=jj

            yp=yparam
            ypi=1.0_r8-yp
            yp2=yp*yp
            ypi2=ypi*ypi
            cyd=3.0_r8*yp2-1.0_r8
            cydi=-3.0_r8*ypi2+1.0_r8

            do v=1,ivec
               i=ii(v)

               sum=hxi(v)*hyi*(
     >              +( fin(1,i,j)  -fin(1,i,j+1))
     >              +(-fin(1,i+1,j)+fin(1,i+1,j+1)))

               sum=sum+sixth*hy*hxi(v)*(
     >              -(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))
     >              +(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(3).eq.1) then
c  evaluate d4f/dx2dy2
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj

               xp=xparam(v)
               xpi=1.0_r8-xp
               yp=yparam
               ypi=1.0_r8-yp

               sum=xpi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+
     >              xp*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(4).eq.1) then
c  evaluate d4f/dxdy3 (not continuous)
            j=jj

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)

               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0_r8*xp2-1.0_r8
               cxdi=-3.0_r8*xpi2+1.0_r8

               sum=hyi*hxi(v)*(
     >              +( fin(2,i,j)  -fin(2,i,j+1))
     >              +(-fin(2,i+1,j)+fin(2,i+1,j+1)))

               sum=sum+sixth*hx(v)*hyi*(
     >              cxdi*(-fin(3,i,j)  +fin(3,i,j+1))+
     >              cxd*(-fin(3,i+1,j) +fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

c-----------------------------------
c  access to 5th derivatives

      else if(ict(1).eq.5) then
         if(ict(2).eq.1) then
c  evaluate d5f/dx3dy2 (not continuous)
            j=jj

            yp=yparam
            ypi=1.0_r8-yp

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)

               sum=hxi(v)*(
     >              -(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))
     >              +(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

         if(ict(3).eq.1) then
c  evaluate d5f/dx2dy3 (not continuous)
            j=jj
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)

               xp=xparam(v)
               xpi=1.0_r8-xp

               sum=hyi*(
     >              xpi*(-fin(3,i,j)  +fin(3,i,j+1))+
     >              xp*(-fin(3,i+1,j)+fin(3,i+1,j+1)))

               fval(v,iadr)=sum
            enddo
         endif

c-----------------------------------
c  access to 6th derivatives

      else if(ict(1).eq.6) then
c  evaluate d6f/dx3dy3 (not continuous)
         iadr=iadr+1
         j=jj
         do v=1,ivec
            i=ii(v)
            sum=hxi(v)*hyi*(
     >              +( fin(3,i,j)  -fin(3,i,j+1))
     >              +(-fin(3,i+1,j)+fin(3,i+1,j+1)))
            fval(v,iadr)=sum
         enddo
      endif

      return
      end

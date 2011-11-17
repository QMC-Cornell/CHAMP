      subroutine orbitals_loc_num_grade(ie,x,orb,dorb,ddorb)
c Written by Cyrus Umrigar and Amit Ghosal
c Calculate orbitals for finite system by 2-dim cubic spline interpolation
      use coefs_mod
      use dim_mod
      use orbital_num_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*),orb(*),dorb(3,*),ddorb(*),splineval(6)

c Warning: At present for simplicity a separate call is made to evaluate
c the splines for each electron and orbital.  However it would be more a bit
c efficient to modify the spline evaluation routine so that all the
c orbitals are evaluated in a single call.  At present spline routine has
c capability to evaluate spline for all electrons in single call, but using
c this would just save the overhead of a subroutine call.

        xi=x(1,ie)
        yi=x(2,ie)

c Select the nearest grid point in case (xi,yi) falls outside the grid.
        xget=xi
        yget=yi
        iflag=0
        if(xi.gt.sizex.or.yi.gt.sizey.or.xi.lt.-sizex.or.yi.lt.-sizey) then
          iflag=1
          if(xi.gt.sizex) xget=sizex
          if(yi.gt.sizey) yget=sizey
          if(xi.lt.-sizex) xget=-sizex
          if(yi.lt.-sizey) yget=-sizey
        endif

c cell indices
        ix=int((xget-xorb_grid(1)-1.d-12)*hxi)+1
        iy=int((yget-yorb_grid(1)-1.d-12)*hyi)+1

c normalized displacement from (xorb_grid(ix),yorb_grid(iy)) corner of cell
c  xscaled=0 @xorb_grid(ix)  xscaled=1 @xorb_grid(ix+1)
c  yscaled=0 @yorb_grid(iy)  yscaled=1 @yorb_grid(iy+1)
c       if (ix.le.0) then
c         write (6,*) 'orbitals_loc_num_grade error'
c         write (6,*) ix, xget, hxi 
c         call flush(6)
c       endif
        xscaled=(xget-xorb_grid(ix))*hxi
        yscaled=(yget-yorb_grid(iy))*hyi

        do 10 iorb=1,norb

c Evaluate spline function at (xget,yget) using the spline coeff's orb_num(1-4,1-ngrid_orbx,1-ngrid_orby,iorb)
          call r8fvbicub(ict,1,1,splineval,ix,iy,xscaled,yscaled,hx,hxi,hy,hyi
     &    ,orb_num(1,1,1,iorb),ngrid_orbx,ngrid_orby)

          if(iflag.eq.1) then

            if(abs(splineval(1)).lt.1.d-99) then
              do n=1,5
                splineval(n)=1.d-99
              enddo

             else

c If the Monte carlo point (x,y) falls outside the grid then assume that the
c wave fn. outside the grid falls off as:
c coeff0*dexp(-(coeff1*(x-coeff3)**2+coeff2*(y-coeff4)**2))
c Get the 5 coefficients by equating the function, its first derivatives and its
c diagonal second derivatives at the nearest grid point (xget,yget).
c Then evaluate above function at Monte Carlo point (xi,yi).

              coeff1=0.5d0*((splineval(2)**2/splineval(1))
     &          -splineval(4))/splineval(1)
              coeff2=0.5d0*((splineval(3)**2/splineval(1))
     &          -splineval(5))/splineval(1)
              coeff3=xget+(splineval(2)*splineval(1))/
     &          (splineval(2)**2-splineval(1)*splineval(4))
              coeff4=yget+(splineval(3)*splineval(1))/
     &          (splineval(3)**2-splineval(1)*splineval(5))
              coeff0=splineval(1)/dexp(-(coeff1*(xget-coeff3)**2
     &          +coeff2*(yget-coeff4)**2))

              splineval(1)=coeff0*dexp(-(coeff1*(xi-coeff3)**2+   !The function
     &          coeff2*(yi-coeff4)**2))                           !and its derivatives
              splineval(2)=-2.d0*coeff1*(xi-coeff3)*splineval(1)  !are evaluated at (xi,yi)
              splineval(3)=-2.d0*coeff2*(yi-coeff4)*splineval(1)
              splineval(4)=((2.d0*coeff1*(xi-coeff3))**2
     &          -2.d0*coeff1)*splineval(1)
              splineval(5)=((2.d0*coeff2*(yi-coeff4))**2
     &          -2.d0*coeff2)*splineval(1)
c             splineval(6)=4.d0*coeff1*coeff2*(xi-coeff3)
c    &          *(yi-coeff4)*splineval(1)

            endif
          endif

          orb(iorb)=splineval(1)
          ddorb(iorb)=splineval(4)+splineval(5)
            do 10 k=1,ndim
              dorb(k,iorb)=splineval(k+1)
  10         continue
      return
      end

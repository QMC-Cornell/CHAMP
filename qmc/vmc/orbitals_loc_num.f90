      subroutine orbitals_loc_num(iel,x,orb,dorb,ddorb)
! Written by Cyrus Umrigar and Amit Ghosal
! Calculate orbitals and derivatives for finite system by 2-dim cubic spline interpolation
      use coefs_mod
      use const_mod
      use dim_mod
      use orbital_num_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*),orb(nelec,*),dorb(3,nelec,*),ddorb(nelec,*),splineval(6)

! Warning: At present for simplicity a separate call is made to evaluate
! the splines for each electron and orbital.  However it would be more a bit
! efficient to modify the spline evaluation routine so that all the
! orbitals are evaluated in a single call.  At present spline routine has
! capability to evaluate spline for all electrons in single call, but using
! this would just save the overhead of a subroutine call.

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
       else
        nelec1=iel
        nelec2=iel
      endif

      do 10 ie=nelec1,nelec2

        xi=x(1,ie)
        yi=x(2,ie)

! Select the nearest grid point in case (xi,yi) falls outside the grid.
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

! cell indices
        ix=int((xget-xorb_grid(1)-1.d-12)*hxi)+1
        iy=int((yget-yorb_grid(1)-1.d-12)*hyi)+1

! normalized displacement from (xorb_grid(ix),yorb_grid(iy)) corner of cell
!  xscaled=0 @xorb_grid(ix)  xscaled=1 @xorb_grid(ix+1)
!  yscaled=0 @yorb_grid(iy)  yscaled=1 @yorb_grid(iy+1)
        xscaled=(xget-xorb_grid(ix))*hxi
        yscaled=(yget-yorb_grid(iy))*hyi

        do 10 iorb=1,norb

! Evaluate spline function at (xget,yget) using the spline coeff's orb_num(1-4,1-ngrid_orbx,1-ngrid_orby,iorb)
          call r8fvbicub(ict,1,1,splineval,ix,iy,xscaled,yscaled,hx,hxi,hy,hyi &
     &    ,orb_num(1,1,1,iorb),ngrid_orbx,ngrid_orby)

          if(iflag.eq.1) then

            if(abs(splineval(1)).lt.1.d-99) then
              do n=1,5
                splineval(n)=1.d-99
              enddo

             else

! If the Monte carlo point (xi,yi) falls outside the grid then assume that the
! wave fn. outside the grid falls off as:
! coeff0*dexp(-(coeff1*(x-coeff3)**2+coeff2*(y-coeff4)**2))
! Get the 5 coefficients by equating the function, its first derivatives and its
! diagonal second derivatives at the nearest grid point (xget,yget).
! Then evaluate above function at Monte Carlo point (xi,yi).

              coeff1=0.5d0*((splineval(2)**2/splineval(1)) &
     &          -splineval(4))/splineval(1)
              coeff2=0.5d0*((splineval(3)**2/splineval(1)) &
     &          -splineval(5))/splineval(1)
              coeff3=xget+(splineval(2)*splineval(1))/ &
     &          (splineval(2)**2-splineval(1)*splineval(4))
              coeff4=yget+(splineval(3)*splineval(1))/ &
     &          (splineval(3)**2-splineval(1)*splineval(5))
              coeff0=splineval(1)/dexp(-(coeff1*(xget-coeff3)**2 &
     &          +coeff2*(yget-coeff4)**2))

! Evaluate function and derivatives at (xi,yi)
              splineval(1)=coeff0*dexp(-(coeff1*(xi-coeff3)**2+ &
     &          coeff2*(yi-coeff4)**2))
              splineval(2)=-2.d0*coeff1*(xi-coeff3)*splineval(1)
              splineval(3)=-2.d0*coeff2*(yi-coeff4)*splineval(1)
              splineval(4)=((2.d0*coeff1*(xi-coeff3))**2 &
     &          -2.d0*coeff1)*splineval(1)
              splineval(5)=((2.d0*coeff2*(yi-coeff4))**2 &
     &          -2.d0*coeff2)*splineval(1)
!             splineval(6)=4.d0*coeff1*coeff2*(xi-coeff3)
!    &          *(yi-coeff4)*splineval(1)

            endif
          endif

          orb(ie,iorb)=splineval(1)
          ddorb(ie,iorb)=splineval(4)+splineval(5)
            do 10 k=1,ndim
   10         dorb(k,ie,iorb)=splineval(k+1)

!       do 20 ie=1,nelec
!  20   write(6,'(''orb in orb_loc_num'',i3,10d12.4)') ie,(orb(ie,iorb),iorb=1,norb)

!     write(6,'(''orb_num1a'',10d12.4)') (orb_num(1,1,1,iorb),iorb=1,norb)

      return
      end
!-----------------------------------------------------------------------

      subroutine orbitals_loc_nume(x,orb)
! Written by Cyrus Umrigar and Amit Ghosal
! Calculate orbitals for finite system by 2-dim cubic spline interpolation
      use coefs_mod
      use orbital_num_mod
      implicit real*8(a-h,o-z)

!     common /dim/ ndim
!     common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

      dimension x(3),orb(*) &
     &,splineval(6)

! Warning: At present for simplicity a separate call is made to evaluate
! the splines for each electron and orbital.  However it would be more a bit
! efficient to modify the spline evaluation routine so that all the
! orbitals are evaluated in a single call.  At present spline routine has
! capability to evaluate spline for all electrons in single call, but using
! this would just save the overhead of a subroutine call.

!     do 10 ie=nelec1,nelec2

        xi=x(1)
        yi=x(2)

! Select the nearest grid point in case (xi,yi) falls outside the grid.
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

! cell indices
        ix=int((xget-xorb_grid(1)-1.d-12)*hxi)+1
        iy=int((yget-yorb_grid(1)-1.d-12)*hyi)+1

! normalized displacement from (xorb_grid(ix),yorb_grid(iy)) corner of cell
!  xscaled=0 @xorb_grid(ix)  xscaled=1 @xorb_grid(ix+1)
!  yscaled=0 @yorb_grid(iy)  yscaled=1 @yorb_grid(iy+1)
        xscaled=(xget-xorb_grid(ix))*hxi
        yscaled=(yget-yorb_grid(iy))*hyi

        do 10 iorb=1,norb

! Evaluate spline function at (xget,yget) using the spline coeff's orb_num(1-4,1-ngrid_orbx,1-ngrid_orby,iorb)
! Warning: We need just the function, so I should reset ict before entering this loop or this routine to save some time.
          call r8fvbicub(ict,1,1,splineval,ix,iy,xscaled,yscaled,hx,hxi,hy,hyi &
     &    ,orb_num(1,1,1,iorb),ngrid_orbx,ngrid_orby)

          if(iflag.eq.1) then

            if(abs(splineval(1)).lt.1.d-99) then
              splineval(1)=1.d-99

             else

! If the Monte carlo point (xi,yi) falls outside the grid then assume that the
! wave fn. outside the grid falls off as:
! coeff0*dexp(-(coeff1*(x-coeff3)**2+coeff2*(y-coeff4)**2))
! Get the 5 coefficients by equating the function, its first derivatives and its
! diagonal second derivatives at the nearest grid point (xget,yget).
! Then evaluate above function at Monte Carlo point (xi,yi).

              coeff1=0.5d0*((splineval(2)**2/splineval(1)) &
     &          -splineval(4))/splineval(1)
              coeff2=0.5d0*((splineval(3)**2/splineval(1)) &
     &          -splineval(5))/splineval(1)
              coeff3=xget+(splineval(2)*splineval(1))/ &
     &          (splineval(2)**2-splineval(1)*splineval(4))
              coeff4=yget+(splineval(3)*splineval(1))/ &
     &          (splineval(3)**2-splineval(1)*splineval(5))
              coeff0=splineval(1)/dexp(-(coeff1*(xget-coeff3)**2 &
     &          +coeff2*(yget-coeff4)**2))

! Evaluate function at (xi,yi)
              splineval(1)=coeff0*dexp(-(coeff1*(xi-coeff3)**2+ &
     &          coeff2*(yi-coeff4)**2))

            endif
          endif

   10     orb(iorb)=splineval(1)

!    write(6,'(''orb in orb_loc_num'',i3,10d12.4)') (orb(iorb),iorb=1,norb)
!     write(6,'(''orb_num1a'',10d12.4)') (orb_num(1,1,1,iorb),iorb=1,norb)

      return
      end

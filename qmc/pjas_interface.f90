!!!! this defines various wrappers to call various f90 functions used in the f77 code.
!!!! should be removed when the f77 is brought up to date.


subroutine pjas_jas_interface (x,rvec, v,d2,div_vj,value)
  use periodic_jastrow_mod
  use pjasen_mod
  use pjasee_mod
  include 'modules.h'
  implicit none

  real(dp)                               :: rvec(3,nelec*(nelec-1)/2)
  real (dp)                              :: x, v, div_vj, value, d2
  dimension x(3,nelec),v(3,nelec),div_vj(nelec)
  integer                                :: i, j, k

!!! all values will be added up (v, d2, div_vj, value)
  call pjas_jas (x,rvec,v,d2,div_vj,value)

  do i=1, nelec
     do k=1, ndim
        fjo (k, i)= v (k, i )
        fijo (k,i,i) =  fijo (k, i, i) + pjasfijo (k, i, i)
     enddo
     fso (i, i) = fso (i, i ) + pjasfso (i, i)
     d2ijo (i,i) = d2ijo (i, i) + pjasd2ijo (i,i)
  enddo

  d2o = d2o + pjasd2o

  !  write(47,*) "fsumo , pjasfsumo", fsumo ,  pjasfsumo
  fsumo = fsumo + pjasfsumo

  do i=2,nelec
     do j=1,i-1
        fso ( i, j ) = fso (i, j) + pjasfso (i, j)
        fso ( j, i ) = fso (j, i) + pjasfso (j, i)
        d2ijo (i,j) = d2ijo (i, j) + pjasd2ijo (i,j)
        d2ijo (j,i) = d2ijo (j, i) + pjasd2ijo (j,i)
        do k=1, ndim
           fijo (k,i,j) =  fijo (k, i, j) + pjasfijo (k, i, j)
           fijo (k,j,i) =  fijo (k, j, i) + pjasfijo (k, j, i)
        enddo
     enddo
  enddo

end subroutine pjas_jas_interface
!-------------------------------------------------------------------------

subroutine  pjas_deriv_jas_interface (x, rvec, v,d2,div_vj,value)
  use periodic_jastrow_mod
  use pjasen_mod
  use pjasee_mod
  include 'modules.h'
  implicit none
  real(dp)                               :: rvec(3,nelec*(nelec-1)/2)
  real (dp)                              :: x, v, value, d2 ,div_vj(nelec)
  dimension x(3,nelec),v(3,nelec)

  integer                                :: i, j, k


!!! all values will be added up (v, d2, div_vj, value)
  call pjas_deriv_jas (x,rvec,v,d2,div_vj, value)

  do i=1, nelec
     do k=1, ndim
        fjo (k, i)= v (k, i )
        fijo (k,i,i) =  fijo (k, i, i) + pjasfijo (k, i, i)
     enddo
     fso (i, i) = fso (i, i ) + pjasfso (i, i)
     d2ijo (i,i) = d2ijo (i, i) + pjasd2ijo (i,i)
  enddo

  d2o = d2o + pjasd2o

  !  write(47,*) "fsumo , pjasfsumo", fsumo ,  pjasfsumo

  fsumo = fsumo + pjasfsumo

  do i=2,nelec
     do j=1,i-1
        fso ( i, j ) = fso (i, j) + pjasfso (i, j)
        fso ( j, i ) = fso (j, i) + pjasfso (j, i)
        d2ijo (i,j) = d2ijo (i, j) + pjasd2ijo (i,j)
        d2ijo (j,i) = d2ijo (j, i) + pjasd2ijo (j,i)
        do k=1, ndim
           fijo (k,i,j) =  fijo (k, i, j) + pjasfijo (k, i, j)
           fijo (k,j,i) =  fijo (k, j, i) + pjasfijo (k, j, i)
        enddo
     enddo
  enddo

end subroutine pjas_deriv_jas_interface
!-------------------------------------------------------------------------

subroutine  pjas_jas_e_interface (iel, x, rvec, v, value)
  use periodic_jastrow_mod
  use pjasen_mod
  use pjasee_mod
  include 'modules.h'
  implicit none
  real(dp)                               :: rvec(3,nelec*(nelec-1)/2)
  integer                                :: iel
  real (dp)                              :: x, v, value
  dimension x(3,nelec),v(3,nelec)
!JT  real (dp)                              :: fsn, fijn, d2ijn, d2n , fsumn, fjn

  call pjas_jas_e (iel, x, rvec, v, value)

end subroutine pjas_jas_e_interface

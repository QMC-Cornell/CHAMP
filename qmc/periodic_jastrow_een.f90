subroutine f_een_cuts (rp, ri, rj, gfunci, gfuncj, gfunc, dgfunci, dgfuncj, d2gfunci,d2gfuncj)
  !---------------------------------------------------------------------------
  ! Description : evaluates the cutoff function used in the een term of the Jastrow function
  !
  ! used to add g (ri) g (rj) to een Jastrow which is used
  ! in periodic calculations
  !
  !!  g(ri)= 1-6 * (r/rp)**2 + 8 * (r/rp)**3 - 3 *(r/rp)**4 for  r < rp
  !!       =  0             r >= rp
  !!
  ! Created     : W. A. Al-Saidi, June 2007
  !---------------------------------------------------------------------------
  implicit none
  integer, parameter                     :: dp = kind(1.0d0)
  real(dp)                               :: ri,rj,rp,rii,rjj,rp2
  real(dp)                               :: ri2, ri3, ri4
  real(dp)                               :: rj2, rj3, rj4
  real(dp)                               :: gfunc, gfunci, gfuncj
  real(dp)                               :: dgfunci , dgfuncj , d2gfunci,d2gfuncj

  rp2 = rp*rp
  rii= ri/rp
  ri2= rii * rii
  ri3 = ri2 * rii
  ri4 = ri3 * rii

  gfunci = 1 - 6 * ri2 + 8 * ri3 - 3 *ri4
  dgfunci = (-12 * rii + 24 * ri2 -12 * ri3)/rp
  d2gfunci = (-12 + 48 * rii - 36 * ri2)/rp2

  rjj=rj/rp
  rj2= rjj * rjj
  rj3 = rj2 * rjj
  rj4 = rj3 * rjj
  gfuncj = 1 -6 * rj2 + 8 * rj3 - 3 *rj4
  dgfuncj = (-12 * rjj + 24 * rj2 -12 * rj3)/rp
  d2gfuncj = (-12 + 48 * rjj - 36 * rj2)/rp2

  gfunc = gfunci * gfuncj

!!$    write(45,*) "ri rj ", ri,rj,rp
!!$    write(45,*) "  gfunci  gfuncj ", gfunci , gfuncj
!!$    write(45,*) "  dgfunci  dgfuncj ", dgfunci , dgfuncj
!!$    write(45,*) "  d2gfunci  d2gfuncj ", d2gfunci , d2gfuncj

end subroutine f_een_cuts


subroutine f_een_cuts_nd (rp, ri, rj, gfunc)
  !---------------------------------------------------------------------------
  ! Description : same as f_een_cuts but calculate no derivatives
  !
  ! Created     : W. A. Al-Saidi, June 2007
  !---------------------------------------------------------------------------
  implicit none
  integer, parameter                     :: dp = kind(1.0d0)
  real(dp)                               :: ri,rj,rp,rii,rjj,rp2
  real(dp)                               :: ri2, ri3, ri4
  real(dp)                               :: rj2, rj3, rj4
  real(dp)                               :: gfunc, gfunci, gfuncj
!  real(dp)                               :: dgfunci, d2gfuncj, d2gfunci, dgfuncj

  rp2 = rp*rp
  rii= ri/rp
  ri2= rii * rii
  ri3 = ri2 * rii
  ri4 = ri3 * rii

  gfunci = 1 - 6 * ri2 + 8 * ri3 - 3 *ri4

  rjj=rj/rp
  rj2= rjj * rjj
  rj3 = rj2 * rjj
  rj4 = rj3 * rjj
  gfuncj = 1 -6 * rj2 + 8 * rj3 - 3 *rj4

  gfunc = gfunci * gfuncj

!!$    write(45,*) "ri rj ", ri,rj,rp
!!$    write(45,*) "  gfunci  gfuncj ", gfunci , gfuncj
!!$    write(45,*) "  dgfunci  dgfuncj ", dgfunci , dgfuncj
!!$    write(45,*) "  d2gfunci  d2gfuncj ", d2gfunci , d2gfuncj

end subroutine f_een_cuts_nd




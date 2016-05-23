subroutine psig(psit_l,velocity,psign_l,d)!,velocity_g)
  use dim_mod, only: ndim
  use const_mod, only: nelec
  use config_mod, only: psig_epsilon
  implicit none 
  real*8,intent(in) :: psit_l,velocity(ndim,nelec)
  real*8, intent(out) :: psign_l
  real*8 f,d,mag_vel2
  integer i,j

!Temporaray calculation of d for printout - should remove later
  mag_vel2 = 0
  do j=1,nelec
     do i=1,ndim
        mag_vel2 = mag_vel2 + velocity(i,j)**2
     enddo
  enddo

  d = 1/sqrt(mag_vel2)

  if(psig_epsilon.gt.0) then
     mag_vel2 = 0
     do j=1,nelec
        do i=1,ndim
           mag_vel2 = mag_vel2 + velocity(i,j)**2
        enddo
     enddo
     
     d = 1/sqrt(mag_vel2)

!  write(9,'(''d,f'',3es12.5)') d,f,psit_l
!  write(8,*)d

     if(d<psig_epsilon) then
        !     write(9,*) 'd<eps'
        f = 0.5*(psig_epsilon/d + d/psig_epsilon)
        !We want f*psit_l, but we're working with log((f*psit_l)**2), so
        !we add 2*log(f)
        psign_l = psit_l+2*log(f)
        ! do j=1,nelec
        !    do i=1,ndim
        !       velocity_g(i,j) = velocity(i,j)*d**2/epsilon
        !    enddo
        ! enddo
     else
        psign_l = psit_l
     endif
  else
     psign_l = psit_l
  endif

end subroutine psig

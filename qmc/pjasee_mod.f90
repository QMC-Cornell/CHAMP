module pjasee_mod

  use all_tools_mod
  use crystal_symmetries_mod
  use eloc_mod
  use pjasen_mod


  implicit none

  real(dp), allocatable                  :: lap_dpsi_pjasee(:,:)
  real(dp), allocatable                  :: grad_dpsi_pjasee (:,:,:)
  real(dp), allocatable                  :: cos_star_ee (:,:)
  real(dp), allocatable                  :: grad_cos_star_ee (:,:,:)


  integer                                :: param_pjasee_nb
  real(dp), allocatable                  :: dpsi_pjasee (:)
  real(dp), allocatable                  :: deloc_pjasee (:)


  real(dp), allocatable                  :: pjasee_parms (:,:)
  real(dp), allocatable                  :: gn_pjasee  (:)

  real (dp), allocatable                 :: cos_pjasee(:,:), sin_pjasee(:,:)
  real (dp), allocatable                 :: cos_pjasee_e (:,:), sin_pjasee_e (:,:)



contains




  subroutine deloc_pjasee_bld
!---------------------------------------------------------------------------
! Description : derivative of local energy wrt periodic jastrow parameters
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
    integer                              :: k,ieta, ie
    real(dp)                             :: sum
!JT    include 'common_vd.h'
    ! header
    if (header_exe) then
       call object_create ('deloc_pjasee',deloc_pjasee_index)
       call object_needed ('nelec')
       call object_needed ('ndim')
       call object_needed ('vd')
       call object_needed ('vj')
       call object_needed ('param_pjasee_nb')
       call object_needed ('grad_dpsi_pjasee')
       call object_needed ('lap_dpsi_pjasee')
    endif

    if (nforce==3) stop "deloc_pjasee is called"

    call object_alloc ('deloc_pjasee', deloc_pjasee, param_pjasee_nb)

    do ieta=1,param_pjasee_nb

       sum = 0
       do ie = 1 , nelec
          sum  = sum + lap_dpsi_pjasee (ie, ieta)

          do k= 1, ndim_pj
             sum  = sum + 2* grad_dpsi_pjasee (k, ie, ieta) * (vd (k, ie)+vj(k,ie))
          enddo

       enddo

       deloc_pjasee (ieta) =-hb * sum

    end do

  end subroutine deloc_pjasee_bld



  subroutine  star_ee_bld  (xvec)

!---------------------------------------------------------------------------
! Description :   building the arrays related to the ee stars
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none


    real(dp)                             :: xvec (3, nelec*(nelec-1)/2)
    integer                              :: ie, ist , ije, je

    call object_alloc ("grad_cos_star_ee ", grad_cos_star_ee, ndim_pj, nelec * (nelec-1)/2,  nstar_ee)

    call object_alloc ("cos_star_ee ", cos_star_ee , nelec * (nelec-1)/2,  nstar_ee)


    do ist= 1,  nstar_ee

       ije = 0

       do ie = 2, nelec

          do je = 1, ie-1

             ije = ije +1

             call star_ee_fac ( ije, ist, cos_star_ee ( ije, ist ),  grad_cos_star_ee (:, ije, ist ))

          enddo

       enddo

    enddo

    call object_modified_by_index (cos_star_ee_index)

    call object_modified_by_index (grad_cos_star_ee_index)


  end subroutine star_ee_bld



  subroutine deriv_nonloc_pjas_ee ( iel, xcoord, value1) !, gn_pjas)
!---------------------------------------------------------------------------
! Description : evaluates value1 = sum_{star} c_star rho(star) for electron iel
!       also creates gn_pjas = d value1/dc_star= rho (s)  for electron iel
!
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
!JT    include "common_jasn.h"
    integer                              :: iel, ist, ist_s, i, j, jj, ij
    real(dp)                             :: xcoord (3,nelec)
    real(dp)                             :: value1, temp , cosstar

    !    call  find_rvec_ee (xcoord, rvec)

    call object_alloc ("gn_pjasee", gn_pjasee, param_pjasee_nb)

    do ist = 1, param_pjasee_nb
       gn_pjasee (ist) = dpsi_pjasee (ist)
    end do

    value1 = 0

    !! find list of cos and sin
    call calc_cos_sin_ee_e (iel, xcoord)

    !! add only updated electron. The old falue is already subtracted
    !! through fso

    do ist = 1, param_pjasee_nb

       ist_s = ist + param_pjasen_nb

       do jj=1,nelec

          if(jj.eq.iel) cycle

          if(jj.lt.iel) then
             i=iel
             j=jj
          else
             i=jj
             j=iel
          endif

          ij=((i-1)*(i-2))/2+j

          !          xvec = find_rvec_ee (xcoord(:,i), xcoord(:,j))

          cosstar= cos_star_sim_fac  (ist, ij)

          temp = pjas_parms(ist_s, iwf) *  cosstar  !! only changed electron

          value1 = value1 + temp

          !          fsn (i,j ) = fsn (i,j) + temp

          gn_pjasee (ist) =  gn_pjasee (ist) + cosstar - peego  ( i, j, ist ) !! subtract the old value

       enddo
    enddo

  end subroutine deriv_nonloc_pjas_ee



  subroutine nonloc_pjas_ee (iel, xcoord, value1)
!---------------------------------------------------------------------------
! Description :evaluates value1 = sum_{star} c_star rho(star) for electron iel
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
!JT    include "common_jasn.h"
    integer                              :: iel, ist, ist_s, i, j, jj, ij
    real(dp)                             :: xcoord (3,nelec)
    real(dp)                             :: value1, temp


    value1 = 0


    !! find list of cos and sin
    call calc_cos_sin_ee_e (iel, xcoord)

    do ist = 1, param_pjasee_nb

       ist_s = ist + param_pjasen_nb

       do jj=1,nelec

          if(jj.eq.iel) cycle

          if(jj.lt.iel) then
             i=iel
             j=jj
          else
             i=jj
             j=iel
          endif

          ij=((i-1)*(i-2))/2+j

          !          xvec = find_rvec_ee (xcoord(:,i), xcoord(:,j))

          temp = pjas_parms(ist_s, iwf) *  cos_star_sim_fac  (ist, ij) !! only changed electron

          !! note: old value is already subtracted with fsno (i,j)
          value1 = value1 + temp

          !          fsn (i, j ) = fsn(i,j) + temp

       enddo

    enddo

  end subroutine nonloc_pjas_ee




!!$  subroutine find_rvec_ee (x, rvecee)
!---------------------------------------------------------------------------
! Description :
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
!!$    implicit  none
!!$
!!$    real (dp)                         :: rvecee(3,nelec*(nelec-1)/2)
!!$
!!$    real(dp)                          :: x(3,nelec), norm
!!$
!!$    integer ij, i,j, k
!!$
!!$    ij=0
!!$    do i=1,nelec
!!$       do j=1,i-1
!!$          ij=ij+1
!!$          do k=1,ndim_pj
!!$             rvecee(k,ij)=x(k,i)-x(k,j)
!!$          enddo
!!$          call find_image3(rvecee(1,ij),norm ,rlatt_sim ,rlatt_sim_inv)
!!$       enddo
!!$    enddo
!!$  end subroutine find_rvec_ee
!!$
!!$

  function find_rvec_ee (x1, x2) result(rvecee)
!---------------------------------------------------------------------------
! Description : compute the vector difference  between electrons x1 and x2
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit  none
    real(dp)                             :: x1(3), x2(3)
    real(dp)                             :: rvecee(3), norm
    integer                              :: k

    do k=1,ndim_pj

       rvecee(k)=x1(k)-x2(k)

    enddo

    call find_image3(rvecee, norm ,rlatt_sim ,rlatt_sim_inv)

  end function find_rvec_ee



  subroutine star_ee_fac (ije, ist,  c_s_fac, grad_c_s_fac)
!---------------------------------------------------------------------------
! Description : building the arrays related to the ee stars
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
    integer                              :: ije
    real(dp)                             :: c_s_fac, grad_c_s_fac (ndim_pj)
    integer                              :: ist
    real(dp)                             :: sindot, cosdot
    integer                              :: i, i1, i2, k

    !! note now basis contains both + and -
    !! sum over the first half only
    i1 = istar_sim (ist)
    i2 = fstar_sim (ist)
    !     i2 = i1+ (mstar_sim (ist))/2 - 1

    grad_c_s_fac  = 0
    c_s_fac = 0

    do i = i1, i2

       cosdot = cos_pjasee (i, ije)
       sindot = sin_pjasee (i, ije)

       c_s_fac =  c_s_fac +  cosdot !! 2 absorbed in def of parms

       do k= 1, ndim_pj
          grad_c_s_fac (k) = grad_c_s_fac (k) - rkv_sim (k,i)* sindot
       enddo

    enddo

  end subroutine star_ee_fac



  function cos_star_sim_fac (ist,ij) result (sum)
!---------------------------------------------------------------------------
! Description : evaluates sum_{star_i,k_star} cos(k_star\dot r)
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none

    integer                              :: ist
    real(dp)                             ::  sum
    integer                              :: i, i1, i2, ij


    !! note now basis contains both + and -
    !! sum over the first half only
    i1 = istar_sim (ist)
    i2 = fstar_sim (ist)
    !     i2 = i1+ (mstar_sim (ist))/2 - 1

    sum  = 0

    do i = i1, i2

       sum = sum + cos_pjasee_e (i,ij)

    enddo

  end function  cos_star_sim_fac


  function grad_cos_star_sim_fac (ist, ij) result (sum)
!---------------------------------------------------------------------------
! Description : evaluates the gradient and the value fo the cosine terms
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
    integer                              :: ist
    real(dp)                             :: sum (ndim_pj)
    integer                              :: i, i1, i2, k , ij

    !! note now basis contains both + and -
    !! sum over the first half only
    i1 = istar_sim (ist)
    i2 = fstar_sim (ist)
    !     i2 = i1+ (mstar_sim (ist))/2 - 1

    sum  = 0

    do i = i1, i2

       do k= 1, ndim_pj
          sum (k) = sum (k) - rkv_sim(k,i)* sin_pjasee_e (i,ij)
       enddo

    enddo

  end function  grad_cos_star_sim_fac




  subroutine  pjasee_deriv_jas (rvec, pjasv, pjasd2, pjasdiv_vj, pjasfsum)
!---------------------------------------------------------------------------
! Description : computes the gradient of the wavefunction with respect to the periodic jastrow
!               parameters, and  the graident and laplacian of this quantity which is used in the
!               computation of the derivates of the local energy with respect to these jastrow parameters
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
    real (dp)                            :: rvec (3, nelec*(nelec-1)/2)
    real (dp)                            :: pjasv (3, nelec)
    real (dp)                            :: pjasfsum, pjasd2
    real (dp)                            :: pjasdiv_vj (nelec)
    real (dp)                            :: gen, fen, gradl (3), lap2, pgradl (3)
    integer                              :: ist, i, ij , j , ist_s

    if( nelec .lt. 2) return

    dpsi_pjasee = 0
    grad_dpsi_pjasee  = 0
    lap_dpsi_pjasee = 0


    !! build cosine and sine table
    call calc_cos_sin_ee (rvec)

    call star_ee_bld (rvec)

    ij=0
    do i = 2,nelec
       do j = 1, i-1

          ij=ij+1

          do ist = 1, param_pjasee_nb

             ist_s = ist + param_pjasen_nb

             fen = pjas_parms (ist_s, iwf) * cos_star_ee (ij, ist)

             gradl =  grad_cos_star_ee (:, ij, ist )

             pjasfso (i,j) = pjasfso(i, j) + fen

             pgradl = pjas_parms (ist_s, iwf) * gradl

             pjasfijo(:, i, j) = pjasfijo (: , i, j) + pgradl

             pjasfijo(:, j, i) = pjasfijo (: , j, i) - pgradl

             lap2 = - sk3_sim (ist) * cos_star_ee (ij, ist )

             pjasd2ijo (i,j) = pjasd2ijo (i,j)  + 2* pjas_parms (ist_s, iwf) * lap2

             gen = cos_star_ee (ij, ist)

             peego(i,j,ist) = peego(i,j,ist) + gen

             dpsi_pjasee (ist) = dpsi_pjasee (ist) + gen

             grad_dpsi_pjasee (:,i,ist) = grad_dpsi_pjasee (:, i,ist) +  gradl

             grad_dpsi_pjasee (:,j,ist) = grad_dpsi_pjasee (:, j,ist) -  gradl

             lap_dpsi_pjasee (i, ist)  = lap_dpsi_pjasee (i, ist) + lap2

             lap_dpsi_pjasee (j, ist)  = lap_dpsi_pjasee (j, ist) + lap2

          enddo

          pjasfsum = pjasfsum + pjasfso(i,j)

          pjasv(:,i) = pjasv(:,i) + pjasfijo(:,i,j)

          pjasv(:,j) = pjasv(:,j) + pjasfijo(:,j,i)

          pjasd2 = pjasd2 + pjasd2ijo(i,j)

          pjasdiv_vj (i)= pjasdiv_vj (i) + pjasd2ijo(i,j)/2

          pjasdiv_vj (j)= pjasdiv_vj (j) + pjasd2ijo(i,j)/2

       enddo

    end do

    call object_modified_by_index (dpsi_pjasee_index)
    call object_modified_by_index (grad_dpsi_pjasee_index)
    call object_modified_by_index (lap_dpsi_pjasee_index)

  end subroutine pjasee_deriv_jas


  subroutine  pjasee_jas (rvec, pjasv, pjasd2, pjasdiv_vj, pjasfsum)
!---------------------------------------------------------------------------
! Description :  same as pjasen_deriv_jas  but without the derivatives of the
!              \nabla J with respect periodic Jastrow parameters.
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
    real (dp)                            :: rvec (3, nelec*(nelec-1)/2)
    real (dp)                            :: pjasv (3, nelec)
    real (dp)                            :: pjasfsum, pjasd2
    real (dp)                            :: pjasdiv_vj (nelec)
    real (dp)                            :: fen, gradl (3), pgradl (3)
    integer                              :: ist, i, ij , j , ist_s


    if( nelec .lt. 2) return

    do ist = 1, param_pjasee_nb
       dpsi_pjasee(ist)=0
    enddo

    !! build cosine and sine table
    call calc_cos_sin_ee (rvec)

    call star_ee_bld (rvec)


    ij=0
    do i = 2,nelec
       do j = 1, i-1

          ij=ij+1

          do ist = 1, param_pjasee_nb

             ist_s = ist + param_pjasen_nb

             fen = pjas_parms (ist_s, iwf)* cos_star_ee (ij, ist)

             gradl =  grad_cos_star_ee (:, ij, ist )

             pjasfso (i,j) = pjasfso(i, j) + fen

             pgradl = pjas_parms (ist_s, iwf) * gradl

             pjasfijo(:, i, j) = pjasfijo (: , i, j) + pgradl

             pjasfijo(:, j, i) = pjasfijo (: , j, i) -  pgradl

             pjasd2ijo (i,j) = pjasd2ijo (i,j)  -  2* fen * sk3_sim (ist)

             dpsi_pjasee (ist) = dpsi_pjasee (ist) +  cos_star_ee (ij, ist)

          enddo

          pjasfsum = pjasfsum + pjasfso(i,j)

          pjasv(:,i) = pjasv(:,i) + pjasfijo(:,i,j)

          pjasv(:,j) = pjasv(:,j) + pjasfijo(:,j,i)

          pjasd2 = pjasd2 + pjasd2ijo(i,j)

          pjasdiv_vj (i)= pjasdiv_vj (i) + pjasd2ijo(i,j)/2

          pjasdiv_vj (j)= pjasdiv_vj (j) + pjasd2ijo(i,j)/2

       enddo

    end do

    !! needed in nonloc
    call object_modified_by_index (dpsi_pjasee_index)

  end subroutine pjasee_jas




  subroutine  pjasee_jas_e (iel, rvec, fsum)
!---------------------------------------------------------------------------
! Description :Same as  pjasee_jas  but when only one electron is modified
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
!JT    include 'common_jasn.h'
    integer                              :: iel
    real (dp)                            :: rvec (3, nelec*(nelec-1)/2)
    real (dp)                            :: fen , fen1(3), temp (3), fsum
    integer                              :: ist, i,j, ij , jj , ist_s

    if( nelec .lt. 2) return

    !! build cosine and sine table
    call calc_cos_sin_ee_e2 (iel, rvec)

    do jj = 1, nelec

       if(jj .eq. iel) cycle

       if(jj .lt. iel) then
          i=iel
          j=jj
       else
          i=jj
          j=iel
       endif

       ij=((i-1)*(i-2))/2+j

       fen = 0

       fen1= 0

       do ist = 1, param_pjasee_nb

          ist_s = ist + param_pjasen_nb

          fen = fen + pjas_parms (ist_s, iwf)* cos_star_sim_fac  (ist,ij)

          temp = grad_cos_star_sim_fac  (ist, ij)

          fen1  = fen1 +  pjas_parms (ist_s, iwf) * temp

       enddo

       fsn (i,j) = fsn (i, j) + fen

       fijn(:,i,j)=fijn(:,i,j) + fen1

       fijn(:,j,i)=fijn(:,j,i) - fen1

       fjn(:,i)= fjn(:,i) + fen1

       fjn(:,j)= fjn(:,j) - fen1

       fsum = fsum + fen

    enddo

  end subroutine pjasee_jas_e



  subroutine  calc_cos_sin_ee_e (iel, rvec)
!---------------------------------------------------------------------------
! Description : calculate the cosine and sine of the star terms when one electron is displaced.
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
    integer                              :: iel
    integer                              :: i, k, n , jj, u, v, ij
    real(dp)                             :: rvec (3, nelec), xvec (3)
    real(dp)                             :: dot, cos_tmp, sin_tmp
    real (dp), allocatable               :: cos_temp(:,:,:),  sin_temp(:,:,:)


    call object_alloc ('cos_pjasee_e',  cos_pjasee_e, nbasis_pw_sim, nelec*nelec/2)
    call object_alloc ('sin_pjasee_e',  sin_pjasee_e, nbasis_pw_sim, nelec*nelec/2)


    allocate( cos_temp (-maxkv_sim:maxkv_sim,ndim, nelec*nelec/2 ))
    allocate( sin_temp (-maxkv_sim:maxkv_sim,ndim, nelec*nelec/2))


    do jj = 1, nelec

       if(jj .eq. iel) cycle

       if(jj .lt. iel) then
          u=iel
          v=jj
       else
          u=jj
          v=iel
       endif

       ij=((u-1)*(u-2))/2+v

       xvec = find_rvec_ee (rvec(:,u), rvec(:,v))

       do i=1,ndim

          dot=0
          do k=1,ndim
             dot = dot + glatt_sim(k,i)* xvec(k)
          enddo

          cos_temp(1,i,ij)=cos(dot)
          sin_temp(1,i,ij)=sin(dot)
          cos_temp(-1,i,ij)=cos_temp(1,i,ij)
          sin_temp(-1,i,ij)=-sin_temp(1,i,ij)
          cos_temp(0,i,ij)=1.d0
          sin_temp(0,i,ij)=0.d0

          do n=2, kvec_pjasee (i)

             cos_temp(n,i, ij)=cos_temp(n-1,i, ij)*cos_temp(1,i, ij)- &
                  & sin_temp(n-1,i, ij)*sin_temp(1,i, ij)

             sin_temp(n,i, ij)=sin_temp(n-1,i, ij)*cos_temp(1,i, ij)+ &
                  & cos_temp(n-1,i, ij)*sin_temp(1,i, ij)

             cos_temp(-n,i, ij)=cos_temp(n,i, ij)

             sin_temp(-n,i, ij)=-sin_temp(n,i, ij)

          enddo

       enddo

    enddo



    do jj = 1, nelec

       if(jj .eq. iel) cycle

       if(jj .lt. iel) then
          u=iel
          v=jj
       else
          u=jj
          v=iel
       endif

       ij=((u-1)*(u-2))/2+v


       do i = 1 , nbasis_pw_sim

          cos_tmp=cos_temp(kv_sim(1,i),1, ij)*cos_temp(kv_sim(2,i),2, ij) &
               &           -sin_temp(kv_sim(1,i),1, ij)*sin_temp(kv_sim(2,i),2, ij)

          sin_tmp=sin_temp(kv_sim(1,i),1, ij)*cos_temp(kv_sim(2,i),2, ij) &
               &           +cos_temp(kv_sim(1,i),1, ij)*sin_temp(kv_sim(2,i),2, ij)

          cos_pjasee_e (i, ij)= cos_tmp*cos_temp(kv_sim(3,i),3, ij) &
               &               -sin_tmp*sin_temp(kv_sim(3,i),3, ij)

          sin_pjasee_e (i, ij)= sin_tmp*cos_temp(kv_sim(3,i),3, ij) &
               &               +cos_tmp*sin_temp(kv_sim(3,i),3, ij)
       enddo

    enddo


    deallocate (cos_temp, sin_temp)

  end subroutine calc_cos_sin_ee_e



  subroutine  calc_cos_sin_ee_e2 ( iel, rvec)
!---------------------------------------------------------------------------
! Description :This is the same as calc_cos_sin_ee except that the interpartilce distance is known
!               so to gain some efficiency I introduced a new sub to handle this case.
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
    integer                              :: iel
    integer                              :: i, k, n , jj, u, v, ij
    real(dp)                             :: rvec (3, nelec*(nelec-1)/2), xvec (3)
    real(dp)                             :: dot, cos_tmp, sin_tmp
    real (dp), allocatable               :: cos_temp(:,:,:),  sin_temp(:,:,:)



    call object_alloc ('cos_pjasee_e',  cos_pjasee_e, nbasis_pw_sim, nelec*nelec/2)
    call object_alloc ('sin_pjasee_e',  sin_pjasee_e, nbasis_pw_sim, nelec*nelec/2)


    allocate( cos_temp (-maxkv_sim:maxkv_sim,ndim, nelec*nelec/2 ))
    allocate( sin_temp (-maxkv_sim:maxkv_sim,ndim, nelec*nelec/2))


    do jj = 1, nelec

       if(jj .eq. iel) cycle

       if(jj .lt. iel) then
          u=iel
          v=jj
       else
          u=jj
          v=iel
       endif

       ij=((u-1)*(u-2))/2+v

       xvec = rvec (:, ij) !find_rvec_ee (xcoord(:,u), xcoord(:,v))

       do i=1,ndim

          dot=0
          do k=1,ndim
             dot = dot + glatt_sim(k,i)* xvec(k)
          enddo

          cos_temp(1,i,ij)=cos(dot)
          sin_temp(1,i,ij)=sin(dot)
          cos_temp(-1,i,ij)=cos_temp(1,i,ij)
          sin_temp(-1,i,ij)=-sin_temp(1,i,ij)
          cos_temp(0,i,ij)=1.d0
          sin_temp(0,i,ij)=0.d0

          do n=2, kvec_pjasee (i)

             cos_temp(n,i, ij)=cos_temp(n-1,i, ij)*cos_temp(1,i, ij)- &
                  & sin_temp(n-1,i, ij)*sin_temp(1,i, ij)

             sin_temp(n,i, ij)=sin_temp(n-1,i, ij)*cos_temp(1,i, ij)+ &
                  & cos_temp(n-1,i, ij)*sin_temp(1,i, ij)

             cos_temp(-n,i, ij)=cos_temp(n,i, ij)

             sin_temp(-n,i, ij)=-sin_temp(n,i, ij)

          enddo

       enddo

    enddo



    do jj = 1, nelec

       if(jj .eq. iel) cycle

       if(jj .lt. iel) then
          u=iel
          v=jj
       else
          u=jj
          v=iel
       endif

       ij=((u-1)*(u-2))/2+v


       do i = 1 , nbasis_pw_sim

          cos_tmp=cos_temp(kv_sim(1,i),1, ij)*cos_temp(kv_sim(2,i),2, ij) &
               &           -sin_temp(kv_sim(1,i),1, ij)*sin_temp(kv_sim(2,i),2, ij)

          sin_tmp=sin_temp(kv_sim(1,i),1, ij)*cos_temp(kv_sim(2,i),2, ij) &
               &           +cos_temp(kv_sim(1,i),1, ij)*sin_temp(kv_sim(2,i),2, ij)

          cos_pjasee_e (i, ij)= cos_tmp*cos_temp(kv_sim(3,i),3, ij) &
               &               -sin_tmp*sin_temp(kv_sim(3,i),3, ij)

          sin_pjasee_e (i, ij)= sin_tmp*cos_temp(kv_sim(3,i),3, ij) &
               &               +cos_tmp*sin_temp(kv_sim(3,i),3, ij)
       enddo
    enddo

    deallocate (cos_temp, sin_temp)

  end subroutine calc_cos_sin_ee_e2



  subroutine  calc_cos_sin_ee (xvec)
!---------------------------------------------------------------------------
! Description : cosine and sine terms
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    include 'modules.h'
    implicit none
    integer                              :: ie, i, k, n, je, ije
    real(dp)                             :: xvec (3,  nelec * (nelec-1)/2 )
    real(dp)                             :: dot, cos_tmp, sin_tmp
    real (dp), allocatable               :: cos_temp(:,:,:),  sin_temp(:,:,:)


    call object_alloc ('cos_pjasee',  cos_pjasee, nbasis_pw_sim, nelec * (nelec-1)/2)
    call object_alloc ('sin_pjasee',  sin_pjasee, nbasis_pw_sim, nelec * (nelec-1)/2)


    allocate( cos_temp (-maxkv_sim:maxkv_sim,ndim, nelec * (nelec-1)/2))
    allocate( sin_temp (-maxkv_sim:maxkv_sim,ndim, nelec * (nelec-1)/2))


    ije = 0

    do ie = 2, nelec

       do je = 1, ie -1

          ije = ije + 1

          do i=1,ndim

             dot=0
             do k=1,ndim
                dot = dot + glatt(k,i)* xvec(k,ije)
             enddo

             cos_temp(1,i,ije)=cos(dot)
             sin_temp(1,i,ije)=sin(dot)
             cos_temp(-1,i,ije)=cos_temp(1,i,ije)
             sin_temp(-1,i,ije)=-sin_temp(1,i,ije)
             cos_temp(0,i,ije)=1.d0
             sin_temp(0,i,ije)=0.d0

             do n=2, kvec_pjasee (i)

                cos_temp(n,i,ije)=cos_temp(n-1,i,ije)*cos_temp(1,i,ije)- &
                     & sin_temp(n-1,i,ije)*sin_temp(1,i,ije)

                sin_temp(n,i,ije)=sin_temp(n-1,i,ije)*cos_temp(1,i,ije)+ &
                     & cos_temp(n-1,i,ije)*sin_temp(1,i,ije)

                cos_temp(-n,i,ije)=cos_temp(n,i,ije)

                sin_temp(-n,i,ije)=-sin_temp(n,i,ije)

             enddo

          enddo
       enddo

    enddo



    do i = 1 , nbasis_pw_sim

       ije = 0

       do ie = 2, nelec

          do je = 1, ie -1

             ije = ije + 1

             cos_tmp=cos_temp(kv_sim(1,i),1,ije)*cos_temp(kv_sim(2,i),2,ije) &
                  &           -sin_temp(kv_sim(1,i),1,ije)*sin_temp(kv_sim(2,i),2,ije)

             sin_tmp=sin_temp(kv_sim(1,i),1,ije)*cos_temp(kv_sim(2,i),2,ije) &
                  &           +cos_temp(kv_sim(1,i),1,ije)*sin_temp(kv_sim(2,i),2,ije)

             cos_pjasee(i, ije)= cos_tmp*cos_temp(kv_sim(3,i),3,ije) &
                  &               -sin_tmp*sin_temp(kv_sim(3,i),3,ije)

             sin_pjasee(i, ije)= sin_tmp*cos_temp(kv_sim(3,i),3,ije) &
                  &               +cos_tmp*sin_temp(kv_sim(3,i),3,ije)

          enddo
       enddo
    enddo

    deallocate (cos_temp, sin_temp)

  end subroutine calc_cos_sin_ee


end module pjasee_mod

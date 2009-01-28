module pjasen_mod

  use all_tools_mod
  use crystal_symmetries_mod
  use eloc_mod

  implicit none

  real(dp)                               :: shift_cos= 0
  integer                                :: n_inv =  1
  integer                                :: ndim_pj = 3
  ! non-loc pot
  real(dp)                               :: psid_pjas
  real(dp), allocatable                  :: dvpsp_pjas (:)
  real(dp), allocatable                  :: deloc_pot_nloc_pjas(:)
  real(dp), allocatable                  :: dvpot_pjas (:,:)

  integer                                :: param_pjas_nb
  integer                                :: param_pjasen_1_nb
  real(dp), allocatable                  :: dpsi_pjas (:)
  real(dp), allocatable                  :: deloc_pjas (:)
  real(dp), allocatable                  :: pjas_parms_sav (:)
  real(dp), allocatable                  :: pjas_parms (:,:)


  real(dp), allocatable                  :: lap_dpsi_pjasen(:,:)
  real(dp), allocatable                  :: grad_dpsi_pjasen (:,:,:)
  real(dp), allocatable                  :: cos_star_en (:,:)
  real(dp), allocatable                  :: star_en (:,:)
  real(dp), allocatable                  :: grad_cos_star_en (:,:,:)
  real(dp), allocatable                  :: grad_star_en (:,:,:)

  integer                                :: param_pjasen_nb
  real(dp), allocatable                  :: dpsi_pjasen (:)
  real(dp), allocatable                  :: deloc_pjasen (:)

  real(dp), allocatable                  :: pjasen_parms (:,:)

  real(dp), allocatable                  :: gn_pjasen  (:)
!!! new

  real (dp) , allocatable                :: pjasfso(:,:), pjasfijo(:,:,:), pjasd2ijo(:,:), pjasfjo(:,:)
  real (dp) , allocatable                :: pjasdiv_vj (:)
  real(dp)                               :: pjasd2o, pjasfsumo
  real (dp) , allocatable                :: pengo (:, :, :), peego (:, :, :)
  real (dp), allocatable                 :: pjasfijn (:,:,:), pjasfsn (:,:)

  real (dp), allocatable                 :: cos_pjasen(:,:), sin_pjasen(:,:)
  real (dp), allocatable                 :: cos_pjasen_e (:), sin_pjasen_e (:)

contains


  subroutine deloc_pjasen_bld
!---------------------------------------------------------------------------
! Description : derivative of local energy wrt periodic jastrow parameters
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    integer                              :: ist, ie , k
    real(dp)                             :: sum
    include 'commons.h'
!JT    include 'common_vd.h'
    ! header
    if (header_exe) then
       call object_create ('deloc_pjasen',deloc_pjasen_index)
       call object_needed ('nelec')
       call object_needed ('ndim')
       call object_needed ('param_pjasen_nb')
       call object_needed ('vd')
       call object_needed ('vj')
       call object_needed ('grad_dpsi_pjasen')
       call object_needed ('lap_dpsi_pjasen')
       return
    endif

    call object_alloc ('deloc_pjasen', deloc_pjasen, param_pjasen_nb)

    do ist=1,param_pjasen_nb

       sum = 0
       do ie = 1 , nelec
          sum  = sum + lap_dpsi_pjasen (ie, ist)

          do k= 1, ndim_pj
             sum  = sum + 2* grad_dpsi_pjasen (k, ie, ist) * (vd (k,ie)+vj (k,ie))
          enddo

       enddo

       deloc_pjasen (ist) =-hb * sum

    end do

  end subroutine deloc_pjasen_bld



  subroutine  star_en_bld  (xvec)
!---------------------------------------------------------------------------
! Description : building the arrays related to the en stars
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    integer                              :: ie, ist,  istt
    real(dp)                             :: xvec (3,  melec )
    real(dp)                             :: c_s_fac (n_inv), grad_c_s_fac (n_inv,ndim_pj)


    call object_alloc ("grad_star_en ", grad_star_en, ndim_pj, nelec, param_pjasen_nb)

    call object_alloc ("star_en ", star_en , nelec, param_pjasen_nb)

    !! calculate cos and sin of g.xvec
    call calc_cos_sin_en (xvec)

    if (inversion) then

       do ist= 1, param_pjasen_nb

          do ie = 1, nelec

             call star_en_fac (ie, ist, c_s_fac, grad_c_s_fac)

             grad_star_en (:, ie, ist ) = grad_c_s_fac (1,:)

             star_en ( ie, ist ) = c_s_fac (1)

          enddo

       enddo

    else

       istt = 1

       do ist= 1, param_pjasen_nb, 2

          do ie = 1, nelec

             call star_en_fac (ie, istt, c_s_fac, grad_c_s_fac)

             grad_star_en (:, ie, ist ) = grad_c_s_fac (1,:)

             grad_star_en (:, ie, ist +1 ) = grad_c_s_fac (2,:)

             star_en ( ie, ist ) = c_s_fac (1)

             star_en ( ie, ist +1  ) = c_s_fac (2)

          enddo

          istt = istt + 1

       enddo

    endif


    call object_modified_by_index (star_en_index)

    call object_modified_by_index (grad_star_en_index)

  end subroutine star_en_bld



  subroutine deriv_nonloc_pjas_en ( iel, xvec, value1)
!---------------------------------------------------------------------------
! Description : derivatives of the potentials terms with respect to the en Jastrow parameters
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    !    include "common_jasn.h"
    integer                              :: iel, ist , istt
    real(dp)                             :: xvec (3), value1, temp, cos_sin(2), cosine


    call object_alloc ("gn_pjasen", gn_pjasen, param_pjasen_nb)

    do ist = 1 , param_pjasen_nb
       gn_pjasen (ist) = dpsi_pjasen (ist)
    end do

    !! add only updated electron. The old falue is already subtracted
    !! through fso

    value1 = 0

    !!
    call calc_cos_sin_en_e (xvec)


    if (inversion) then

       do ist = 1, param_pjasen_nb

          cosine = cos_star_fac  (ist) !! only changed electron

          temp =  pjas_parms(ist, iwf) * cosine

          value1 = value1 + temp

!!$          fsn(iel,iel) = fsn(iel,iel) + temp

          !! update the gradient.
          !! add new value and subtract old value stored in pengo (iel, iel, ist)
          gn_pjasen (ist) =  gn_pjasen (ist) + cosine  - pengo (iel, iel, ist)

       enddo

    else

       istt = 1

       do ist = 1, param_pjasen_nb, 2

          cos_sin = cos_sin_star_fac  (istt) !! only changed electron

          temp  =  pjas_parms(ist, iwf) * cos_sin (1) +  pjas_parms(ist+1, iwf) * cos_sin (2)

          value1 = value1 +  temp

!!$          fsn(iel,iel) = fsn(iel,iel) + temp

          !! update the gradient.
          !! add new value and subtract old value stored in pengo (iel, iel, ist)
          gn_pjasen (ist) =  gn_pjasen (ist) + cos_sin (1) - pengo (iel, iel, ist)

          gn_pjasen (ist+1) =  gn_pjasen (ist+1) + cos_sin (2) - pengo (iel, iel, ist+1)

          istt = istt + 1

       enddo

    endif

  end subroutine deriv_nonloc_pjas_en




  subroutine nonloc_pjas_en ( iel, xvec, value1)
!---------------------------------------------------------------------------
! Description : compute only the non-local contribution
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    integer                              :: iel, ist , istt
    real(dp)                             :: xvec (3), value1, temp, temp1 (2)

    value1 = 0

    !! add only updated electron. The old falue is already subtracted
    !! through fso

    call calc_cos_sin_en_e (xvec)

    if (inversion) then

       do ist = 1 , param_pjasen_nb

          temp = pjas_parms(ist, iwf) *  cos_star_fac  (ist) !! only changed electron

          value1 = value1 + temp

!!$          fsn(iel,iel) = fsn (iel,iel) + temp

       enddo

    else

       istt= 1

       do ist = 1 , param_pjasen_nb, 2

          temp1 = cos_sin_star_fac  (istt ) !! only changed electron

          temp = pjas_parms(ist, iwf) * temp1 (1) + pjas_parms(ist+1, iwf) * temp1 (2)

          value1 = value1 + temp

!!$          fsn(iel,iel) = fsn (iel,iel) + temp  !!! update fsn (i, i)

          istt = istt +1

       enddo

    endif

  end subroutine nonloc_pjas_en


  function cos_star_fac (ist) result (sum)
!---------------------------------------------------------------------------
! Description : evaluates sum_{star_i,k_star} cos(k_star\dot r)
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    integer                              :: ist
    real(dp)                             :: sum, phas, cosdot
    integer                              :: i, i1, i2

    !! note now basis contains both + and -
    !! sum over the first half only
    i1=istar (ist)
    i2=fstar (ist)
    !    i2=i1+ (mstar (ist))/2 - 1

    sum  = 0

    do i = i1, i2

       phas = phase (i) !! is real for now

       cosdot = cos_pjasen_e (i)

       sum = sum + phas * cosdot !! 2 absorbed in def of parms

    enddo

  end function  cos_star_fac


  function cos_sin_star_fac (ist) result (sum)
!---------------------------------------------------------------------------
! Description :  evaluates the cosine and sine terms
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    integer                              :: ist
    real(dp)                             :: sum (2), phas , cosdot, sindot
    integer                              :: i, i1, i2

    !! note now basis contains both + and -
    !! sum over the first half only
    i1=istar (ist)
    i2=fstar (ist)
    !    i2=i1+ (mstar (ist))/2 - 1

    sum  = 0

    do i = i1, i2
       phas = phase (i) !! is real for now

       cosdot = cos_pjasen_e (i)
       sindot = sin_pjasen_e (i)

       sum (1) = sum (1) + phas * cosdot !! 2 absorbed in def of parms
       sum (2) = sum (2) + phas * sindot !! 2 absorbed in def of parms
    enddo

  end function  cos_sin_star_fac


  function grad_cos_star_fac (ist) result (sum)
!---------------------------------------------------------------------------
! Description : compute the gradient of the cosine star
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    integer                              :: ist
    real(dp)                             :: sum (ndim_pj), phas, sindot
    integer                              :: i, i1, i2, k

    !! note now basis contains both + and -
    !! sum over the first half only
    i1=istar (ist)
    i2=fstar (ist)
    !    i2=i1+ (mstar (ist))/2 - 1

    sum  = 0

    do i = i1, i2

       phas = phase (i) !! is real for now

       sindot = sin_pjasen_e (i)

       do k= 1, ndim_pj
          sum (k) = sum (k) - rkv (k,i)* phas *  sindot
       enddo
    enddo

  end function  grad_cos_star_fac


  function grad_cos_sin_star_fac (ist) result (sum)
!---------------------------------------------------------------------------
! Description : computes the gradient of the cosine and sine terms
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    integer                              :: ist
    real(dp)                             :: sum (2,ndim_pj), phas, sindot, cosdot
    integer                              :: i, i1, i2, k

    !! note now basis contains both + and -
    !! sum over the first half only
    i1=istar (ist)
    i2=fstar (ist)
    !    i2=i1+ (mstar (ist))/2 - 1

    sum  = 0

    do i = i1, i2

       phas = phase (i) !! is real for now

       cosdot = cos_pjasen_e (i)
       sindot = sin_pjasen_e (i)

       do k= 1, ndim_pj
          sum (1,k) = sum (1,k) - rkv (k,i)* phas *  sindot
          sum (2,k) = sum (2,k) + rkv (k,i)* phas *  cosdot
       enddo
    enddo

  end function  grad_cos_sin_star_fac




  subroutine star_en_fac ( ie, ist, c_s_fac, grad_c_s_fac)
!---------------------------------------------------------------------------
! Description : bulid the stars for the en periodic jastrow
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"

    real(dp)                             :: c_s_fac (n_inv), grad_c_s_fac (n_inv,ndim_pj)
    integer                              :: ist, ie
    real(dp)                             :: phas, cost (2)
    integer                              :: i, i1, i2, k

    !! note now basis contains both + and -
    !! sum over the first half only
    i1=istar (ist)
    i2=fstar (ist)
    !    i2=i1+ (mstar (ist))/2 - 1

    grad_c_s_fac  = 0
    c_s_fac = 0

    if (inversion) then

       do i = i1, i2

          phas = phase (i) !! is real for now


          cost = (/ cos_pjasen (i,ie),  sin_pjasen (i,ie) /)

          c_s_fac (1)  =  c_s_fac (1) +  phas * cost (1) !! 2 absorbed in def of parms

          do k= 1, ndim_pj
             grad_c_s_fac (1,k) = grad_c_s_fac (1,k) - rkv (k,i)* phas *  cost (2)
          enddo
       enddo

    else

       do i = i1, i2

          phas = phase (i) !! is real for now

          cost = (/ cos_pjasen (i,ie),  sin_pjasen (i,ie) /)


          c_s_fac (1)  =  c_s_fac (1) +  phas * cost (1) !! 2 absorbed in def of parms

          c_s_fac (2)  =  c_s_fac (2) +  phas * cost (2) !! 2 absorbed in def of parms

          do k= 1, ndim_pj
             grad_c_s_fac (1,k) = grad_c_s_fac (1,k) - rkv (k,i)* phas *  cost (2)

             grad_c_s_fac (2,k) = grad_c_s_fac (2,k) + rkv (k,i)* phas *  cost (1)
          enddo

       enddo

    endif

  end subroutine star_en_fac



  subroutine star_en_fac_e ( ie, ist, c_s_fac, grad_c_s_fac)
!---------------------------------------------------------------------------
! Description : Build the stars terms where only one electron is displaced
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"

    real(dp)                             :: c_s_fac (n_inv), grad_c_s_fac (n_inv,ndim_pj)
    integer                              :: ist, ie
    real(dp)                             :: phas, cost (2)
    integer                              :: i, i1, i2, k

    !! note now basis contains both + and -
    !! sum over the first half only
    i1=istar (ist)
    i2=fstar (ist)
    !    i2=i1+ (mstar (ist))/2 - 1

    grad_c_s_fac  = 0
    c_s_fac = 0

    if (inversion) then

       do i = i1, i2

          phas = phase (i) !! is real for now


          cost = (/ cos_pjasen_e (i),  sin_pjasen_e (i) /)

          !          write(88,'(4F15.8)') cost (1) , cos_pjasen (i,ie), cost (2) , sin_pjasen (i,ie)

          c_s_fac (1)  =  c_s_fac (1) +  phas * cost (1) !! 2 absorbed in def of parms

          do k= 1, ndim_pj
             grad_c_s_fac (1,k) = grad_c_s_fac (1,k) - rkv (k,i)* phas *  cost (2)
          enddo
       enddo

    else

       do i = i1, i2

          phas = phase (i) !! is real for now

          cost = (/ cos_pjasen_e (i),  sin_pjasen_e (i) /)

          c_s_fac (1)  =  c_s_fac (1) +  phas * cost (1) !! 2 absorbed in def of parms

          c_s_fac (2)  =  c_s_fac (2) +  phas * cost (2) !! 2 absorbed in def of parms

          do k= 1, ndim_pj
             grad_c_s_fac (1,k) = grad_c_s_fac (1,k) - rkv (k,i)* phas *  cost (2)

             grad_c_s_fac (2,k) = grad_c_s_fac (2,k) + rkv (k,i)* phas *  cost (1)
          enddo

       enddo

    endif

  end subroutine star_en_fac_e



  subroutine  pjasen_deriv_jas (xvec, pjasv, pjasd2, pjasdiv_vj,  pjasfsum)
!---------------------------------------------------------------------------
! Description : computes the gradient of the wavefunction with respect to the periodic jastrow
!               parameters, and  the graident and laplacian of this quantity which is used in the
!               computation of the derivates of the local energy with respect to these jastrow parameters
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    real (dp)                            :: xvec (3, melec)
    real (dp)                            :: pjasv (3,melec)
    real (dp)                            :: pjasfsum
    real (dp)                            :: pjasdiv_vj (melec), pjasd2
    real (dp)                            :: gen,  gradt (3)
    integer                              :: ist, i

    call star_en_bld (xvec)

    do ist =1, param_pjasen_nb
       dpsi_pjasen(ist)=0
       do i =1, nelec
          lap_dpsi_pjasen (i,ist)=0
          grad_dpsi_pjasen(:,i,ist)=0
       enddo
    enddo

    ! e-n terms
    do i = 1, nelec

       do ist = 1, param_pjasen_nb

          gen = star_en (i, ist)

          gradt = grad_star_en (:, i, ist)

          pjasfso (i,i) = pjasfso(i, i) +  pjas_parms (ist, iwf) * gen

          pjasfijo(:, i,i) = pjasfijo (:,i, i) +  &
               &    pjas_parms (ist, iwf) * gradt

          pjasd2ijo (i,i) = pjasd2ijo (i,i)  -  &
               &    pjas_parms (ist, iwf) *  sk3 (ist) * gen

          pengo(i,i,ist) = pengo(i,i,ist)+ gen

          dpsi_pjasen (ist) = dpsi_pjasen (ist) + gen

          grad_dpsi_pjasen (:,i,ist)  = grad_dpsi_pjasen (:, i, ist) + gradt

          lap_dpsi_pjasen (i, ist) =  lap_dpsi_pjasen (i, ist) - &
               & sk3 (ist) * gen

       enddo

       pjasfsum = pjasfsum + pjasfso(i,i)

       pjasfjo(:,i) = pjasfjo(:,i) + pjasfijo(:,i,i)

       pjasdiv_vj (i)= pjasdiv_vj (i) + pjasd2ijo(i,i)

       pjasd2 = pjasd2 + pjasd2ijo(i,i)

    end do

    do i=1,nelec
       pjasv (:, i) = pjasv (:, i) + pjasfjo (:,i)
    enddo

    call object_modified_by_index (dpsi_pjasen_index)
    call object_modified_by_index (grad_dpsi_pjasen_index)
    call object_modified_by_index (lap_dpsi_pjasen_index)

  end subroutine pjasen_deriv_jas



  subroutine  pjasen_jas (xvec, pjasv, pjasd2, pjasdiv_vj,  pjasfsum)
!---------------------------------------------------------------------------
! Description : same as pjasen_deriv_jas  but without the derivatives of the
!               \nabla J with respect periodic Jastrow parameters.
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    real (dp)                            :: xvec (3, melec)
    real (dp)                            :: pjasv (3,melec)
    real (dp)                            :: pjasfsum
    real (dp)                            :: pjasdiv_vj (melec), pjasd2
    real (dp)                            :: fen
    integer                              :: ist, i

    do ist = 1, param_pjasen_nb
       dpsi_pjasen(ist)=0
    enddo

    call star_en_bld (xvec)


    do i = 1, nelec

       do ist = 1, param_pjasen_nb

          fen = star_en (i, ist)

          pjasfso (i,i) = pjasfso(i, i) +  pjas_parms (ist, iwf)* fen

          pjasfijo(:, i,i) = pjasfijo (:,i, i) + &
               & pjas_parms (ist, iwf) * grad_star_en (:, i, ist)

          pjasd2ijo (i,i)  = pjasd2ijo (i,i)  - &
               & pjas_parms (ist, iwf) *  sk3 (ist) * fen

          dpsi_pjasen (ist) = dpsi_pjasen (ist) + fen

       enddo

       pjasfsum = pjasfsum + pjasfso(i,i)

       pjasfjo(:,i) = pjasfjo(:,i) + pjasfijo(:,i,i)

       pjasdiv_vj (i)= pjasdiv_vj (i) + pjasd2ijo(i,i)

       pjasd2 = pjasd2 + pjasd2ijo(i,i)

    end do


    do i=1,nelec
       pjasv (:, i) = pjasv (:, i) + pjasfjo (:,i)
    enddo

    !! needed in nonloc
    call object_modified_by_index (dpsi_pjasen_index)

  end subroutine pjasen_jas



  subroutine  pjasen_jas_e (iel, xvec, fsum)
!---------------------------------------------------------------------------
! Description : Same as  pjasen_jas  but when only one electron is modified
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
!JT    include 'common_jasn.h'
    integer                              :: iel
    real (dp)                            :: xvec (3, melec)
    real (dp)                            :: fen , fen1 (3), fsum
    integer                              :: ist, istt
    real(dp)                             :: c_s_fac (n_inv), grad_c_s_fac (n_inv,ndim_pj)

    fen = 0
    fen1 = 0

    call calc_cos_sin_en_e (xvec(:,iel))

    if (inversion) then

       do ist = 1, param_pjasen_nb

          call star_en_fac_e (iel, ist, c_s_fac, grad_c_s_fac)

          fen = fen +  pjas_parms (ist, iwf)* c_s_fac (1)

          fen1= fen1 + pjas_parms (ist, iwf)* grad_c_s_fac (1,:)

       enddo

    else

       istt = 1

       do ist = 1, param_pjasen_nb, 2

          call star_en_fac_e (iel, istt, c_s_fac, grad_c_s_fac)

          fen = fen +  pjas_parms (ist, iwf)* c_s_fac (1) + &
               & pjas_parms (ist+1, iwf)* c_s_fac (2)

          fen1= fen1 + pjas_parms (ist, iwf)* grad_c_s_fac (1,:) + &
               & pjas_parms (ist+1, iwf)* grad_c_s_fac (2,:)

          istt = istt + 1

       enddo
    endif

    fsn(iel,iel) = fsn (iel,iel) + fen

    fijn(:,iel, iel) = fijn (:,iel, iel) + fen1

!!! add and subtract old
!!! old is already subtracted with short range jastrow jastrow
!!! warning: problem if the short range jastrow is removed
    fjn(:, iel ) = fjn (:, iel) + fen1

    fsum = fsum + fen

  end subroutine pjasen_jas_e




  subroutine  calc_cos_sin_en_e (xvec)
!---------------------------------------------------------------------------
! Description : calculate the cosine and sine of the star terms when one electron is displaced.
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    integer                              :: i, k, n
    real(dp)                             :: xvec (3)
    real(dp)                             :: dot, cos_tmp, sin_tmp
    real (dp), allocatable               :: cos_temp(:,:),  sin_temp(:,:)


    call object_alloc ('cos_pjasen_e',  cos_pjasen_e, nbasis_pw)
    call object_alloc ('sin_pjasen_e',  sin_pjasen_e, nbasis_pw)


    allocate( cos_temp (-maxkv:maxkv,ndim ))
    allocate( sin_temp (-maxkv:maxkv,ndim))

    do i=1,ndim
       dot=0
       do k=1,ndim
          dot = dot + glatt(k,i)* xvec(k)
       enddo
       cos_temp(1,i)=cos(dot)
       sin_temp(1,i)=sin(dot)
       cos_temp(-1,i)=cos_temp(1,i)
       sin_temp(-1,i)=-sin_temp(1,i)
       cos_temp(0,i)=1.d0
       sin_temp(0,i)=0.d0

       do n=2, kvec_pjasen (i)

          cos_temp(n,i)=cos_temp(n-1,i)*cos_temp(1,i)- &
               & sin_temp(n-1,i)*sin_temp(1,i)

          sin_temp(n,i)=sin_temp(n-1,i)*cos_temp(1,i)+ &
               & cos_temp(n-1,i)*sin_temp(1,i)

          cos_temp(-n,i)=cos_temp(n,i)

          sin_temp(-n,i)=-sin_temp(n,i)

       enddo

    enddo

    do i = 1 , nbasis_pw

       cos_tmp=cos_temp(kv(1,i),1)*cos_temp(kv(2,i),2) &
            &           -sin_temp(kv(1,i),1)*sin_temp(kv(2,i),2)

       sin_tmp=sin_temp(kv(1,i),1)*cos_temp(kv(2,i),2) &
            &           +cos_temp(kv(1,i),1)*sin_temp(kv(2,i),2)

       cos_pjasen_e (i)= cos_tmp*cos_temp(kv(3,i),3) &
            &               -sin_tmp*sin_temp(kv(3,i),3)

       sin_pjasen_e (i)= sin_tmp*cos_temp(kv(3,i),3) &
            &               +cos_tmp*sin_temp(kv(3,i),3)
    enddo

    deallocate (cos_temp, sin_temp)

  end subroutine calc_cos_sin_en_e



  subroutine  calc_cos_sin_en (xvec)
!---------------------------------------------------------------------------
! Description : calculate the cosine and sine of the star terms.
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    integer                              :: ie, i, k, n
    real(dp)                             :: xvec (3,  melec )
    real(dp)                             :: dot, cos_tmp, sin_tmp
    real (dp), allocatable               :: cos_temp(:,:,:),  sin_temp(:,:,:)


    call object_alloc ('cos_pjasen',  cos_pjasen, nbasis_pw, nelec)
    call object_alloc ('sin_pjasen',  sin_pjasen, nbasis_pw, nelec)


    allocate( cos_temp (-maxkv:maxkv,ndim, nelec))
    allocate( sin_temp (-maxkv:maxkv,ndim, nelec))


    do ie=1,nelec

       do i=1,ndim

          dot=0
          do k=1,ndim
             dot = dot + glatt(k,i)* xvec(k,ie)
          enddo

          cos_temp(1,i,ie)=cos(dot)
          sin_temp(1,i,ie)=sin(dot)
          cos_temp(-1,i,ie)=cos_temp(1,i,ie)
          sin_temp(-1,i,ie)=-sin_temp(1,i,ie)
          cos_temp(0,i,ie)=1.d0
          sin_temp(0,i,ie)=0.d0

          do n=2, kvec_pjasen (i)

             cos_temp(n,i,ie)=cos_temp(n-1,i,ie)*cos_temp(1,i,ie)- &
                  & sin_temp(n-1,i,ie)*sin_temp(1,i,ie)

             sin_temp(n,i,ie)=sin_temp(n-1,i,ie)*cos_temp(1,i,ie)+ &
                  & cos_temp(n-1,i,ie)*sin_temp(1,i,ie)

             cos_temp(-n,i,ie)=cos_temp(n,i,ie)

             sin_temp(-n,i,ie)=-sin_temp(n,i,ie)

          enddo

       enddo
    enddo



    do i = 1 , nbasis_pw

       do  ie =1 , nelec

          cos_tmp=cos_temp(kv(1,i),1,ie)*cos_temp(kv(2,i),2,ie) &
               &           -sin_temp(kv(1,i),1,ie)*sin_temp(kv(2,i),2,ie)

          sin_tmp=sin_temp(kv(1,i),1,ie)*cos_temp(kv(2,i),2,ie) &
               &           +cos_temp(kv(1,i),1,ie)*sin_temp(kv(2,i),2,ie)

          cos_pjasen(i, ie)= cos_tmp*cos_temp(kv(3,i),3,ie) &
               &               -sin_tmp*sin_temp(kv(3,i),3,ie)

          sin_pjasen(i, ie)= sin_tmp*cos_temp(kv(3,i),3,ie) &
               &               +cos_tmp*sin_temp(kv(3,i),3,ie)

       enddo
    enddo

    deallocate (cos_temp, sin_temp)

  end subroutine calc_cos_sin_en

!!$
!!$
!!$
!!$   subroutine plot_jasen
!---------------------------------------------------------------------------
! Description :
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
!!$     implicit none
!!$     include 'commons.h'
!!$     integer                          :: i,k,be,j,npoints
!!$     real(dp)                         :: r (3) ,dr,maxdis
!!$
!!$     logical,save                     :: first=.true.
!!$     logical,save                     :: first1=.true.
!!$
!!$!!! only main processor will enter
!!$    if (idtask .ne. 0) return
!!$
!!$    write(*,*) "calling plot_phi"
!!$
!!$    if (first) then
!!$       open(181,file="plot_jasen",status="replace")
!!$       first=.false.
!!$    endif
!!$
!!$
!!$    r=0.01
!!$    dr=0.05
!!$
!!$    maxdis=cutjas_en
!!$
!!$    npoints=maxdis/dr+10
!!$
!!$
!!$    if (inversion) then
!!$
!!$       do i=1,npoints
!!$
!!$          r(1)= i * dr
!!$
!!$          r(2)= i * dr
!!$
!!$          r(3)= i * dr
!!$
!!$          sum1 = 0
!!$
!!$          do ist = 1,  param_pjasen_nb
!!$
!!$             cos  = cos_star_fac  (ist, xvec (1:ndim_pj,ie))
!!$
!!$             sum1 = sum1 + pjas_parms (ist, iwf) * cos_sin (1)
!!$
!!$             if (inversion) sum2 = sum2 + cos_sin (2)
!!$
!!$          enddo
!!$
!!$       write(181,*) r, sum1, sum2
!!$
!!$
!!$    enddo
!!$
!!$    else
!!$
!!$
!!$
!!$       do i=1,npoints
!!$
!!$          r(1)= i * dr
!!$
!!$          r(2)= i * dr
!!$
!!$          r(3)= i * dr
!!$
!!$          sum1 = 0
!!$          do ist = 1, nstar_en
!!$
!!$             cos_sin = cos_sin_star_fac (ist, r)
!!$
!!$             sum1 = sum1 + pjas_parms (ist, iwf) * cos_sin (1)
!!$
!!$             if (inversion) sum2 = sum2 + cos_sin (2)
!!$
!!$          enddo
!!$
!!$       write(181,*) r, sum1, sum2
!!$
!!$
!!$    enddo
!!$
!!$
!!$ endif
!!$
!!$
!!$    write(181,*) "&"
!!$
!!$  end subroutine plot_jasen
!!$


end module pjasen_mod



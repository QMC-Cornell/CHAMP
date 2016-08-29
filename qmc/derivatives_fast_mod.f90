module derivatives_fast_mod

  use all_tools_mod
  use electrons_mod
  use psi_mod

! Declaration of global variables and default values
! ex_info_bld 
  real(dp)                       :: coef_ref
  integer                        :: det_ref_up
  integer                        :: det_ref_dn
  integer                        :: orb_occ_up_nb
  integer                        :: orb_occ_dn_nb
  integer                        :: orb_vir_up_nb
  integer                        :: orb_vir_dn_nb
  integer, allocatable           :: orb_occ_lab_up(:)
  integer, allocatable           :: orb_occ_lab_dn(:)
  integer, allocatable           :: orb_vir_lab_up(:)
  integer, allocatable           :: orb_vir_lab_dn(:)
  integer, allocatable           :: orb_occ_lab_rev_up(:)
  integer, allocatable           :: orb_occ_lab_rev_dn(:)
  integer, allocatable           :: orb_vir_lab_rev_up(:)
  integer, allocatable           :: orb_vir_lab_rev_dn(:)
  integer                        :: ex_order_up_max
  integer                        :: ex_order_dn_max
  integer, allocatable           :: ex_from_up(:,:)
  integer, allocatable           :: ex_from_dn(:,:)
  integer, allocatable           :: ex_to_up(:,:)
  integer, allocatable           :: ex_to_dn(:,:)
  integer, allocatable           :: ex_order_up(:)
  integer, allocatable           :: ex_order_dn(:)

! ainv_atilde_bld 
  real(dp), allocatable          :: ainv_atilde_up(:,:)
  real(dp), allocatable          :: ainv_atilde_dn(:,:)

! alphaI_det_inv_bld
  real(dp), allocatable          :: alphaI_inv_up(:,:,:)
  real(dp), allocatable          :: alphaI_inv_dn(:,:,:)
  real(dp), allocatable          :: alphaI_det_up(:)
  real(dp), allocatable          :: alphaI_det_dn(:)

! ymat_ainv_and_yKmat_bld
  real(dp), allocatable          :: ymat_ainv_up(:,:)
  real(dp), allocatable          :: ymat_ainv_dn(:,:)
  real(dp), allocatable          :: yKmat_up(:,:,:)
  real(dp), allocatable          :: yKmat_dn(:,:,:)

! myphi_from_yKmat_bld 
  real(dp)                       :: myphi

! gamma_bld
  real(dp), allocatable          :: gamma_up(:,:)
  real(dp), allocatable          :: gamma_dn(:,:)

! grd_det_over_det_fast_bld
  real(dp), allocatable          :: grd_det_over_det_fast(:,:,:)
  real(dp), allocatable          :: grd_det_over_det_legacy(:,:)

! lap_det_over_det_fast_bld
  real(dp), allocatable          :: lap_det_over_det_fast(:)
  real(dp), allocatable          :: lap_det_over_det_fast_av(:)
  real(dp), allocatable          :: lap_det_over_det_fast_av_err(:)
  real(dp)                       :: lap_det_over_det_legacy

  contains

! ==============================================================================
  subroutine ex_info_bld
! ------------------------------------------------------------------------------
! Description   : In an expansion 
! Description   :  Phi = \sum_I C_I |CSF>
! Description   :      = \sum_I \sum_kI C_I c_kI |Dup>|Ddn>
! Description   :  -  find the main determinant if not given by the user
! Description   :  -  compute maps for the occupied and virtual orbitals
! Description   :     wrt the main determinant
! Description   :  -  compute the excitation information for all determinant
! Description   :     ("from" and "to" orbitals as well as order of the excitation
! Description   :      wrt the main determinant)
!
! Created       : B. Mussard, Aug. 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'ex_info_bld'
  integer :: csf_i, det_i, orb_i, det_in_csf_i, ex_i, from, to
  integer :: det_up, det_dn
  integer :: det_unq_up, det_unq_dn
  integer :: ex_max

! header
  if (header_exe) then

   call object_create ('coef_ref')
   call object_create ('det_ref_up')
   call object_create ('det_ref_dn')

   call object_create ('orb_occ_up_nb')
   call object_create ('orb_occ_dn_nb')
   call object_create ('orb_vir_up_nb')
   call object_create ('orb_vir_dn_nb')
   call object_create ('orb_occ_lab_up')
   call object_create ('orb_occ_lab_dn')
   call object_create ('orb_vir_lab_up')
   call object_create ('orb_vir_lab_dn')
   call object_create ('orb_occ_lab_rev_up')
   call object_create ('orb_occ_lab_rev_dn')
   call object_create ('orb_vir_lab_rev_up')
   call object_create ('orb_vir_lab_rev_dn')

   call object_create ('ex_order_up_max')
   call object_create ('ex_order_dn_max')
   call object_create ('ex_from_up')
   call object_create ('ex_from_dn')
   call object_create ('ex_to_up')
   call object_create ('ex_to_dn')
   call object_create ('ex_order_up')
   call object_create ('ex_order_dn')

   call object_needed('ncsf')
   call object_needed('ndet_in_csf')
   call object_needed('iwdet_in_csf')
   call object_needed('det_to_det_unq_up')
   call object_needed('det_to_det_unq_dn')
   call object_needed('csf_coef')
   call object_needed('cdet_in_csf')
   call object_needed('orb_tot_nb')
   call object_needed('nup')
   call object_needed('ndn')
   call object_needed('orb_occ_in_det_unq_up')
   call object_needed('orb_occ_in_det_unq_dn')
   call object_needed('ndetup')
   call object_needed('ndetdn')

   return

  endif

! find main determinant as the largest C_I x c_kI
! (or it could have been imposed by user)
  if (det_ref_up.eq.0) then
    coef_ref=0.d0
    do csf_i=1, ncsf
      do det_in_csf_i=1, ndet_in_csf(csf_i)
        det_i=iwdet_in_csf(det_in_csf_i, csf_i)
        det_up=det_to_det_unq_up(det_i)
        det_dn=det_to_det_unq_dn(det_i)
        if (csf_coef(csf_i, 1)*cdet_in_csf(det_in_csf_i, csf_i).gt.coef_ref) then
          coef_ref=csf_coef(csf_i, 1)*cdet_in_csf(det_in_csf_i, csf_i)
          det_ref_up=det_up
          det_ref_dn=det_dn
        endif
      enddo
    enddo
  endif
  write(6,*) 'main_det', det_ref_up, det_ref_dn, coef_ref

! find occupied and virtual orbitals
! wrt the main determinant
  orb_occ_up_nb=0
  orb_vir_up_nb=0
  orb_occ_dn_nb=0
  orb_vir_dn_nb=0
  call object_alloc ('orb_occ_lab_up',     orb_occ_lab_up,     nup)
  call object_alloc ('orb_vir_lab_up',     orb_vir_lab_up,     orb_tot_nb-nup)
  call object_alloc ('orb_occ_lab_dn',     orb_occ_lab_dn,     ndn)
  call object_alloc ('orb_vir_lab_dn',     orb_vir_lab_dn,     orb_tot_nb-ndn)
  call object_alloc ('orb_occ_lab_rev_up', orb_occ_lab_rev_up, orb_tot_nb)
  call object_alloc ('orb_vir_lab_rev_up', orb_vir_lab_rev_up, orb_tot_nb)
  call object_alloc ('orb_occ_lab_rev_dn', orb_occ_lab_rev_dn, orb_tot_nb)
  call object_alloc ('orb_vir_lab_rev_dn', orb_vir_lab_rev_dn, orb_tot_nb)

  do orb_i=1, orb_tot_nb
!   up
    if (orb_occ_in_det_unq_up(orb_i, det_ref_up)) then
      orb_occ_up_nb=orb_occ_up_nb+1
      orb_occ_lab_up(orb_occ_up_nb)=orb_i
      orb_occ_lab_rev_up(orb_i)=orb_occ_up_nb
    else
      orb_vir_up_nb=orb_vir_up_nb+1
      orb_vir_lab_up(orb_vir_up_nb)=orb_i
      orb_vir_lab_rev_up(orb_i)=orb_vir_up_nb
    endif

!   dn
    if (orb_occ_in_det_unq_dn(orb_i, det_ref_dn)) then
      orb_occ_dn_nb=orb_occ_dn_nb+1
      orb_occ_lab_dn(orb_occ_dn_nb)=orb_i
      orb_occ_lab_rev_dn(orb_i)=orb_occ_dn_nb
    else
      orb_vir_dn_nb=orb_vir_dn_nb+1
      orb_vir_lab_dn(orb_vir_dn_nb)=orb_i
      orb_vir_lab_rev_dn(orb_i)=orb_vir_dn_nb
    endif
  enddo
  write(6,*) 'orb_occ_up',     orb_occ_lab_up
  write(6,*) 'orb_vir_up',     orb_vir_lab_up
  write(6,*) 'orb_occ_dn',     orb_occ_lab_dn
  write(6,*) 'orb_vir_dn',     orb_vir_lab_dn
  write(6,*) 'orb_occ_rev_up', orb_occ_lab_rev_up
  write(6,*) 'orb_vir_rev_up', orb_vir_lab_rev_up
  write(6,*) 'orb_occ_rev_dn', orb_occ_lab_rev_dn
  write(6,*) 'orb_vir_rev_dn', orb_vir_lab_rev_dn

! gather excitation information:
!  -  "from" and "to" orbitals
!  -  order of the excitation
  ex_max=0
  call object_alloc ('ex_order_up', ex_order_up, ndetup)
  call object_alloc ('ex_order_dn', ex_order_dn, ndetdn)
! these arrays will be adjusted for actual size at the end
  ex_order_up_max=5
  ex_order_dn_max=5
  call object_alloc ('ex_from_up',  ex_from_up,  ndetup, ex_order_up_max)
  call object_alloc ('ex_from_dn',  ex_from_dn,  ndetdn, ex_order_dn_max)
  call object_alloc ('ex_to_up',    ex_to_up,    ndetup, ex_order_up_max)
  call object_alloc ('ex_to_dn',    ex_to_dn,    ndetdn, ex_order_dn_max)

! up
  do det_unq_up=1, ndetup
    ex_i=0
    do orb_i=1, orb_tot_nb
      if     ((     orb_occ_in_det_unq_up(orb_i, det_ref_up))&
         .and.(.not.orb_occ_in_det_unq_up(orb_i, det_unq_up))) then
        from=orb_i
      elseif ((.not.orb_occ_in_det_unq_up(orb_i, det_ref_up))&
         .and.(     orb_occ_in_det_unq_up(orb_i, det_unq_up))) then
        to=orb_i
      endif
      if ((from.ne.0).and.(to.ne.0)) then
        ex_i=ex_i+1
        if (ex_i>ex_order_up_max) call die(lhere, "Recompile with greater value of 'ex_order_up_max'")
        if (ex_i>ex_max) ex_max=ex_i
        ex_from_up(det_unq_up, ex_i)=from
        ex_to_up  (det_unq_up, ex_i)=to
        from=0
        to=0
      endif
    enddo
    ex_order_up(det_unq_up)=ex_i
  enddo
  ex_order_up_max=ex_max
  write(6,*) 'ex_order_up_max', ex_order_up_max

! adjust size of the arrays
  call object_alloc ('ex_from_dn', ex_from_dn, ndetdn, ex_order_dn_max)
  call object_alloc ('ex_to_dn',   ex_to_dn,   ndetdn, ex_order_dn_max)
  do det_unq_up=1, ndetup
  write(6,*) 'det_unq_up', det_unq_up, ex_order_up(det_unq_up), ex_from_up(det_unq_up,:), ex_to_up(det_unq_up,:)
  enddo

! dn
  do det_unq_dn=1, ndetdn
    ex_i=0
    do orb_i=1, orb_tot_nb
      if     ((     orb_occ_in_det_unq_dn(orb_i, det_ref_dn))&
         .and.(.not.orb_occ_in_det_unq_dn(orb_i, det_unq_dn))) then
        from=orb_i
      elseif ((.not.orb_occ_in_det_unq_dn(orb_i, det_ref_dn))&
         .and.(     orb_occ_in_det_unq_dn(orb_i, det_unq_dn))) then
        to=orb_i
      endif
      if ((from.ne.0).and.(to.ne.0)) then
        ex_i=ex_i+1
        if (ex_i>ex_order_dn_max) call die(lhere, "Recompile with greater value of 'ex_order_dn_max'")
        if (ex_i>ex_max) ex_max=ex_i
        ex_from_dn(det_unq_dn, ex_i)=from
        ex_to_dn  (det_unq_dn, ex_i)=to
        from=0
        to=0
      endif
    enddo
    ex_order_dn(det_unq_dn)=ex_i
  enddo
  ex_order_dn_max=ex_max
  write(6,*) 'ex_order_dn_max', ex_order_dn_max

! adjust size of the arrays
  call object_alloc ('ex_from_up', ex_from_up, ndetup, ex_order_up_max)
  call object_alloc ('ex_to_up',   ex_to_up,   ndetup, ex_order_up_max)
  do det_unq_dn=1, ndetup
  write(6,*) 'det_unq_dn', det_unq_dn, ex_order_dn(det_unq_dn), ex_from_dn(det_unq_dn,:), ex_to_dn(det_unq_dn,:)
  enddo

  end subroutine ex_info_bld

! ==============================================================================
  subroutine ainv_atilde_bld
! ------------------------------------------------------------------------------
! Description   : Compute A^-1.Atilde (UP and DN)
!
! Created       : B. Mussard, Aug. 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer :: orb_i, orb_j, elec_i, elec_k, link_k, orb_lab

! header
  if (header_exe) then

   call object_create ('ainv_atilde_up')
   call object_create ('ainv_atilde_dn')

   call object_needed('orb_occ_up_nb')
   call object_needed('orb_occ_dn_nb')
   call object_needed('orb_vir_up_nb')
   call object_needed('orb_vir_dn_nb')
   call object_needed('orb_vir_lab_up')
   call object_needed('orb_vir_lab_dn')
   call object_needed('orb')
   call object_needed('slater_mat_trans_inv_up')
   call object_needed('slater_mat_trans_inv_dn')

   return

  endif

  write(6,*) 'orb:'
  do orb_i=1, orb_tot_nb
  write(6,*) (orb(elec_i, orb_i), elec_i=1, nelec)
  enddo
  write(6,*) 'ainv_up:'
  do elec_i=1, orb_occ_up_nb
  write(6,*) (slater_mat_trans_inv_up(orb_j, elec_i, det_ref_up), orb_j=1, orb_occ_up_nb)
  enddo
  write(6,*) 'ainv_dn:'
  do elec_i=1, orb_occ_dn_nb
  write(6,*) (slater_mat_trans_inv_dn(orb_j, elec_i, det_ref_dn), orb_j=1, orb_occ_dn_nb)
  enddo

! atilde_up
  call object_alloc ('ainv_atilde_up', ainv_atilde_up, orb_occ_up_nb, orb_vir_up_nb)
  ainv_atilde_up=0.d0
  do elec_i=1, orb_occ_up_nb
    do orb_j=1, orb_vir_up_nb
      orb_lab=orb_vir_lab_up(orb_j)
      ainv_atilde_up(elec_i, orb_j)=orb(elec_i, orb_lab)
    enddo
  enddo
  write(6,*) 'atilde_up:'
  do elec_i=1, orb_occ_up_nb
  write(6,*) (ainv_atilde_up(elec_i, orb_j), orb_j=1, orb_vir_up_nb)
  enddo

! ainv_atilde_up
  do elec_i=1, orb_occ_up_nb
    do orb_j=1, orb_vir_up_nb
      do link_k=1, orb_occ_up_nb
        ainv_atilde_up(elec_i, orb_j)=slater_mat_trans_inv_up(link_k, elec_i, det_ref_up)*ainv_atilde_up(link_k, orb_j)
      enddo
    enddo
  enddo
  write(6,*) 'ainv_atilde_up:'
  do elec_i=1, orb_occ_up_nb
  write(6,*) (ainv_atilde_up(elec_i, orb_j), orb_j=1, orb_vir_up_nb)
  enddo

! atilde_dn
  call object_alloc ('ainv_atilde_dn', ainv_atilde_dn, orb_occ_dn_nb, orb_vir_dn_nb)
  ainv_atilde_dn=0.d0
  do elec_i=1, orb_occ_dn_nb
    do orb_j=1, orb_vir_dn_nb
      orb_lab=orb_vir_lab_dn(orb_j)
      ainv_atilde_dn(elec_i, orb_j)=orb(elec_i+nup, orb_lab)
    enddo
  enddo
  write(6,*) 'atilde_dn:'
  do elec_i=1, orb_occ_dn_nb
  write(6,*) (ainv_atilde_dn(elec_i, orb_j), orb_j=1, orb_vir_dn_nb)
  enddo

! ainv_atilde_dn
  do elec_i=1, orb_occ_dn_nb
    do orb_j=1, orb_vir_dn_nb
      do link_k=1, orb_occ_dn_nb
        ainv_atilde_dn(elec_i, orb_j)=slater_mat_trans_inv_dn(link_k, elec_i, det_ref_dn)*ainv_atilde_dn(link_k, orb_j)
      enddo
    enddo
  enddo
  write(6,*) 'ainv_atilde_dn:'
  do elec_i=1, orb_occ_dn_nb
  write(6,*) (ainv_atilde_dn(elec_i, orb_j), orb_j=1, orb_vir_dn_nb)
  enddo

  end subroutine ainv_atilde_bld

! ==============================================================================
  subroutine alphaI_det_inv_bld
! ------------------------------------------------------------------------------
! Description   : Compute all det(alpha_I) and alpha_I^-1 (UP and DN)
! Description   : (alpha_I are portions of A^-1.Atilde corresponding to a particular excitation)
!
! Created       : B. Mussard, Aug. 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer :: det_unq_up, det_unq_dn
  integer :: ex_i, ex_j, n, info, elec_i, orb_j
  !real(8), allocatable :: work(:)
  !integer, allocatable :: ipiv(:)

! header
  if (header_exe) then

   call object_create ('alphaI_inv_up')
   call object_create ('alphaI_inv_dn')
   call object_create ('alphaI_det_up')
   call object_create ('alphaI_det_dn')

   call object_needed('ndetup')
   call object_needed('ndetdn')
   call object_needed('ex_order_up_max')
   call object_needed('ex_order_dn_max')
   call object_needed('ex_order_up')
   call object_needed('ex_order_dn')
   call object_needed('orb_occ_lab_rev_up')
   call object_needed('orb_occ_lab_rev_dn')
   call object_needed('orb_vir_lab_rev_up')
   call object_needed('orb_vir_lab_rev_dn')
   call object_needed('ex_from_up')
   call object_needed('ex_from_dn')
   call object_needed('ex_to_up')
   call object_needed('ex_to_dn')
   call object_needed('ainv_atilde_up')
   call object_needed('ainv_atilde_dn')

   return

  endif

! up
  call object_alloc('alphaI_inv_up', alphaI_inv_up, ndetup, ex_order_up_max, ex_order_up_max)
  call object_alloc('alphaI_det_up', alphaI_det_up, ndetup)
  do det_unq_up=1, ndetup
    n=ex_order_up(det_unq_up)
    if (n.eq.0) then
      alphaI_inv_up(det_unq_up, 1, 1)=1.d0
      alphaI_det_up(det_unq_up)=1.d0
    else
      do ex_i=1, n
        do ex_j=1, n
          elec_i=orb_occ_lab_rev_up(ex_from_up(det_unq_up, ex_i))
          orb_j=orb_vir_lab_rev_up(ex_to_up(det_unq_up, ex_j))
          alphaI_inv_up(det_unq_up, ex_i, ex_j)=ainv_atilde_up(elec_i, orb_j)
        enddo
      enddo
      
      ![WEB]determinant
      alphaI_det_up(det_unq_up)=alphaI_inv_up(det_unq_up, 1, 1)

      ![WEB]inverse
      alphaI_inv_up(det_unq_up, 1, 1)=1.d0/alphaI_inv_up(det_unq_up, 1, 1)
    endif
  enddo

! dn
  call object_alloc('alphaI_inv_dn', alphaI_inv_dn, ndetdn, ex_order_dn_max, ex_order_dn_max)
  call object_alloc('alphaI_det_dn', alphaI_det_dn, ndetdn)
  do det_unq_dn=1, ndetdn
    n=ex_order_dn(det_unq_dn)
    if (n.eq.0) then
      alphaI_inv_dn(det_unq_dn, 1, 1)=1.d0
      alphaI_det_dn(det_unq_dn)=1.d0
    else
      do ex_i=1, n
        do ex_j=1, n
          elec_i=orb_occ_lab_rev_dn(ex_from_dn(det_unq_dn, ex_i))
          orb_j=orb_vir_lab_rev_dn(ex_to_dn(det_unq_dn, ex_j))
          alphaI_inv_dn(det_unq_dn, ex_i, ex_j)=ainv_atilde_dn(elec_i, orb_j)
        enddo
      enddo
      
      ![WEB]determinant
      alphaI_det_dn(det_unq_dn)=alphaI_inv_dn(det_unq_dn, 1, 1)

      ![WEB]inverse
      alphaI_inv_dn(det_unq_dn, 1, 1)=1.d0/alphaI_inv_dn(det_unq_dn, 1, 1)
    endif
  enddo

  end subroutine alphaI_det_inv_bld

! ==============================================================================
  subroutine ymat_ainv_and_yKmat_bld
! ------------------------------------------------------------------------------
! Description   : Compute the matrices Y_K (UP and DN) for each order of excitation
! Description   :  and Y.A^-1 (UP and DN)
!
! Created       : B. Mussard, Aug. 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer :: csf_i, det_i, det_in_csf_i, det_up, det_dn, ex_up, ex_dn
  integer :: ex_i, ex_j, elec_j, elec_k, orb_i
  real(dp):: common_coef

! header
  if (header_exe) then

   call object_create ('ymat_ainv_up')
   call object_create ('ymat_ainv_dn')
   call object_create ('yKmat_up')
   call object_create ('yKmat_dn')

   call object_needed('orb_occ_up_nb')
   call object_needed('orb_occ_dn_nb')
   call object_needed('orb_vir_up_nb')
   call object_needed('orb_vir_dn_nb')
   call object_needed('ex_order_up_max')
   call object_needed('ex_order_dn_max')
   call object_needed('ncsf')
   call object_needed('ndet_in_csf')
   call object_needed('iwdet_in_csf')
   call object_needed('det_to_det_unq_up')
   call object_needed('det_to_det_unq_dn')
   call object_needed('ex_order_up')
   call object_needed('ex_order_dn')
   call object_needed('csf_coef')
   call object_needed('cdet_in_csf')
   call object_needed('detu')
   call object_needed('detd')
   call object_needed('alphaI_det_up')
   call object_needed('alphaI_det_dn')
   call object_needed('orb_occ_lab_rev_up')
   call object_needed('orb_occ_lab_rev_dn')
   call object_needed('orb_vir_lab_rev_up')
   call object_needed('orb_vir_lab_rev_dn')
   call object_needed('alphaI_inv_up')
   call object_needed('alphaI_inv_dn')

   return

  endif

! allocations
  call object_alloc ('ymat_ainv_up',  ymat_ainv_up,  orb_vir_up_nb,   orb_occ_up_nb)
  call object_alloc ('ymat_ainv_dn',  ymat_ainv_dn,  orb_vir_dn_nb,   orb_occ_dn_nb)
  call object_alloc ('yKmat_up', yKmat_up, ex_order_up_max, orb_vir_up_nb, orb_occ_up_nb)
  call object_alloc ('yKmat_dn', yKmat_dn, ex_order_dn_max, orb_vir_dn_nb, orb_occ_dn_nb)
  ymat_ainv_up=0.d0
  ymat_ainv_dn=0.d0
  yKmat_up=0.d0
  yKmat_dn=0.d0

! yKmat
  do csf_i=1, ncsf
    do det_in_csf_i=1, ndet_in_csf(csf_i)
      det_i=iwdet_in_csf(det_in_csf_i, csf_i)
      det_up=det_to_det_unq_up(det_i)
      det_dn=det_to_det_unq_dn(det_i)
      ex_up=ex_order_up(det_up)
      ex_dn=ex_order_dn(det_dn)
      if ((ex_up.ne.0).and.(ex_dn.eq.0)) cycle

      common_coef=csf_coef(csf_i, 1)*cdet_in_csf(det_in_csf_i, csf_i) &
                 *detu(det_up)*detd(det_dn) &
                 *alphaI_det_up(det_up)*alphaI_det_dn(det_dn)
      write(6,*) 'common_coef', common_coef

      !up
      do ex_i=1, ex_up
        do ex_j=1, ex_up
          orb_i=orb_vir_lab_rev_up(ex_to_up  (det_up, ex_j))
          elec_j=orb_occ_lab_rev_up(ex_from_up(det_up, ex_i))
          yKmat_up(ex_up, orb_i, elec_j)= yKmat_up(ex_up, orb_i, elec_j) &
                                        +common_coef*alphaI_inv_up(det_up, ex_i, ex_j)
        enddo
      enddo
      !dn
      do ex_i=1, ex_dn
        do ex_j=1, ex_dn
          orb_i=orb_vir_lab_rev_dn(ex_to_dn(det_dn, ex_j))
          elec_j=orb_occ_lab_rev_dn(ex_from_dn(det_dn, ex_i))
          yKmat_dn(ex_dn, orb_i, elec_j)= yKmat_dn(ex_dn, orb_i, elec_j) &
                                        +common_coef*alphaI_inv_dn(det_dn, ex_i, ex_j)
        enddo
      enddo

    enddo
  enddo

  write(6,*) 'yKmat_up'
  do ex_up=1, ex_order_up_max
  do orb_i=1, orb_vir_up_nb
  write(6,*) (yKmat_up(ex_up, orb_i, elec_j), elec_j=1, orb_occ_up_nb)
  enddo
  enddo
  write(6,*) 'yKmat_dn'
  do ex_dn=1, ex_order_dn_max
  do orb_i=1, orb_vir_dn_nb
  write(6,*) (yKmat_dn(ex_dn, orb_i, elec_j), elec_j=1, orb_occ_dn_nb)
  enddo
  enddo

! ymat
  do ex_up=1, ex_order_up_max
    ymat_ainv_up=ymat_ainv_up+yKmat_up(ex_up,:,:)
  enddo
! [WEB: product]
  do orb_i=1,orb_vir_up_nb
  do elec_j= 1, orb_occ_up_nb
  do elec_k= 1, orb_occ_up_nb
    ymat_ainv_up(orb_i,elec_j)=ymat_ainv_up(orb_i,elec_j)+ymat_ainv_up(orb_i,elec_k)*slater_mat_trans_inv_up(elec_k,elec_j,det_ref_up)
  enddo
  enddo
  enddo

  do ex_dn=1, ex_order_dn_max
    ymat_ainv_dn=ymat_ainv_dn+yKmat_dn(ex_dn,:,:)
  enddo
! [WEB: product]
  do orb_i=1,orb_vir_dn_nb
  do elec_j= 1, orb_occ_dn_nb
  do elec_k= 1, orb_occ_dn_nb
    ymat_ainv_dn(orb_i,elec_j)=ymat_ainv_dn(orb_i,elec_j)+ymat_ainv_dn(orb_i,elec_k)*slater_mat_trans_inv_dn(elec_k,elec_j,det_ref_dn)
  enddo
  enddo
  enddo

  write(6,*) 'ymat_ainv_up'
  do orb_i=1, orb_vir_up_nb
  write(6,*) (ymat_ainv_up(orb_i, elec_j), elec_j=1, orb_occ_up_nb)
  enddo
  write(6,*) 'ymat_ainv_dn'
  do orb_i=1, orb_vir_dn_nb
  write(6,*) (ymat_ainv_dn(orb_i, elec_j), elec_j=1, orb_occ_dn_nb)
  enddo

  end subroutine ymat_ainv_and_yKmat_bld

! ==============================================================================
  subroutine myphi_from_yKmat_bld
! ------------------------------------------------------------------------------
! Description   : Compute Phi from Y_K (UP and DN) - this should give the same result
!
! Created       : B. Mussard, Aug. 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer :: i, j, k
  real(dp):: phi_up, phi_dn

! header
  if (header_exe) then

   call object_create ('myphi')

   call object_needed('ex_order_up_max')
   call object_needed('ex_order_dn_max')
   call object_needed('orb_occ_up_nb')
   call object_needed('orb_occ_dn_nb')
   call object_needed('orb_vir_up_nb')
   call object_needed('orb_vir_dn_nb')
   call object_needed('ainv_atilde_up')
   call object_needed('ainv_atilde_dn')
   call object_needed('yKmat_up')
   call object_needed('yKmat_dn')
   call object_needed('coef_ref')
   call object_needed('det_ref_up')
   call object_needed('det_ref_dn')
   call object_needed('detu')
   call object_needed('detd')

   return

  endif

! up
  phi_up=0.d0
  do k=1, ex_order_up_max
    do i=1, orb_occ_up_nb
      do j=1, orb_vir_up_nb
        phi_up=phi_up+ainv_atilde_up(i, j)*yKmat_up(k, j, i)/k
      enddo
    enddo
  enddo
  phi_up=phi_up+coef_ref*detu(det_ref_up)*detd(det_ref_dn)

! dn
  phi_dn=0.d0
  do k=1, ex_order_dn_max
    do i=1, orb_occ_dn_nb
      do j=1, orb_vir_dn_nb
        phi_dn=phi_dn+ainv_atilde_dn(i, j)*yKmat_dn(k, j, i)/k
      enddo
    enddo
  enddo
  phi_dn=phi_dn+coef_ref*detu(det_ref_up)*detd(det_ref_dn)

  write(6,*) 'phi-phi=', phi_up-phi_dn

  myphi=phi_up

  end subroutine myphi_from_yKmat_bld

! ==============================================================================
  subroutine gamma_bld
! ------------------------------------------------------------------------------
! Description   : Compute the matrices Gamma (UP and DN) 
!
! Created       : B. Mussard, Aug. 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer :: i, j, k

! header
  if (header_exe) then

   call object_create ('gamma_up')
   call object_create ('gamma_dn')

   call object_needed('orb_tot_nb')
   call object_needed('orb_occ_up_nb')
   call object_needed('orb_occ_dn_nb')
   call object_needed('orb_vir_up_nb')
   call object_needed('orb_vir_dn_nb')
   call object_needed('myphi')
   call object_needed('slater_mat_trans_inv_up')
   call object_needed('slater_mat_trans_inv_dn')
   call object_needed('det_ref_up')
   call object_needed('det_ref_dn')
   call object_needed('ainv_atilde_up')
   call object_needed('ainv_atilde_dn')
   call object_needed('ymat_ainv_up')
   call object_needed('ymat_ainv_dn')

   return

  endif

! up
  call object_alloc('gamma_up', gamma_up, orb_tot_nb, orb_occ_up_nb)
  do i=1, orb_occ_up_nb
    do j=1, orb_occ_up_nb
      gamma_up(i, j)=myphi*slater_mat_trans_inv_up(j, i, det_ref_up)
      do k=1, orb_vir_up_nb
        gamma_up(i, j)=gamma_up(i, j)-ainv_atilde_up(i, k)*ymat_ainv_up(k, j)
      enddo
    enddo
  enddo
  gamma_up(orb_occ_up_nb+1:,:)=ymat_ainv_up
  write(6,*) 'gamma_up'
  do i=1, orb_tot_nb
  write(6,*) (gamma_up(i,j),j=1,orb_occ_up_nb)
  enddo

! dn
  call object_alloc('gamma_dn', gamma_dn, orb_tot_nb, orb_occ_dn_nb)
  do i=1, orb_occ_dn_nb
    do j=1, orb_occ_dn_nb
      gamma_dn(i, j)=myphi*slater_mat_trans_inv_dn(j, i, det_ref_dn)
      do k=1, orb_vir_dn_nb
        gamma_dn(i, j)=gamma_dn(i, j)-ainv_atilde_dn(i, k)*ymat_ainv_dn(k, j)
      enddo
    enddo
  enddo
  gamma_dn(orb_occ_dn_nb+1:,:)=ymat_ainv_dn
  write(6,*) 'gamma_dn'
  do i=1, orb_tot_nb
  write(6,*) (gamma_dn(i,j),j=1,orb_occ_dn_nb)
  enddo

  end subroutine gamma_bld




















! ==============================================================================
  subroutine grd_det_over_det_fast_bld
! ------------------------------------------------------------------------------
! Description   : Grad(Psi)/Psi
!
! Created       : B. Mussard, Jan 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

!! local
!  integer dim_i, elec_i, elec_k
!
!! header
!  if (header_exe) then
!
!   call object_create ('grd_det_over_det_fast')
!
!   call object_needed('ndim')
!   call object_needed('nelec')
!   !call object_needed('nwalk')
!   call object_needed('dorb')
!   call object_needed('det_unq_orb_lab_srt_up')
!   call object_needed('det_unq_orb_lab_srt_dn')
!   call object_needed('slater_mat_trans_inv_up')
!   call object_needed('slater_mat_trans_inv_dn')
!   call object_needed('grd_det_over_det_legacy')
!
!   return
!
!  endif
!
!! begin
!
!! allocations
!  call object_alloc ('grd_det_over_det_fast', grd_det_over_det_fast, ndim, nelec, 1)
!
!  ! for all dimensions and electrons:
!  ! Grad(DET)/DET = -1/2 tr(A^-1.B)
!  ! where A is Slater determinant
!  ! and   B_ij=Grad_(dim, i)(\phi_j(r_i))
!  grd_det_over_det_fast = 0.d0
!  do dim_i=1, ndim
!    do elec_i=1, nup
!      do elec_k=1, nup
!        grd_det_over_det_fast(dim_i, elec_i, 1)=grd_det_over_det_fast(dim_i, elec_i, 1) &
!                                           & +slater_mat_trans_inv_up(elec_i, elec_k, 1)*dorb(dim_i, elec_i, det_unq_orb_lab_srt_up(elec_k, 1))
!      enddo
!    enddo
!    do elec_i=1, ndn
!      do elec_k=1, ndn
!        grd_det_over_det_fast(dim_i, nup+elec_i, 1)=grd_det_over_det_fast(dim_i, nup+elec_i, 1) &
!                                               & +slater_mat_trans_inv_dn(elec_i, elec_k, 1)*dorb(dim_i, nup+elec_i, det_unq_orb_lab_srt_dn(elec_k, 1))
!      enddo
!    enddo
!  enddo
!
!  call is_equal_or_die(grd_det_over_det_legacy, grd_det_over_det_fast(:,:, 1), 10.d-10, .false.)

  end subroutine grd_det_over_det_fast_bld

! ==============================================================================
  subroutine lap_det_over_det_fast_bld
! ------------------------------------------------------------------------------
! Description   : Kinetic local energy
!
! Created       : B. Mussard, Jan 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

!! local
!  integer elec_i, orb_j
!
!! header
!  if (header_exe) then
!
!   call object_create ('lap_det_over_det_fast')
!
!   call object_needed('nup')
!   call object_needed('ndn')
!   call object_needed('nwalk')
!   call object_needed('ddorb')
!   call object_needed('det_unq_orb_lab_srt_up')
!   call object_needed('det_unq_orb_lab_srt_dn')
!   call object_needed('slater_mat_trans_inv_up')
!   call object_needed('slater_mat_trans_inv_dn')
!   call object_needed('lap_det_over_det_legacy')
!
!   return
!
!  endif
!
!! begin
!
!! allocations
!  call object_alloc ('lap_det_over_det_fast', lap_det_over_det_fast, 1)
!  call object_alloc ('lap_det_over_det_fast_av', lap_det_over_det_fast_av, 1)
!  call object_alloc ('lap_det_over_det_fast_av_err', lap_det_over_det_fast_av_err, 1)
!  call object_associate ('lap_det_over_det_fast', lap_det_over_det_fast, 1)
!  call object_associate ('lap_det_over_det_fast_av', lap_det_over_det_fast_av, 1)
!  call object_associate ('lap_det_over_det_fast_av_err', lap_det_over_det_fast_av_err, 1)
!
!  ! T(DET)/DET = -1/2 tr(A^-1.B)
!  ! where A is Slater determinant
!  ! and   B_ij=\laplacian\phi_j(r_i)
!  lap_det_over_det_fast = 0.d0
!  do elec_i=1, nup
!    do orb_j=1, nup
!      lap_det_over_det_fast(1)=lap_det_over_det_fast(1) &
!                            & +slater_mat_trans_inv_up(elec_i, orb_j, 1)*ddorb(elec_i, det_unq_orb_lab_srt_up(orb_j, 1))
!    enddo
!  enddo
!  do elec_i=1, ndn
!    do orb_j=1, ndn
!      lap_det_over_det_fast(1)=lap_det_over_det_fast(1) &
!                            & +slater_mat_trans_inv_dn(elec_i, orb_j, 1)*ddorb(nup+elec_i, det_unq_orb_lab_srt_dn(orb_j, 1))
!    enddo
!  enddo
!  lap_det_over_det_fast(1)=-0.5d0*lap_det_over_det_fast(1)
!
!  if (abs(lap_det_over_det_fast(1)-lap_det_over_det_legacy).ge.1.d-10) then
!    write(6, '(a, f15.10)') 'is_equal_or_die_scalar:   checking equality ', lap_det_over_det_fast(1)-lap_det_over_det_legacy
!  endif

  end subroutine lap_det_over_det_fast_bld

end module derivatives_fast_mod

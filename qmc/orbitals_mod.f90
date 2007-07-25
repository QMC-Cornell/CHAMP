module orbitals_mod

  use all_tools_mod
  use basis_mod

! Declaration of global variables and default values
  integer                             :: orb_tot_nb

  logical                             :: l_cusp_en     = .false.
  logical                             :: l_cusp_en_occ = .false.
  logical                             :: l_cusp_en_opt = .false.
  logical                             :: l_approx_orb_rot = .false.
  logical                             :: l_ortho_orb_now = .false.
  logical                             :: l_ortho_orb_opt = .false.

  integer, allocatable                :: det_unq_orb_lab_srt_up (:,:)
  integer, allocatable                :: det_unq_orb_lab_srt_dn (:,:)
  integer, allocatable                :: det_to_det_unq_up (:)
  integer, allocatable                :: det_to_det_unq_dn (:)

  integer                             :: orb_occ_in_wf_nb  = 0
  integer                             :: orb_cls_in_wf_nb  = 0
  integer                             :: orb_act_in_wf_nb  = 0
  integer                             :: orb_opn_in_wf_nb  = 0
  integer                             :: orb_vir_in_wf_nb  = 0

  logical, allocatable                :: orb_occ_in_det_unq_up (:,:)
  logical, allocatable                :: orb_occ_in_det_unq_dn (:,:)
  integer, allocatable                :: orb_pos_in_det_unq_up (:,:)
  integer, allocatable                :: orb_pos_in_det_unq_dn (:,:)
  logical, allocatable                :: orb_occ_in_wf (:)
  logical, allocatable                :: orb_cls_in_wf (:)
  logical, allocatable                :: orb_act_in_wf (:)
  logical, allocatable                :: orb_opn_in_wf (:)
  logical, allocatable                :: orb_vir_in_wf (:)
  integer, allocatable                :: orb_occ_in_wf_lab (:)
  integer, allocatable                :: orb_cls_in_wf_lab (:)
  integer, allocatable                :: orb_act_in_wf_lab (:)
  integer, allocatable                :: orb_opn_in_wf_lab (:)
  integer, allocatable                :: orb_vir_in_wf_lab (:)

  integer                             :: orb_opt_nb  = 0
  integer, allocatable                :: orb_opt_lab (:)
  integer                             :: orb_opt_occ_nb  = 0
  integer, allocatable                :: orb_opt_occ_lab (:)
  integer                             :: orb_opt_cls_nb  = 0
  integer, allocatable                :: orb_opt_cls_lab (:)
  integer                             :: orb_opt_opn_nb  = 0
  integer, allocatable                :: orb_opt_opn_lab (:)
  integer                             :: orb_opt_act_nb  = 0
  integer, allocatable                :: orb_opt_act_lab (:)
  integer                             :: orb_opt_vir_nb  = 0
  integer, allocatable                :: orb_opt_vir_lab (:)
  integer                             :: orb_opt_last_lab

  character(len=max_string_len), allocatable  :: orb_sym_lab (:)
  real(dp), allocatable               :: orb_energies (:)
  real(dp), allocatable               :: orb_ovlp (:,:)
  real(dp), allocatable               :: orb_cls_ovlp (:,:)
  real(dp), allocatable               :: orb_occ_ovlp (:,:)
  real(dp), allocatable               :: orb_act_ovlp (:,:)
  real(dp), allocatable               :: orb_vir_ovlp (:,:)
  real(dp), allocatable               :: orb_cls_ovlp_inv (:,:)
  real(dp), allocatable               :: orb_occ_ovlp_inv (:,:)
  real(dp), allocatable               :: orb_cls_ovlp_eigvec (:,:)
  real(dp), allocatable               :: orb_act_ovlp_eigvec (:,:)
  real(dp), allocatable               :: orb_vir_ovlp_eigvec (:,:)
  real(dp), allocatable               :: orb_cls_ovlp_eigval (:)
  real(dp), allocatable               :: orb_act_ovlp_eigval (:)
  real(dp), allocatable               :: orb_vir_ovlp_eigval (:)
  real(dp), allocatable               :: orb_cls_ovlp_m12 (:,:)
  real(dp), allocatable               :: orb_act_ovlp_m12 (:,:)
  real(dp), allocatable               :: orb_vir_ovlp_m12 (:,:)
  real(dp), allocatable               :: coef_orb_on_norm_basis (:,:,:)
  real(dp), allocatable               :: coef_orb_on_ortho_basis (:,:,:)
  real(dp), allocatable               :: coef_sav (:,:)
  real(dp), allocatable               :: coef_orb_on_norm_basis_sav (:,:)
  real(dp), allocatable               :: coef_orb_on_ortho_basis_sav (:,:)
  logical, allocatable                :: is_orb_s (:)

  logical, allocatable                :: orb_ex_forbidden (:,:)

! orbital files
  character (len=max_string_len_file) :: file_orbitals_pw_in     = 'orbitals_pw'
  character (len=max_string_len_file) :: file_orbitals_pw_out    = 'orbitals_pw.out'
  character (len=max_string_len_file) :: file_orbitals_pw_tm_in  = 'orbitals_pw_tm'
  character (len=max_string_len_file) :: file_orbitals_pw_tm_out = 'orbitals_pw_tm.out'

  contains

!===========================================================================
  subroutine orbitals_menu
!---------------------------------------------------------------------------
! Description : menu for orbitals
!
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'orbitals_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for orbitals menu:'
   write(6,'(a)') ' orbitals'
   write(6,'(a)') '  energies -2.4 -0.2 0.44 7.4 end: orbital energies for perturbative optimization method'
   write(6,'(a)') '  symmetry A1G A2U EG EU end: symmetry labels for orbitals, default: no symmetry'
   write(6,'(a)') '  orthonormalize = [bool] orthonormalize orbitals once of current wave function? (default=false)'
   write(6,'(a)') '  opt 1 2 3 4 end : which (occupied and virtual) orbitals to consider in the optimization, default: all orbitals'
   write(6,'(a)') '  orb_opt_nb = [integer] : total number of orbitals to consider in the optimization, default: total number of orbitals'
   write(6,'(a)') '  excitations_forbidden  6 7   10 23 end : list of forbidden orbital excitations for orbital optimization'
   write(6,'(a)') '  approx_orb_rot = [bool] update orbitals by approximate first-order rotation (default=false)'
   write(6,'(a)') '  cusp  ... end : menu to impose e-N cusp on orbitals (does not work presently)'
   write(6,'(a)') '  file_orbitals_pw_in  = [string] : input file for real orbitals on plane waves (default=orbitals_pw)'
   write(6,'(a)') '  file_orbitals_pw_out = [string] : output file for real orbitals on plane waves (default=orbitals_pw.out)'
   write(6,'(a)') '  file_orbitals_pw_tm_in  = [string] : input file for orbitals on plane waves in TM format (default=orbitals_pw_tm)'
   write(6,'(a)') '  file_orbitals_pw_tm_out = [string] : output file for orbitals on plane waves in TM format (default=orbitals_pw_tm.out)'
   write(6,'(a)') ' end'
   write(6,*)

  case ('energies')
   call orb_energies_rd

  case ('symmetry')
   call orb_sym_lab_rd

  case ('orthonormalize')
   call get_next_value (l_ortho_orb_now)

  case ('cusp')
   call orb_cusp_menu

  case ('opt')
   call orb_opt_lab_rd

  case ('orb_opt_nb')
   call get_next_value (orb_opt_nb)
   call object_modified ('orb_opt_nb')

  case ('excitations_forbidden')
   call orb_ex_forbidden_rd

  case ('file_orbitals_pw_in')
   call get_next_value (file_orbitals_pw_in)

  case ('file_orbitals_pw_out')
   call get_next_value (file_orbitals_pw_out)

  case ('file_orbitals_pw_tm_in')
   call get_next_value (file_orbitals_pw_tm_in)

  case ('file_orbitals_pw_tm_out')
   call get_next_value (file_orbitals_pw_tm_out)

  case ('approx_orb_rot')
   call get_next_value (l_approx_orb_rot)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

! energy-invariant orthonormalization the orbitals
  if (l_ortho_orb_now) then
    call ortho_orb
  endif

  end subroutine orbitals_menu

!===========================================================================
  subroutine orb_cusp_menu
!---------------------------------------------------------------------------
! Description : menu for e-N cusp on orbitals
!
! Created     : J. Toulouse, 13 Jan 2006
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

  character*(max_string_len_rout),save :: lhere = 'orb_cusp_menu'
  character*(max_string_len) line, word

! begin
  l_cusp_en = .true.
  write(6,*) trim(lhere), ': l_cusp_en=', l_cusp_en

! loop over menu lines
  do
  call get_next_word (word)

  if(trim(word) == 'help') then
   write(6,*) 'cusp-h: menu for e-N cusp'
   write(6,*) 'cusp-h: cusp'
   write(6,*) 'cusp-h:  occ: impose the e-N cusp conditions on occupied orbitals only'
   write(6,*) 'cusp-h:  opt: impose the e-N cusp conditions during the orbital optimization'
   write(6,*) 'cusp-h: end'

  elseif(trim(word) == 'occ') then
   l_cusp_en_occ = .true.

  elseif(trim(word) == 'opt') then
   l_cusp_en_opt = .true.

  elseif(trim(word) == 'end') then
   exit

  else

   write(6,*) trim(lhere),': unknown word = ',trim(word)
   call die(lhere)

  endif

  enddo ! end loop over menu lines

! impose the e-N cusp condition
  write(6,*) trim(lhere), ': call cusp_en_orb'
  call cusp_en_orb

  end subroutine orb_cusp_menu

! ==============================================================================
  subroutine det_unq_orb_lab_srt_bld
! ------------------------------------------------------------------------------
! Description   : unique spin-up and dn determinants represented as sorted list of orbital labels
!
! Created       : J. Toulouse, 28 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_unq_up_i, det_unq_dn_i
  integer elec_up_i, elec_dn_i

! header
  if (header_exe) then

   call object_create ('det_unq_orb_lab_srt_up')
   call object_create ('det_unq_orb_lab_srt_dn')
   call object_create ('det_to_det_unq_up')
   call object_create ('det_to_det_unq_dn')

   call object_needed ('ndet')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('iworbdup')
   call object_needed ('iworbddn')
   call object_needed ('iwdetup')
   call object_needed ('iwdetdn')

   return

  endif

! begin

! allocations
  call object_alloc ('det_unq_orb_lab_srt_up', det_unq_orb_lab_srt_up, nup, ndetup)
  call object_alloc ('det_unq_orb_lab_srt_dn', det_unq_orb_lab_srt_dn, ndn, ndetdn)
  call object_alloc ('det_to_det_unq_up', det_to_det_unq_up, ndet)
  call object_alloc ('det_to_det_unq_dn', det_to_det_unq_dn, ndet)

! unique spin-up determinants
  det_unq_orb_lab_srt_up (1:nup, 1:ndetup) = iworbdup (1:nup, 1:ndetup)
  det_unq_orb_lab_srt_dn (1:ndn, 1:ndetdn) = iworbddn (1:ndn, 1:ndetdn)

! from determinants to unique determinants
  det_to_det_unq_up (1:ndet) = iwdetup (1:ndet)
  det_to_det_unq_dn (1:ndet) = iwdetdn (1:ndet)

! imposing that the determinants must be given already sorted to avoid taking care of permutations elsewehere
  do det_unq_up_i = 1, ndetup
    if (.not. is_sorted (det_unq_orb_lab_srt_up (:, det_unq_up_i))) then
        write(6,'(a,i,a,100i)') 'orbitals in unique spin-up determinant # ', det_unq_up_i,' are not in increasing order:',det_unq_orb_lab_srt_up (:, det_unq_up_i)
        call die (here)
    endif
  enddo
  do det_unq_dn_i = 1, ndetdn
    if (.not. is_sorted (det_unq_orb_lab_srt_dn (:, det_unq_dn_i))) then
        write(6,'(a,i,a,100i)') 'orbitals in unique spin-dn determinant # ', det_unq_dn_i,' are not in increasing order:',det_unq_orb_lab_srt_dn (:, det_unq_dn_i)
        call die (here)
    endif
  enddo

  end subroutine det_unq_orb_lab_srt_bld

! ==============================================================================
  subroutine orb_occupations_bld
! ------------------------------------------------------------------------------
! Description   : build objects giving the orbitals occupation
!
! Created       : J. Toulouse, 20 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_i, det_unq_up_i, det_unq_dn_i, ielec, orb_i
  integer elec_up_i, elec_dn_i
  integer orb_occ_in_wf_i, orb_cls_in_wf_i, orb_act_in_wf_i, orb_opn_in_wf_i, orb_vir_in_wf_i
  integer orb_up_lab, orb_dn_lab

! header
  if (header_exe) then

   call object_create ('orb_occ_in_det_unq_up')
   call object_create ('orb_occ_in_det_unq_dn')
   call object_create ('orb_pos_in_det_unq_up')
   call object_create ('orb_pos_in_det_unq_dn')
   call object_create ('orb_occ_in_wf')
   call object_create ('orb_cls_in_wf')
   call object_create ('orb_act_in_wf')
   call object_create ('orb_opn_in_wf')
   call object_create ('orb_vir_in_wf')
   call object_create ('orb_occ_in_wf_nb')
   call object_create ('orb_cls_in_wf_nb')
   call object_create ('orb_act_in_wf_nb')
   call object_create ('orb_opn_in_wf_nb')
   call object_create ('orb_vir_in_wf_nb')
   call object_create ('orb_occ_in_wf_lab')
   call object_create ('orb_cls_in_wf_lab')
   call object_create ('orb_act_in_wf_lab')
   call object_create ('orb_opn_in_wf_lab')
   call object_create ('orb_vir_in_wf_lab')
   call object_create ('orb_occ_last_in_wf_lab')

   call object_needed ('ndet')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('orb_tot_nb')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('orb_occ_in_det_unq_up', orb_occ_in_det_unq_up, orb_tot_nb, ndetup)
  call object_alloc ('orb_occ_in_det_unq_dn', orb_occ_in_det_unq_dn, orb_tot_nb, ndetdn)
  call object_alloc ('orb_pos_in_det_unq_up', orb_pos_in_det_unq_up, orb_tot_nb, ndetup)
  call object_alloc ('orb_pos_in_det_unq_dn', orb_pos_in_det_unq_dn, orb_tot_nb, ndetdn)
  call object_alloc ('orb_occ_in_wf', orb_occ_in_wf, orb_tot_nb)
  call object_alloc ('orb_cls_in_wf', orb_cls_in_wf, orb_tot_nb)
  call object_alloc ('orb_act_in_wf', orb_act_in_wf, orb_tot_nb)
  call object_alloc ('orb_opn_in_wf', orb_opn_in_wf, orb_tot_nb)
  call object_alloc ('orb_vir_in_wf', orb_vir_in_wf, orb_tot_nb)

  orb_occ_in_det_unq_up = .false.
  orb_occ_in_det_unq_dn = .false.
  orb_pos_in_det_unq_up = 0
  orb_pos_in_det_unq_dn = 0

! loop over unique spin-up determinants
  do det_unq_up_i = 1, ndetup
    do elec_up_i = 1, nup
      orb_up_lab = det_unq_orb_lab_srt_up (elec_up_i, det_unq_up_i)
      orb_occ_in_det_unq_up (orb_up_lab, det_unq_up_i) = .true.
      orb_pos_in_det_unq_up (orb_up_lab, det_unq_up_i) = elec_up_i
    enddo ! elec_up_i
  enddo ! det_unq_up_i

! loop over unique spin-dn determinants
  do det_unq_dn_i = 1, ndetdn
    do elec_dn_i = 1, ndn
      orb_dn_lab = det_unq_orb_lab_srt_dn (elec_dn_i, det_unq_dn_i)
      orb_occ_in_det_unq_dn (orb_dn_lab, det_unq_dn_i) = .true.
      orb_pos_in_det_unq_dn (orb_dn_lab, det_unq_dn_i) = elec_dn_i
    enddo ! elec_dn_i
  enddo ! det_unq_dn_i

! determine:
! occ orbitals in at least one determinant
! cls orbitals in all determinants
! opn orbitals in at least one determinant
! vir orbitals in all determinants

  orb_occ_in_wf = .false.
  orb_cls_in_wf = .true.
  orb_opn_in_wf = .false.
  orb_vir_in_wf = .true.

! loop over all determinants
  do det_i = 1, ndet

!  unique spin-up and spin-dn determinants
   det_unq_up_i = det_to_det_unq_up (det_i)
   det_unq_dn_i = det_to_det_unq_dn (det_i)

!  loop over all orbitals
   do orb_i = 1, orb_tot_nb

      if (.not. orb_occ_in_det_unq_up (orb_i, det_unq_up_i) .or. .not. orb_occ_in_det_unq_dn (orb_i, det_unq_dn_i)) then
        orb_cls_in_wf (orb_i) = .false.
      endif

      if (orb_occ_in_det_unq_up (orb_i, det_unq_up_i) .or. orb_occ_in_det_unq_dn (orb_i, det_unq_dn_i)) then
        orb_occ_in_wf (orb_i) = .true.
        orb_vir_in_wf (orb_i)  = .false.
      endif

      if ( (orb_occ_in_det_unq_up (orb_i, det_unq_up_i) .and. .not. orb_occ_in_det_unq_dn (orb_i, det_unq_dn_i) )  &
         .or.  (.not. orb_occ_in_det_unq_up (orb_i, det_unq_up_i) .and. orb_occ_in_det_unq_dn (orb_i, det_unq_dn_i) ) ) then
        orb_opn_in_wf (orb_i) = .true.
      endif

    enddo ! orb_i
  enddo ! det_i

  orb_occ_in_wf_nb =  0
  orb_cls_in_wf_nb =  0
  orb_opn_in_wf_nb =  0
  orb_vir_in_wf_nb =  0

  do orb_i = 1, orb_tot_nb

    if (orb_occ_in_wf (orb_i)) then
      orb_occ_in_wf_nb = orb_occ_in_wf_nb + 1
    endif
    if (orb_cls_in_wf (orb_i)) then
      orb_cls_in_wf_nb = orb_cls_in_wf_nb + 1
    endif
    if (orb_opn_in_wf (orb_i)) then
      orb_opn_in_wf_nb = orb_opn_in_wf_nb + 1
    endif
    if (orb_vir_in_wf (orb_i)) then
      orb_vir_in_wf_nb = orb_vir_in_wf_nb + 1
    endif

  enddo ! orb_i

! determine:
! act orbitals = not cls orbitals, occ in at least one determinant
  orb_act_in_wf_nb =  0
  orb_act_in_wf = .false.
  do orb_i = 1, orb_tot_nb
    if (.not. orb_cls_in_wf (orb_i) .and. orb_occ_in_wf (orb_i)) then
       orb_act_in_wf (orb_i) = .true.
       orb_act_in_wf_nb = orb_act_in_wf_nb + 1
    endif
  enddo

! labels of orbitals
  call object_alloc ('orb_occ_in_wf_lab', orb_occ_in_wf_lab, orb_occ_in_wf_nb)
  call object_alloc ('orb_cls_in_wf_lab', orb_cls_in_wf_lab, orb_cls_in_wf_nb)
  call object_alloc ('orb_act_in_wf_lab', orb_act_in_wf_lab, orb_act_in_wf_nb)
  call object_alloc ('orb_opn_in_wf_lab', orb_opn_in_wf_lab, orb_opn_in_wf_nb)
  call object_alloc ('orb_vir_in_wf_lab', orb_vir_in_wf_lab, orb_vir_in_wf_nb)

  orb_occ_in_wf_i = 0
  orb_cls_in_wf_i = 0
  orb_act_in_wf_i = 0
  orb_opn_in_wf_i = 0
  orb_vir_in_wf_i = 0

  do orb_i = 1, orb_tot_nb

    if (orb_occ_in_wf (orb_i)) then
     orb_occ_in_wf_i = orb_occ_in_wf_i + 1
     orb_occ_in_wf_lab (orb_occ_in_wf_i) = orb_i
     orb_occ_last_in_wf_lab = orb_i
    endif

    if (orb_cls_in_wf (orb_i)) then
     orb_cls_in_wf_i = orb_cls_in_wf_i + 1
     orb_cls_in_wf_lab (orb_cls_in_wf_i) = orb_i
    endif

    if (orb_act_in_wf (orb_i)) then
     orb_act_in_wf_i = orb_act_in_wf_i + 1
     orb_act_in_wf_lab (orb_act_in_wf_i) = orb_i
    endif

    if (orb_opn_in_wf (orb_i)) then
     orb_opn_in_wf_i = orb_opn_in_wf_i + 1
     orb_opn_in_wf_lab (orb_opn_in_wf_i) = orb_i
    endif

    if (orb_vir_in_wf (orb_i)) then
     orb_vir_in_wf_i = orb_vir_in_wf_i + 1
     orb_vir_in_wf_lab (orb_vir_in_wf_i) = orb_i
    endif

  enddo ! orb_i

  if (orb_occ_in_wf_i /= orb_occ_in_wf_nb) then
    call die (here, 'orb_occ_in_wf_i='+orb_occ_in_wf_i+' /= orb_occ_in_wf_nb='+orb_occ_in_wf_nb)
  endif

  if (orb_cls_in_wf_i /= orb_cls_in_wf_nb) then
    call die (here, 'orb_cls_in_wf_i='+orb_cls_in_wf_i+' /= orb_cls_in_wf_nb='+orb_cls_in_wf_nb)
  endif

  if (orb_act_in_wf_i /= orb_act_in_wf_nb) then
    call die (here, 'orb_act_in_wf_i='+orb_act_in_wf_i+' /= orb_act_in_wf_nb='+orb_act_in_wf_nb)
  endif

  if (orb_opn_in_wf_i /= orb_opn_in_wf_nb) then
    call die (here, 'orb_opn_in_wf_i='+orb_opn_in_wf_i+' /= orb_opn_in_wf_nb='+orb_opn_in_wf_nb)
  endif

  if (orb_vir_in_wf_i /= orb_vir_in_wf_nb) then
    call die (here, 'orb_vir_in_wf_i='+orb_vir_in_wf_i+' /= orb_vir_in_wf_nb='+orb_vir_in_wf_nb)
  endif

  write(6,*)
  write(6,'(a)') 'Orbital occupation information:'
  write(6,'(a,i3,a      )') 'There are ',       orb_tot_nb,' total    orbitals'
  write(6,'(a,i3,a,500i4)') 'There are ', orb_occ_in_wf_nb,' occupied orbitals of labels:', orb_occ_in_wf_lab
  write(6,'(a,i3,a,500i4)') 'There are ', orb_cls_in_wf_nb,' closed   orbitals of labels:', orb_cls_in_wf_lab
  write(6,'(a,i3,a,500i4)') 'There are ', orb_act_in_wf_nb,' active   orbitals of labels:', orb_act_in_wf_lab
  write(6,'(a,i3,a,500i4)') 'There are ', orb_opn_in_wf_nb,' open     orbitals of labels:', orb_opn_in_wf_lab
  write(6,'(a,i3,a,500i4)') 'There are ', orb_vir_in_wf_nb,' virtual  orbitals of labels:', orb_vir_in_wf_lab

  if (orb_occ_in_wf_nb <= 0) then
    call die (here, 'number of occupied orbitals ='+ orb_occ_in_wf_nb+' is <= 0')
  endif

! label of last occupied orbital
  write(6,'(a,i3)') 'Index of last occupied orbital: ', orb_occ_last_in_wf_lab

  if (orb_occ_last_in_wf_lab <= 0) then
    call die (here, 'label of last occupied orbital ='+ orb_occ_last_in_wf_lab+' is <= 0')
  endif

  end subroutine orb_occupations_bld

!===========================================================================
  subroutine orb_energies_rd
!---------------------------------------------------------------------------
! Description : read orbitals eigenvalues
!
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'orb_energies_rd'
  integer orb_energies_nb, orb_i

! begin
  call get_next_value_list ('orb_energies', orb_energies, orb_energies_nb)

  call object_provide ('orb_tot_nb')
  if (orb_energies_nb /= orb_tot_nb) then
   call die (lhere, 'orb_energies_nb='+orb_energies_nb+ ' /= orb_tot_nb='+orb_tot_nb)
  endif

!  do orb_i = 1, orb_tot_nb
!   write(6,*) trim(here),': orb_i=',orb_i,' orb_energies=',orb_energies(orb_i)
!  enddo

  call object_modified ('orb_energies')

  end subroutine orb_energies_rd

!===========================================================================
  subroutine orb_sym_lab_rd
!---------------------------------------------------------------------------
! Description : read orbitals symmetry labels
!
! Created     : J. Toulouse, 01 Jan 2006
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'orb_sym_lab_rd'
  integer orb_sym_lab_nb, orb_i

! begin
  orb_sym_lab_nb  = 0
  call get_next_value_list ('orb_sym_lab', orb_sym_lab, orb_sym_lab_nb)

  call object_provide ('orb_tot_nb')
  if (orb_sym_lab_nb /= orb_tot_nb) then
   call die (lhere, 'orb_sym_lab_nb='+orb_sym_lab_nb+ ' /= orb_tot_nb='+orb_tot_nb)
  endif

  call object_modified ('orb_sym_lab')

  end subroutine orb_sym_lab_rd

! ==============================================================================
  subroutine orb_sym_lab_default_bld
! ------------------------------------------------------------------------------
! Description   : build default orbital symmetry labels (A A A A A A A ...)
! Description   : if the actual symmetry labels are not specified in the input
!
! Created       : J. Toulouse, 01 Jan 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('orb_sym_lab')

   call object_needed ('orb_tot_nb')

   return

  endif

! begin

! requirements
  if (orb_tot_nb <= 0) then
   call die (here, 'orb_tot_nb'+orb_tot_nb+' <= 0')
  endif

! allocations
  call object_alloc ('orb_sym_lab', orb_sym_lab, orb_tot_nb)

  orb_sym_lab (:) = 'A'

  write (6,'(2a)') trim(here), ': warning: not using symmetry for orbitals'

  end subroutine orb_sym_lab_default_bld

! ==============================================================================
  subroutine orb_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : total overlap matrix of all orbitals
!
! Created       : J. Toulouse, 17 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer bas_i, bas_j, orb_i, orb_j

! header
  if (header_exe) then

   call object_create ('orb_ovlp')

   call object_needed ('orb_tot_nb')
   call object_needed ('nbasis')
   call object_needed ('coef')
   call object_needed ('basis_ovlp')

   return

  endif

! begin

! requirements
  if (inum_orb /= 0 ) then
   call die (here, 'implemented only for analytical orbitals, i.e inum_orb=0')
  endif
! end requirements

! allocations
  call object_alloc ('orb_ovlp', orb_ovlp, orb_tot_nb, orb_tot_nb)

  do orb_i = 1, orb_tot_nb
    do orb_j = orb_i, orb_tot_nb

     orb_ovlp (orb_i, orb_j) = 0.d0

      do bas_i = 1, nbasis
        do bas_j = 1, nbasis
          orb_ovlp (orb_i, orb_j) = orb_ovlp (orb_i, orb_j) + coef (bas_i, orb_i, 1) * coef (bas_j, orb_j, 1) * basis_ovlp (bas_i, bas_j)
         enddo ! bas_j
       enddo ! bas_i

      orb_ovlp (orb_j, orb_i) = orb_ovlp (orb_i, orb_j)

    enddo ! orb_j
  enddo ! orb_i

  end subroutine orb_ovlp_bld

! ==============================================================================
  subroutine orb_cls_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap matrix of closed orbitals
!
! Created       : J. Toulouse, 01 Jun 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, orb_j, orb_cls_in_wf_i, orb_cls_in_wf_j, bas_i, bas_j

! header
  if (header_exe) then

   call object_create ('orb_cls_ovlp')

   call object_needed ('orb_cls_in_wf_nb')
   call object_needed ('orb_cls_in_wf_lab')
   call object_needed ('nbasis')
   call object_needed ('coef')
   call object_needed ('basis_ovlp')

   return

  endif

! begin

! requirements
  if (inum_orb /= 0 ) then
   call die (here, 'implemented only for analytical orbitals, i.e inum_orb=0')
  endif
! end requirements

! allocations
  call object_alloc ('orb_cls_ovlp', orb_cls_ovlp, orb_cls_in_wf_nb, orb_cls_in_wf_nb)

  do orb_cls_in_wf_i = 1, orb_cls_in_wf_nb
     orb_i = orb_cls_in_wf_lab (orb_cls_in_wf_i)

    do orb_cls_in_wf_j = orb_cls_in_wf_i, orb_cls_in_wf_nb
       orb_j = orb_cls_in_wf_lab (orb_cls_in_wf_j)

       orb_cls_ovlp (orb_cls_in_wf_i, orb_cls_in_wf_j) = 0.d0

       do bas_i = 1, nbasis
        do bas_j = 1, nbasis
          orb_cls_ovlp (orb_cls_in_wf_i, orb_cls_in_wf_j) = orb_cls_ovlp (orb_cls_in_wf_i, orb_cls_in_wf_j) + coef (bas_i, orb_i, 1) * coef (bas_j, orb_j, 1) * basis_ovlp (bas_i, bas_j)
         enddo ! bas_j
       enddo ! bas_i

       orb_cls_ovlp (orb_cls_in_wf_j, orb_cls_in_wf_i) = orb_cls_ovlp (orb_cls_in_wf_i, orb_cls_in_wf_j)

    enddo ! orb_cls_in_wf_j
  enddo ! orb_cls_in_wf_i

  end subroutine orb_cls_ovlp_bld

! ==============================================================================
  subroutine orb_occ_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap matrix of occupied orbitals
!
! Created       : J. Toulouse, 22 May 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, orb_j, orb_occ_in_wf_i, orb_occ_in_wf_j, bas_i, bas_j

! header
  if (header_exe) then

   call object_create ('orb_occ_ovlp')

   call object_needed ('orb_occ_in_wf_nb')
   call object_needed ('orb_occ_in_wf_lab')
   call object_needed ('nbasis')
   call object_needed ('coef')
   call object_needed ('basis_ovlp')

   return

  endif

! begin

! requirements
  if (inum_orb /= 0 ) then
   call die (here, 'implemented only for analytical orbitals, i.e inum_orb=0')
  endif
! end requirements

! allocations
  call object_alloc ('orb_occ_ovlp', orb_occ_ovlp, orb_occ_in_wf_nb, orb_occ_in_wf_nb)

  do orb_occ_in_wf_i = 1, orb_occ_in_wf_nb
     orb_i = orb_occ_in_wf_lab (orb_occ_in_wf_i)

    do orb_occ_in_wf_j = orb_occ_in_wf_i, orb_occ_in_wf_nb
       orb_j = orb_occ_in_wf_lab (orb_occ_in_wf_j)

       orb_occ_ovlp (orb_occ_in_wf_i, orb_occ_in_wf_j) = 0.d0

       do bas_i = 1, nbasis
        do bas_j = 1, nbasis
          orb_occ_ovlp (orb_occ_in_wf_i, orb_occ_in_wf_j) = orb_occ_ovlp (orb_occ_in_wf_i, orb_occ_in_wf_j) + coef (bas_i, orb_i, 1) * coef (bas_j, orb_j, 1) * basis_ovlp (bas_i, bas_j)
         enddo ! bas_j
       enddo ! bas_i

       orb_occ_ovlp (orb_occ_in_wf_j, orb_occ_in_wf_i) = orb_occ_ovlp (orb_occ_in_wf_i, orb_occ_in_wf_j)

    enddo ! orb_occ_in_wf_j
  enddo ! orb_occ_in_wf_i

  end subroutine orb_occ_ovlp_bld

! ==============================================================================
  subroutine orb_act_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap matrix of active orbitals
!
! Created       : J. Toulouse, 23 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, orb_j, orb_act_in_wf_i, orb_act_in_wf_j, bas_i, bas_j

! header
  if (header_exe) then

   call object_create ('orb_act_ovlp')

   call object_needed ('orb_act_in_wf_nb')
   call object_needed ('orb_act_in_wf_lab')
   call object_needed ('nbasis')
   call object_needed ('coef')
   call object_needed ('basis_ovlp')

   return

  endif

! begin

! requirements
  if (inum_orb /= 0 ) then
   call die (here, 'implemented only for analytical orbitals, i.e inum_orb=0')
  endif
! end requirements

! allocations
  call object_alloc ('orb_act_ovlp', orb_act_ovlp, orb_act_in_wf_nb, orb_act_in_wf_nb)

  do orb_act_in_wf_i = 1, orb_act_in_wf_nb
     orb_i = orb_act_in_wf_lab (orb_act_in_wf_i)

    do orb_act_in_wf_j = orb_act_in_wf_i, orb_act_in_wf_nb
       orb_j = orb_act_in_wf_lab (orb_act_in_wf_j)

       orb_act_ovlp (orb_act_in_wf_i, orb_act_in_wf_j) = 0.d0

       do bas_i = 1, nbasis
        do bas_j = 1, nbasis
          orb_act_ovlp (orb_act_in_wf_i, orb_act_in_wf_j) = orb_act_ovlp (orb_act_in_wf_i, orb_act_in_wf_j) + coef (bas_i, orb_i, 1) * coef (bas_j, orb_j, 1) * basis_ovlp (bas_i, bas_j)
         enddo ! bas_j
       enddo ! bas_i

       orb_act_ovlp (orb_act_in_wf_j, orb_act_in_wf_i) = orb_act_ovlp (orb_act_in_wf_i, orb_act_in_wf_j)

    enddo ! orb_act_in_wf_j
  enddo ! orb_act_in_wf_i

  end subroutine orb_act_ovlp_bld

! ==============================================================================
  subroutine orb_vir_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap matrix of virtual orbitals
!
! Created       : J. Toulouse, 01 Jun 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, orb_j, orb_vir_in_wf_i, orb_vir_in_wf_j, bas_i, bas_j

! header
  if (header_exe) then

   call object_create ('orb_vir_ovlp')

   call object_needed ('orb_vir_in_wf_nb')
   call object_needed ('orb_vir_in_wf_lab')
   call object_needed ('nbasis')
   call object_needed ('coef')
   call object_needed ('basis_ovlp')

   return

  endif

! begin

! requirements
  if (inum_orb /= 0 ) then
   call die (here, 'implemented only for analytical orbitals, i.e inum_orb=0')
  endif
! end requirements

! allocations
  call object_alloc ('orb_vir_ovlp', orb_vir_ovlp, orb_vir_in_wf_nb, orb_vir_in_wf_nb)

  do orb_vir_in_wf_i = 1, orb_vir_in_wf_nb
     orb_i = orb_vir_in_wf_lab (orb_vir_in_wf_i)

    do orb_vir_in_wf_j = orb_vir_in_wf_i, orb_vir_in_wf_nb
       orb_j = orb_vir_in_wf_lab (orb_vir_in_wf_j)

       orb_vir_ovlp (orb_vir_in_wf_i, orb_vir_in_wf_j) = 0.d0

       do bas_i = 1, nbasis
        do bas_j = 1, nbasis
          orb_vir_ovlp (orb_vir_in_wf_i, orb_vir_in_wf_j) = orb_vir_ovlp (orb_vir_in_wf_i, orb_vir_in_wf_j) + coef (bas_i, orb_i, 1) * coef (bas_j, orb_j, 1) * basis_ovlp (bas_i, bas_j)
         enddo ! bas_j
       enddo ! bas_i

       orb_vir_ovlp (orb_vir_in_wf_j, orb_vir_in_wf_i) = orb_vir_ovlp (orb_vir_in_wf_i, orb_vir_in_wf_j)

    enddo ! orb_vir_in_wf_j
  enddo ! orb_vir_in_wf_i

  end subroutine orb_vir_ovlp_bld

! ==============================================================================
  subroutine orb_cls_ovlp_inv_bld
! ------------------------------------------------------------------------------
! Description   : inverse of overlap matrix of closed orbitals
!
! Created       : J. Toulouse, 23 July 2007
! ------------------------------------------------------------------------------
  implicit none

! local
  real(dp) :: threshold

! header
  if (header_exe) then

   call object_create ('orb_cls_ovlp_inv')

   call object_needed ('orb_cls_ovlp')
   call object_needed ('orb_cls_in_wf_nb')

   return

  endif

! allocations
  call object_alloc ('orb_cls_ovlp_inv', orb_cls_ovlp_inv, orb_cls_in_wf_nb, orb_cls_in_wf_nb)

  threshold = 1.d-10
  call inverse_by_svd (orb_cls_ovlp, orb_cls_ovlp_inv, orb_cls_in_wf_nb, threshold)

  end subroutine orb_cls_ovlp_inv_bld

! ==============================================================================
  subroutine orb_occ_ovlp_inv_bld
! ------------------------------------------------------------------------------
! Description   : inverse of overlap matrix of occupied orbitals
!
! Created       : J. Toulouse, 22 May 2007
! ------------------------------------------------------------------------------
  implicit none

! local
  real(dp) :: threshold

! header
  if (header_exe) then

   call object_create ('orb_occ_ovlp_inv')

   call object_needed ('orb_occ_ovlp')
   call object_needed ('orb_occ_in_wf_nb')

   return

  endif

! allocations
  call object_alloc ('orb_occ_ovlp_inv', orb_occ_ovlp_inv, orb_occ_in_wf_nb, orb_occ_in_wf_nb)

  threshold = 1.d-10
  call inverse_by_svd (orb_occ_ovlp, orb_occ_ovlp_inv, orb_occ_in_wf_nb, threshold)

  end subroutine orb_occ_ovlp_inv_bld

! ==============================================================================
  subroutine orb_cls_ovlp_eig_bld
! ------------------------------------------------------------------------------
! Description   : Eigensystem of closed orbital overlap matrix
!
! Created       : J. Toulouse, 01 Jun 2007
! ------------------------------------------------------------------------------
  implicit none

! local
  integer orb_i, orb_j, lin_dep_nb
  real(dp) lin_dep_thres

! header
  if (header_exe) then

   call object_create ('orb_cls_ovlp_eigvec')
   call object_create ('orb_cls_ovlp_eigval')

   call object_needed ('orb_cls_in_wf_nb')
   call object_needed ('orb_cls_ovlp')

   return

  endif

! begin

! allocation
  call object_alloc ('orb_cls_ovlp_eigvec', orb_cls_ovlp_eigvec, orb_cls_in_wf_nb, orb_cls_in_wf_nb)
  call object_alloc ('orb_cls_ovlp_eigval', orb_cls_ovlp_eigval, orb_cls_in_wf_nb)

  if (orb_cls_in_wf_nb == 0) return

! diagonalization
  call eigensystem (orb_cls_ovlp, orb_cls_ovlp_eigvec, orb_cls_ovlp_eigval, orb_cls_in_wf_nb)

!  write(6,'(a)') 'Eigenvalues of overlap matrix of closed orbitals:'
!  do orb_i = 1, orb_cls_in_wf_nb
!     write(6,'(a,i3,a,e)') 'eigenvalue # ',orb_i,' : ', orb_cls_ovlp_eigval (orb_i)
!  enddo ! bas_i
!  write(6,'(a)') 'Eigenvectors:'
!  do orb_i = 1, orb_cls_in_wf_nb
!     write(6,'(a,i3,a,100f12.6)') 'eigenvector # ', orb_i,' :', (orb_cls_ovlp_eigvec (orb_i, orb_j), orb_j = 1, orb_cls_in_wf_nb)
!  enddo ! bas_i

! check linear dependancies
  lin_dep_thres = 1.d-12
  lin_dep_nb = 0
  do orb_i = 1, orb_cls_in_wf_nb
     if (orb_cls_ovlp_eigval (orb_i) < lin_dep_thres) then
       lin_dep_nb = lin_dep_nb + 1
     endif
  enddo ! orb_i
  if (lin_dep_nb > 0) then
   write(6,'(a,i3,a,e,a)') 'Warning: there are ',lin_dep_nb,' eigenvalues < ',lin_dep_thres,' in the overlap of the closed orbitals.'
  endif

  end subroutine orb_cls_ovlp_eig_bld

! ==============================================================================
  subroutine orb_act_ovlp_eig_bld
! ------------------------------------------------------------------------------
! Description   : Eigensystem of active orbital overlap matrix
!
! Created       : J. Toulouse, 23 Jul 2007
! ------------------------------------------------------------------------------
  implicit none

! local
  integer orb_i, orb_j, lin_dep_nb
  real(dp) lin_dep_thres

! header
  if (header_exe) then

   call object_create ('orb_act_ovlp_eigvec')
   call object_create ('orb_act_ovlp_eigval')

   call object_needed ('orb_act_in_wf_nb')
   call object_needed ('orb_act_ovlp')

   return

  endif

! begin

! allocation
  call object_alloc ('orb_act_ovlp_eigvec', orb_act_ovlp_eigvec, orb_act_in_wf_nb, orb_act_in_wf_nb)
  call object_alloc ('orb_act_ovlp_eigval', orb_act_ovlp_eigval, orb_act_in_wf_nb)

  if (orb_act_in_wf_nb == 0) return

! diagonalization
  call eigensystem (orb_act_ovlp, orb_act_ovlp_eigvec, orb_act_ovlp_eigval, orb_act_in_wf_nb)

!  write(6,'(a)') 'Eigenvalues of overlap matrix of active orbitals:'
!  do orb_i = 1, orb_act_in_wf_nb
!     write(6,'(a,i3,a,e)') 'eigenvalue # ',orb_i,' : ', orb_act_ovlp_eigval (orb_i)
!  enddo ! bas_i
!  write(6,'(a)') 'Eigenvectors:'
!  do orb_i = 1, orb_act_in_wf_nb
!     write(6,'(a,i3,a,100f12.6)') 'eigenvector # ', orb_i,' :', (orb_act_ovlp_eigvec (orb_i, orb_j), orb_j = 1, orb_act_in_wf_nb)
!  enddo ! bas_i

! check linear dependancies
  lin_dep_thres = 1.d-12
  lin_dep_nb = 0
  do orb_i = 1, orb_act_in_wf_nb
     if (orb_act_ovlp_eigval (orb_i) < lin_dep_thres) then
       lin_dep_nb = lin_dep_nb + 1
     endif
  enddo ! orb_i
  if (lin_dep_nb > 0) then
   write(6,'(a,i3,a,e,a)') 'Warning: there are ',lin_dep_nb,' eigenvalues < ',lin_dep_thres,' in the overlap of the active orbitals.'
  endif

  end subroutine orb_act_ovlp_eig_bld

! ==============================================================================
  subroutine orb_vir_ovlp_eig_bld
! ------------------------------------------------------------------------------
! Description   : Eigensystem of closed orbital overlap matrix
!
! Created       : J. Toulouse, 01 Jun 2007
! ------------------------------------------------------------------------------
  implicit none

! local
  integer orb_i, orb_j, lin_dep_nb
  real(dp) lin_dep_thres

! header
  if (header_exe) then

   call object_create ('orb_vir_ovlp_eigvec')
   call object_create ('orb_vir_ovlp_eigval')

   call object_needed ('orb_vir_in_wf_nb')
   call object_needed ('orb_vir_ovlp')

   return

  endif

! begin

! allocation
  call object_alloc ('orb_vir_ovlp_eigvec', orb_vir_ovlp_eigvec, orb_vir_in_wf_nb, orb_vir_in_wf_nb)
  call object_alloc ('orb_vir_ovlp_eigval', orb_vir_ovlp_eigval, orb_vir_in_wf_nb)

  if (orb_vir_in_wf_nb == 0) return

! diagonalization
  call eigensystem (orb_vir_ovlp, orb_vir_ovlp_eigvec, orb_vir_ovlp_eigval, orb_vir_in_wf_nb)

!  write(6,'(a)') 'Eigenvalues of overlap matrix of virtual orbitals:'
!  do orb_i = 1, orb_vir_in_wf_nb
!     write(6,'(a,i3,a,e)') 'eigenvalue # ',orb_i,' : ', orb_vir_ovlp_eigval (orb_i)
!  enddo ! bas_i
!  write(6,'(a)') 'Eigenvectors:'
!  do orb_i = 1, orb_vir_in_wf_nb
!     write(6,'(a,i3,a,100f12.6)') 'eigenvector # ', orb_i,' :', (orb_vir_ovlp_eigvec (orb_i, orb_j), orb_j = 1, orb_vir_in_wf_nb)
!  enddo ! bas_i

! check linear dependancies
  lin_dep_thres = 1.d-12
  lin_dep_nb = 0
  do orb_i = 1, orb_vir_in_wf_nb
     if (orb_vir_ovlp_eigval (orb_i) < lin_dep_thres) then
       lin_dep_nb = lin_dep_nb + 1
     endif
  enddo ! orb_i
  if (lin_dep_nb > 0) then
   write(6,'(a,i3,a,e,a)') 'Warning: there are ',lin_dep_nb,' eigenvalues < ',lin_dep_thres,' in the overlap of the virtual orbitals.'
  endif

  end subroutine orb_vir_ovlp_eig_bld

! ==============================================================================
  subroutine orb_cls_ovlp_m12_bld
! ------------------------------------------------------------------------------
! Description   : closed orbital overlap matrix to the power -1/2 for symmetric orthonormalization
!
! Created       : J. Toulouse, 01 Jan 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, orb_j, orb_k

! header
  if (header_exe) then

   call object_create ('orb_cls_ovlp_m12')

   call object_needed ('orb_cls_in_wf_nb')
   call object_needed ('orb_cls_ovlp_eigvec')
   call object_needed ('orb_cls_ovlp_eigval')

   return

  endif

! begin

! allocation
  call object_alloc ('orb_cls_ovlp_m12', orb_cls_ovlp_m12, orb_cls_in_wf_nb, orb_cls_in_wf_nb)

  do orb_i = 1, orb_cls_in_wf_nb
   do orb_j = 1, orb_cls_in_wf_nb
     orb_cls_ovlp_m12 (orb_i, orb_j) = 0.d0
     do orb_k = 1, orb_cls_in_wf_nb
      orb_cls_ovlp_m12 (orb_i, orb_j) = orb_cls_ovlp_m12 (orb_i, orb_j) + orb_cls_ovlp_eigvec (orb_i, orb_k) * (1.d0/dsqrt(orb_cls_ovlp_eigval (orb_k))) * orb_cls_ovlp_eigvec (orb_j, orb_k)
     enddo ! orb_k
   enddo ! orb_j
  enddo ! orb_i

  end subroutine orb_cls_ovlp_m12_bld

! ==============================================================================
  subroutine orb_act_ovlp_m12_bld
! ------------------------------------------------------------------------------
! Description   : active orbital overlap matrix to the power -1/2 for symmetric orthonormalization
!
! Created       : J. Toulouse, 23 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, orb_j, orb_k

! header
  if (header_exe) then

   call object_create ('orb_act_ovlp_m12')

   call object_needed ('orb_act_in_wf_nb')
   call object_needed ('orb_act_ovlp_eigvec')
   call object_needed ('orb_act_ovlp_eigval')

   return

  endif

! begin

! allocation
  call object_alloc ('orb_act_ovlp_m12', orb_act_ovlp_m12, orb_act_in_wf_nb, orb_act_in_wf_nb)

  do orb_i = 1, orb_act_in_wf_nb
   do orb_j = 1, orb_act_in_wf_nb
     orb_act_ovlp_m12 (orb_i, orb_j) = 0.d0
     do orb_k = 1, orb_act_in_wf_nb
      orb_act_ovlp_m12 (orb_i, orb_j) = orb_act_ovlp_m12 (orb_i, orb_j) + orb_act_ovlp_eigvec (orb_i, orb_k) * (1.d0/dsqrt(orb_act_ovlp_eigval (orb_k))) * orb_act_ovlp_eigvec (orb_j, orb_k)
     enddo ! orb_k
   enddo ! orb_j
  enddo ! orb_i

  end subroutine orb_act_ovlp_m12_bld

! ==============================================================================
  subroutine orb_vir_ovlp_m12_bld
! ------------------------------------------------------------------------------
! Description   : virtual orbital overlap matrix to the power -1/2 for symmetric orthonormalization
!
! Created       : J. Toulouse, 01 Jan 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, orb_j, orb_k

! header
  if (header_exe) then

   call object_create ('orb_vir_ovlp_m12')

   call object_needed ('orb_vir_in_wf_nb')
   call object_needed ('orb_vir_ovlp_eigvec')
   call object_needed ('orb_vir_ovlp_eigval')

   return

  endif

! begin

! allocation
  call object_alloc ('orb_vir_ovlp_m12', orb_vir_ovlp_m12, orb_vir_in_wf_nb, orb_vir_in_wf_nb)

  do orb_i = 1, orb_vir_in_wf_nb
   do orb_j = 1, orb_vir_in_wf_nb
     orb_vir_ovlp_m12 (orb_i, orb_j) = 0.d0
     do orb_k = 1, orb_vir_in_wf_nb
      orb_vir_ovlp_m12 (orb_i, orb_j) = orb_vir_ovlp_m12 (orb_i, orb_j) + orb_vir_ovlp_eigvec (orb_i, orb_k) * (1.d0/dsqrt(orb_vir_ovlp_eigval (orb_k))) * orb_vir_ovlp_eigvec (orb_j, orb_k)
     enddo ! orb_k
   enddo ! orb_j
  enddo ! orb_i

  end subroutine orb_vir_ovlp_m12_bld

! ==============================================================================
  subroutine ortho_orb_cls
! ------------------------------------------------------------------------------
! Description   : orthonormalize closed orbitals
!
! Created       : J. Toulouse, 01 Jan 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_cls_in_wf_i, orb_cls_i
  integer orb_cls_in_wf_j, orb_cls_j
  real(dp), allocatable :: coef_new (:,:)

! begin
  call object_provide ('orb_cls_in_wf_nb')
  call object_provide ('orb_cls_in_wf_lab')
  call object_provide ('orb_cls_ovlp_m12')
  call object_provide ('nbasis')
  call object_provide ('coef')

  call alloc ('coef_new', coef_new, nbasis, orb_cls_in_wf_nb)
  coef_new (:,:) = 0.d0

  do orb_cls_in_wf_i = 1, orb_cls_in_wf_nb
     orb_cls_i = orb_cls_in_wf_lab (orb_cls_in_wf_i)
     do orb_cls_in_wf_j = 1, orb_cls_in_wf_nb
        orb_cls_j = orb_cls_in_wf_lab (orb_cls_in_wf_j)
         coef_new (1:nbasis, orb_cls_in_wf_i) = coef_new (1:nbasis, orb_cls_in_wf_i) + orb_cls_ovlp_m12 (orb_cls_in_wf_i, orb_cls_in_wf_j) * coef (1:nbasis, orb_cls_j, 1)
     enddo ! orb_cls_in_wf_j
  enddo ! orb_cls_in_wf_i

  do orb_cls_in_wf_i = 1, orb_cls_in_wf_nb
     orb_cls_i = orb_cls_in_wf_lab (orb_cls_in_wf_i)
     coef (1:nbasis, orb_cls_i, 1) = coef_new (1:nbasis, orb_cls_in_wf_i)
  enddo ! orb_cls_in_wf_i

  call object_modified ('coef')

  end subroutine ortho_orb_cls

! ==============================================================================
  subroutine ortho_orb_act_to_orb_cls
! ------------------------------------------------------------------------------
! Description   : orthogonalize active (or open) orbitals to closed orbitals
!
! Created       : J. Toulouse, 23 July 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_act_in_wf_i, orb_act_i
  integer orb_cls_in_wf_i, orb_cls_i
  integer orb_cls_in_wf_j, orb_cls_j

! begin
  call object_provide ('orb_cls_in_wf_nb')
  call object_provide ('orb_cls_in_wf_lab')
  call object_provide ('orb_act_in_wf_nb')
  call object_provide ('orb_act_in_wf_lab')
  call object_provide ('orb_cls_ovlp_inv')
  call object_provide ('orb_ovlp')
  call object_provide ('nbasis')
  call object_provide ('coef')

  do orb_act_in_wf_i = 1, orb_act_in_wf_nb
     orb_act_i = orb_act_in_wf_lab (orb_act_in_wf_i)
     do orb_cls_in_wf_i = 1, orb_cls_in_wf_nb
        orb_cls_i = orb_cls_in_wf_lab (orb_cls_in_wf_i)
        do orb_cls_in_wf_j = 1, orb_cls_in_wf_nb
           orb_cls_j = orb_cls_in_wf_lab (orb_cls_in_wf_j)
            coef (1:nbasis, orb_act_i, 1) = coef (1:nbasis, orb_act_i, 1) - orb_cls_ovlp_inv (orb_cls_in_wf_i, orb_cls_in_wf_j) * orb_ovlp (orb_cls_j, orb_act_i) * coef (1:nbasis, orb_cls_i, 1)
        enddo ! orb_cls_in_wf_j
     enddo ! orb_cls_in_wf_i
  enddo ! orb_act_in_wf_i

  call object_modified ('coef')

  end subroutine ortho_orb_act_to_orb_cls

! ==============================================================================
  subroutine ortho_orb_act
! ------------------------------------------------------------------------------
! Description   : orthonormalize active (or open) orbitals
! Description   : This is energy invariant only for single-determinant or CASSCF wave function
!
! Created       : J. Toulouse, 23 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_act_in_wf_i, orb_act_i
  integer orb_act_in_wf_j, orb_act_j
  real(dp), allocatable :: coef_new (:,:)

! begin
  call object_provide ('orb_act_in_wf_nb')
  call object_provide ('orb_act_in_wf_lab')
  call object_provide ('orb_act_ovlp_m12')
  call object_provide ('nbasis')
  call object_provide ('coef')

  call alloc ('coef_new', coef_new, nbasis, orb_act_in_wf_nb)
  coef_new (:,:) = 0.d0

  do orb_act_in_wf_i = 1, orb_act_in_wf_nb
     orb_act_i = orb_act_in_wf_lab (orb_act_in_wf_i)
     do orb_act_in_wf_j = 1, orb_act_in_wf_nb
        orb_act_j = orb_act_in_wf_lab (orb_act_in_wf_j)
         coef_new (1:nbasis, orb_act_in_wf_i) = coef_new (1:nbasis, orb_act_in_wf_i) + orb_act_ovlp_m12 (orb_act_in_wf_i, orb_act_in_wf_j) * coef (1:nbasis, orb_act_j, 1)
     enddo ! orb_act_in_wf_j
  enddo ! orb_act_in_wf_i

  do orb_act_in_wf_i = 1, orb_act_in_wf_nb
     orb_act_i = orb_act_in_wf_lab (orb_act_in_wf_i)
     coef (1:nbasis, orb_act_i, 1) = coef_new (1:nbasis, orb_act_in_wf_i)
  enddo ! orb_act_in_wf_i

  call object_modified ('coef')

  end subroutine ortho_orb_act

! ==============================================================================
  subroutine ortho_orb_vir_to_orb_occ
! ------------------------------------------------------------------------------
! Description   : orthogonalize virtual orbitals to occupied orbitals
!
! Created       : J. Toulouse, 22 May 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_vir_in_wf_i, orb_vir_i
  integer orb_occ_in_wf_i, orb_occ_i
  integer orb_occ_in_wf_j, orb_occ_j

! begin
  call object_provide ('orb_occ_in_wf_nb')
  call object_provide ('orb_occ_in_wf_lab')
  call object_provide ('orb_vir_in_wf_nb')
  call object_provide ('orb_vir_in_wf_lab')
  call object_provide ('orb_occ_ovlp_inv')
  call object_provide ('orb_ovlp')
  call object_provide ('nbasis')
  call object_provide ('coef')

  do orb_vir_in_wf_i = 1, orb_vir_in_wf_nb
     orb_vir_i = orb_vir_in_wf_lab (orb_vir_in_wf_i)
     do orb_occ_in_wf_i = 1, orb_occ_in_wf_nb
        orb_occ_i = orb_occ_in_wf_lab (orb_occ_in_wf_i)
        do orb_occ_in_wf_j = 1, orb_occ_in_wf_nb
           orb_occ_j = orb_occ_in_wf_lab (orb_occ_in_wf_j)
            coef (1:nbasis, orb_vir_i, 1) = coef (1:nbasis, orb_vir_i, 1) - orb_occ_ovlp_inv (orb_occ_in_wf_i, orb_occ_in_wf_j) * orb_ovlp (orb_occ_j, orb_vir_i) * coef (1:nbasis, orb_occ_i, 1)
        enddo ! orb_occ_in_wf_j
     enddo ! orb_occ_in_wf_i
  enddo ! orb_vir_in_wf_i

  call object_modified ('coef')

  end subroutine ortho_orb_vir_to_orb_occ

! ==============================================================================
  subroutine ortho_orb_vir
! ------------------------------------------------------------------------------
! Description   : orthonormalize virtual orbitals
!
! Created       : J. Toulouse, 01 Jan 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_vir_in_wf_i, orb_vir_i
  integer orb_vir_in_wf_j, orb_vir_j
  real(dp), allocatable :: coef_new (:,:)

! begin
  call object_provide ('orb_vir_in_wf_nb')
  call object_provide ('orb_vir_in_wf_lab')
  call object_provide ('orb_vir_ovlp_m12')
  call object_provide ('nbasis')
  call object_provide ('coef')

  call alloc ('coef_new', coef_new, nbasis, orb_vir_in_wf_nb)
  coef_new (:,:) = 0.d0

  do orb_vir_in_wf_i = 1, orb_vir_in_wf_nb
     orb_vir_i = orb_vir_in_wf_lab (orb_vir_in_wf_i)
     do orb_vir_in_wf_j = 1, orb_vir_in_wf_nb
        orb_vir_j = orb_vir_in_wf_lab (orb_vir_in_wf_j)
         coef_new (1:nbasis, orb_vir_in_wf_i) = coef_new (1:nbasis, orb_vir_in_wf_i) + orb_vir_ovlp_m12 (orb_vir_in_wf_i, orb_vir_in_wf_j) * coef (1:nbasis, orb_vir_j, 1)
     enddo ! orb_vir_in_wf_j
  enddo ! orb_vir_in_wf_i

  do orb_vir_in_wf_i = 1, orb_vir_in_wf_nb
     orb_vir_i = orb_vir_in_wf_lab (orb_vir_in_wf_i)
     coef (1:nbasis, orb_vir_i, 1) = coef_new (1:nbasis, orb_vir_in_wf_i)
  enddo ! orb_vir_in_wf_i

  call object_modified ('coef')

  end subroutine ortho_orb_vir

! ==============================================================================
  subroutine ortho_orb
! ------------------------------------------------------------------------------
! Description   : orthonormalize orbitals
! Description   : This seems to break the optimization when optimizing the exponents?
!
! Created       : J. Toulouse, 23 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! orthonormalize orbitals
  write (6,'(a)') 'Orthonormalizing orbitals.'
  call ortho_orb_cls
  call ortho_orb_act_to_orb_cls
  if (ndet == 1 .or. l_casscf) then
   call ortho_orb_act
  endif
  call ortho_orb_vir_to_orb_occ
  call ortho_orb_vir

  call coef_orb_on_norm_basis_from_coef (1)
  if (trim(basis_functions_varied) == 'unnormalized') then
   call coef_orb_on_ortho_basis_from_coef (1)
  endif

  end subroutine ortho_orb

! ==============================================================================
  subroutine coef_orb_on_norm_basis_from_coef (iwf_from)
! ------------------------------------------------------------------------------
! Description   : compute orbital coefficients on nornalized basis
! Description   : from orbital coefficients on unnormalized basis
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  integer, intent(in) :: iwf_from

! local
  integer bas_i

! begin

! objects needed
  call object_provide ('nbasis')
  call object_provide ('orb_tot_nb')
  call object_provide ('norm_basis')
  call object_provide ('coef')

! allocations
  call object_alloc ('coef_orb_on_norm_basis', coef_orb_on_norm_basis, nbasis, orb_tot_nb, MFORCE)

  do bas_i = 1, nbasis
    coef_orb_on_norm_basis (bas_i, 1:orb_tot_nb, iwf) = coef (bas_i, 1:orb_tot_nb, iwf_from) / norm_basis (bas_i)
  enddo

  call object_modified ('coef_orb_on_norm_basis')

  end subroutine coef_orb_on_norm_basis_from_coef

! ==============================================================================
  subroutine coef_from_coef_orb_on_norm_basis (iwf_from)
! ------------------------------------------------------------------------------
! Description   : compute orbital coefficients on unnormalized basis
! Description   : from orbital coefficients on nornalized basis
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  integer, intent(in) :: iwf_from

! local
  integer bas_i

! begin

! object_needed
  call object_provide ('nbasis')
  call object_provide ('orb_tot_nb')
  call object_provide ('norm_basis')
  call object_provide ('coef_orb_on_norm_basis')

  do bas_i = 1, nbasis
    coef (bas_i, 1:orb_tot_nb, iwf) = coef_orb_on_norm_basis (bas_i, 1:orb_tot_nb, iwf_from) * norm_basis (bas_i)
  enddo

  call object_modified ('coef')

  end subroutine coef_from_coef_orb_on_norm_basis

! ==============================================================================
  subroutine coef_orb_on_ortho_basis_from_coef (iwf_from)
! ------------------------------------------------------------------------------
! Description   : compute orbital coefficients on orthornalized basis
! Description   : from orbital coefficients on unnormalized basis
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  integer, intent(in) :: iwf_from

! local
  character(len=max_string_len_rout), save :: lhere = 'coef_orb_on_ortho_basis_from_coef'
  integer bas_i, bas_k, orb_i

! begin
# if defined (DEBUG)
   call routine_enter (lhere)
# endif

! allocations
  call object_alloc ('coef_orb_on_ortho_basis', coef_orb_on_ortho_basis, nbasis, orb_tot_nb, MFORCE)

  call object_provide ('orb_tot_nb')
  call object_provide ('nbasis')
  call object_provide ('coef')
  call object_provide ('basis_ovlp_12')

  coef_orb_on_ortho_basis (:,:,:) = 0.d0

  do orb_i = 1, orb_tot_nb
    do bas_k = 1, nbasis
      do bas_i = 1, nbasis
       coef_orb_on_ortho_basis (bas_k, orb_i, iwf) = coef_orb_on_ortho_basis (bas_k, orb_i, iwf) + coef (bas_i, orb_i, iwf_from) * basis_ovlp_12 (bas_i, bas_k)
      enddo ! bas_k
    enddo ! bas_i
  enddo ! orb_i

  call object_modified ('coef_orb_on_ortho_basis')

# if defined (DEBUG)
   call routine_exit (lhere)
# endif

  end subroutine coef_orb_on_ortho_basis_from_coef

! ==============================================================================
  subroutine coef_from_coef_orb_on_ortho_basis (iwf_from)
! ------------------------------------------------------------------------------
! Description   : compute orbital coefficients on unnormalized basis
! Description   : from orbital coefficients on orthornalized basis
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  integer, intent(in) :: iwf_from

! local
  character(len=max_string_len_rout), save :: lhere = 'coef_from_coef_orb_on_ortho_basis'
  integer bas_i, bas_k, orb_i

! begin
# if defined (DEBUG)
   call routine_enter (lhere)
# endif

  call object_provide ('orb_tot_nb')
  call object_provide ('nbasis')
  call object_provide ('coef_orb_on_ortho_basis')
  call object_provide ('basis_ovlp_m12')

  do orb_i = 1, orb_tot_nb
    do bas_k = 1, nbasis
      coef (bas_k, orb_i, iwf) = 0.d0
      do bas_i = 1, nbasis
       coef (bas_k, orb_i, iwf) = coef (bas_k, orb_i, iwf) + coef_orb_on_ortho_basis (bas_i, orb_i, iwf_from) * basis_ovlp_m12 (bas_i, bas_k)
      enddo ! bas_k
    enddo ! bas_i
  enddo ! orb_i

  call object_modified ('coef')

# if defined (DEBUG)
   call routine_exit (lhere)
# endif

  end subroutine coef_from_coef_orb_on_ortho_basis

!===========================================================================
  subroutine orb_opt_lab_rd
!---------------------------------------------------------------------------
! Description : read orbital space to be optimized
!
! Created     : J. Toulouse, 25 Oct 2005
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'orb_opt_lab_rd'
  integer orb_i, orb_j

! begin
  call get_next_value_list ('orb_opt_lab', orb_opt_lab, orb_opt_nb)

! checking
  do orb_i = 1, orb_opt_nb
    do orb_j = 1, orb_opt_nb
      if (orb_i == orb_j) cycle
      if (orb_opt_lab (orb_i) == orb_opt_lab (orb_j) ) then
        write(6,'(2a,i3,a)') trim(lhere),': optimized orbital ', orb_opt_lab (orb_i), ' appears more than once'
        write(6,'(2a,i3,a,i3)') trim(lhere),': optimized orbital # ', orb_j, ' is ', orb_opt_lab (orb_j)
        write(6,'(2a,i3,a,i3)') trim(lhere),': optimized orbital # ', orb_j, ' is ', orb_opt_lab (orb_j)
        call die (lhere)
      endif
    enddo
  enddo

  if (orb_opt_nb > orb_tot_nb ) then
   write(6,'(2a,i3,a,i3)') trim(here),': orb_opt_nb=',orb_opt_nb,' > orb_tot_nb=',orb_tot_nb
   call die (lhere)
  endif

! label of last orbital in the optimization space
  orb_opt_last_lab = orb_opt_lab (orb_opt_nb)

  call object_associate ('orb_opt_nb', orb_opt_nb)
  call object_associate ('orb_opt_lab', orb_opt_lab, orb_opt_nb)
  call object_modified ('orb_opt_nb')
  call object_modified ('orb_opt_lab')
  call object_modified ('orb_opt_last_lab')

!  write(6,'(2a,i3)') trim(here),': number of orbitals in the optimization space =',orb_opt_nb
!  write(6,'(2a,100i3)') trim(here),': orbital labels =',orb_opt_lab

  end subroutine orb_opt_lab_rd

! ==============================================================================
  subroutine orb_opt_lab_bld
! ------------------------------------------------------------------------------
! Description   : build list of orbital labels to be optimized if not read
!
! Created       : J. Toulouse, 09 Jan 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i

! header
  if (header_exe) then

   call object_create ('orb_opt_nb')
   call object_create ('orb_opt_lab')
   call object_create ('orb_opt_last_lab')

   call object_needed ('orb_occ_in_wf_nb')

   return

  endif

! begin

! if orb_opt_nb not specified, take all the orbitals
   if (.not. object_valid ('orb_opt_nb')) then
    call object_provide ('orb_tot_nb')
    orb_opt_nb = orb_tot_nb
   endif

!  check orb_opt_nb
   call require ('orb_opt_nb > 0', orb_opt_nb > 0)
   if (orb_opt_nb <= orb_occ_in_wf_nb .and. (ndet == 1 .or. l_casscf)) then
     write(6,*) trim(here),': orb_opt_nb=',orb_opt_nb,' <= orb_occ_in_wf_nb=',orb_occ_in_wf_nb
     write(6,*) trim(here),': number of orbitals in optimization space must > number of occupied orbitals'
     call die (here)
   endif

!  orbital space for optimization
   call object_alloc ('orb_opt_lab', orb_opt_lab, orb_opt_nb)
   do orb_i = 1, orb_opt_nb
    orb_opt_lab (orb_i) = orb_i
    orb_opt_last_lab = orb_i
   enddo

!  write(6,'(2a,i3)') trim(here),': number of orbitals in the optimization space =',orb_opt_nb
!  write(6,'(2a,100i3)') trim(here),': orbital labels =',orb_opt_lab

 end subroutine orb_opt_lab_bld

!===========================================================================
  subroutine orb_ex_forbidden_rd
!---------------------------------------------------------------------------
! Description : read forbidden excitations for orbital optimization
!
! Created     : J. Toulouse, 25 Oct 2006
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'orb_ex_forbidden_rd'
  integer orb_i
  integer, allocatable :: orb_ex_forbidden_flat (:)
  integer orb_ex_forbidden_flat_nb, orb_ex_forbidden_nb

! begin

  call get_next_value_list ('orb_ex_forbidden_flat', orb_ex_forbidden_flat, orb_ex_forbidden_flat_nb)

  if (mod(orb_ex_forbidden_flat_nb,2) /= 0) then
    write(6,'(2a,i3,a)') trim(lhere),': the card "forbidden_excitations ... end" contains ', orb_ex_forbidden_flat_nb,' elements; it must be even number!'
    call die (lhere)
  endif

  orb_ex_forbidden_nb = orb_ex_forbidden_flat_nb/2

  call object_provide ('orb_tot_nb')

  do orb_i = 1, orb_ex_forbidden_flat_nb
    if (orb_ex_forbidden_flat (orb_i) < 1 .or. orb_ex_forbidden_flat (orb_i) > orb_tot_nb) then
      write(6,'(2a,i3,a,i3)') trim(lhere),': orbital label', orb_ex_forbidden_flat (orb_i),' is not between 1 and orb_tot_nb=',orb_tot_nb
      call die(lhere)
    endif
  enddo

  call object_alloc ('orb_ex_forbidden', orb_ex_forbidden, orb_tot_nb, orb_tot_nb)

  orb_ex_forbidden (:,:) = .false.

  do orb_i = 1, orb_ex_forbidden_nb
   orb_ex_forbidden (orb_ex_forbidden_flat(orb_i), orb_ex_forbidden_flat(orb_i+1)) = .true.
   write(6,'(2a,i3,a,i3,a)') trim(lhere),': orbital excitation ',orb_ex_forbidden_flat(orb_i),' -> ', orb_ex_forbidden_flat(orb_i+1),' is marked as forbidden'
  enddo


  call object_modified ('orb_ex_forbidden')

  end subroutine orb_ex_forbidden_rd

! ==============================================================================
  subroutine orb_ex_forbidden_bld
! ------------------------------------------------------------------------------
! Description   : build list of forbidden orbital excitation if not read
!
! Created       : J. Toulouse, 25 Oct 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('orb_ex_forbidden')

   call object_needed ('orb_tot_nb')

   return

  endif

! begin

  call object_alloc ('orb_ex_forbidden', orb_ex_forbidden, orb_tot_nb, orb_tot_nb)

  orb_ex_forbidden (:,:) = .false.

 end subroutine orb_ex_forbidden_bld

! ==============================================================================
  subroutine orb_optimized_bld
! ------------------------------------------------------------------------------
! Description   : orbital space to be optimized
!
! Created       : J. Toulouse, 25 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
!  include 'commons.h'

! local
  integer orb_i
  integer orb_opt_i, orb_opt
  integer orb_opt_occ_i, orb_opt_cls_i, orb_opt_opn_i, orb_opt_act_i, orb_opt_vir_i

! header
  if (header_exe) then

   call object_create ('orb_opt_occ_nb')
   call object_create ('orb_opt_occ_lab')
   call object_create ('orb_opt_cls_nb')
   call object_create ('orb_opt_cls_lab')
   call object_create ('orb_opt_opn_nb')
   call object_create ('orb_opt_opn_lab')
   call object_create ('orb_opt_act_nb')
   call object_create ('orb_opt_act_lab')
   call object_create ('orb_opt_vir_nb')
   call object_create ('orb_opt_vir_lab')

   call object_needed ('orb_opt_nb')
   call object_needed ('orb_opt_lab')
   call object_needed ('orb_occ_in_wf')
   call object_needed ('orb_cls_in_wf')
   call object_needed ('orb_opn_in_wf')
   call object_needed ('orb_act_in_wf')
   call object_needed ('orb_vir_in_wf')

   return

  endif

! begin
  orb_opt_occ_nb = 0
  orb_opt_cls_nb = 0
  orb_opt_opn_nb = 0
  orb_opt_act_nb = 0
  orb_opt_vir_nb = 0

  do orb_opt_i = 1, orb_opt_nb
   orb_opt = orb_opt_lab (orb_opt_i)

   if (orb_occ_in_wf (orb_opt)) then
    orb_opt_occ_nb = orb_opt_occ_nb + 1
   endif

   if (orb_cls_in_wf (orb_opt)) then
    orb_opt_cls_nb = orb_opt_cls_nb + 1
   endif

   if (orb_opn_in_wf (orb_opt)) then
    orb_opt_opn_nb = orb_opt_opn_nb + 1
   endif

   if (orb_act_in_wf (orb_opt)) then
    orb_opt_act_nb = orb_opt_act_nb + 1
   endif

   if (orb_vir_in_wf (orb_opt)) then
    orb_opt_vir_nb = orb_opt_vir_nb + 1
   endif

  enddo

! allocations
  call object_alloc ('orb_opt_occ_lab', orb_opt_occ_lab, orb_opt_occ_nb)
  call object_alloc ('orb_opt_cls_lab', orb_opt_cls_lab, orb_opt_cls_nb)
  call object_alloc ('orb_opt_opn_lab', orb_opt_opn_lab, orb_opt_opn_nb)
  call object_alloc ('orb_opt_act_lab', orb_opt_act_lab, orb_opt_act_nb)
  call object_alloc ('orb_opt_vir_lab', orb_opt_vir_lab, orb_opt_vir_nb)

  orb_opt_occ_i = 0
  orb_opt_cls_i = 0
  orb_opt_opn_i = 0
  orb_opt_act_i = 0
  orb_opt_vir_i = 0

  do orb_opt_i = 1, orb_opt_nb

   orb_opt = orb_opt_lab (orb_opt_i)

   if (orb_occ_in_wf (orb_opt)) then
    orb_opt_occ_i = orb_opt_occ_i + 1
    orb_opt_occ_lab (orb_opt_occ_i) = orb_opt
   endif

   if (orb_cls_in_wf (orb_opt)) then
    orb_opt_cls_i = orb_opt_cls_i + 1
    orb_opt_cls_lab (orb_opt_cls_i) = orb_opt
   endif

   if (orb_opn_in_wf (orb_opt)) then
    orb_opt_opn_i = orb_opt_opn_i + 1
    orb_opt_opn_lab (orb_opt_opn_i) = orb_opt
   endif

   if (orb_act_in_wf (orb_opt)) then
    orb_opt_act_i = orb_opt_act_i + 1
    orb_opt_act_lab (orb_opt_act_i) = orb_opt
   endif

   if (orb_vir_in_wf (orb_opt)) then
    orb_opt_vir_i = orb_opt_vir_i + 1
    orb_opt_vir_lab (orb_opt_vir_i) = orb_opt
   endif
  enddo

  if (orb_opt_occ_i /= orb_opt_occ_nb) then
    call die (here, 'orb_opt_occ_i='+orb_opt_occ_i+' /= orb_opt_occ_nb='+orb_opt_occ_nb)
  endif
  if (orb_opt_cls_i /= orb_opt_cls_nb) then
    call die (here, 'orb_opt_cls_i='+orb_opt_cls_i+' /= orb_opt_cls_nb='+orb_opt_cls_nb)
  endif
  if (orb_opt_opn_i /= orb_opt_opn_nb) then
    call die (here, 'orb_opt_opn_i='+orb_opt_opn_i+' /= orb_opt_opn_nb='+orb_opt_opn_nb)
  endif
  if (orb_opt_act_i /= orb_opt_act_nb) then
    call die (here, 'orb_opt_act_i='+orb_opt_act_i+' /= orb_opt_act_nb='+orb_opt_act_nb)
  endif
  if (orb_opt_vir_i /= orb_opt_vir_nb) then
    call die (here, 'orb_opt_vir_i='+orb_opt_vir_i+' /= orb_opt_vir_nb='+orb_opt_vir_nb)
  endif

  write(6,'(a,i3,a,1000i3)') 'There are ',orb_opt_nb,' orbitals in the optimization space of labels:', orb_opt_lab

 end subroutine orb_optimized_bld

! ==============================================================================
  subroutine is_orb_s_bld
! ------------------------------------------------------------------------------
! Description   : determine whether the orbitals are of the s-type
! Description   : same information that "lo" array
!
! Created       : J. Toulouse, 06 Jan 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, bas_i, l_i

! header
  if (header_exe) then

   call object_create ('is_orb_s')

   call object_needed ('nbasis')
   call object_needed ('orb_tot_nb')

   return

  endif

! begin

! preconditions
  if (inum_orb /= 0) then
   write(6,*) trim(here),': inum_orb=',inum_orb,' /= 0'
   write(6,*) trim(here),': implemented only for no numerical orbitals'
   call die (here)
  endif

! allocations
  call alloc ('is_orb_s', is_orb_s, orb_tot_nb)

  call object_provide ('coef')

  is_orb_s = .false.

  do orb_i = 1, orb_tot_nb
   do bas_i = 1, nbasis
     l_i   = l_bas (bas_i)

!    if current orbital has non-zero component on s-type basis function
     if (l_i == 0) then
       if (coef(bas_i,orb_i,1) /= 0.d0) then
         is_orb_s (orb_i) = .true.
         cycle
       endif
     endif

   enddo
  enddo

  write(6,*) 'is_orb_s=',is_orb_s

 end subroutine is_orb_s_bld

! ==============================================================================
  subroutine lo_bld
! ------------------------------------------------------------------------------
! Description   : determine whether the orbitals are of the s-type
! Description   : lo= 0 is s-type orbital, 1 otherwise
!
! Created       : J. Toulouse, 06 Jan 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, bas_i, l_i

! header
  if (header_exe) then

   call object_create ('lo')

   call object_needed ('nbasis')
   call object_needed ('orb_tot_nb')

   return

  endif

! begin

! preconditions
  if (inum_orb /= 0) then
   write(6,*) trim(here),': inum_orb=',inum_orb,' /= 0'
   write(6,*) trim(here),': implemented only for no numerical orbitals'
   call die (here)
  endif

  call object_provide ('coef')


  do orb_i = 1, orb_tot_nb
   lo (orb_i) = 1
   do bas_i = 1, nbasis
     l_i   = l_bas (bas_i)

!    if current orbital has non-zero component on s-type basis function
     if (l_i == 0) then
       if (coef(bas_i,orb_i,1) /= 0.d0) then
         lo (orb_i) = 0
         cycle
       endif
     endif

   enddo
  enddo

!  write(6,*) trim(here),': lo=',lo

 end subroutine lo_bld

! ==============================================================================
!  subroutine coef_cusp_bld
!! ------------------------------------------------------------------------------
!! Description   : coef of orbitals with e-N cusp condition imposed
!!
!! Created       : J. Toulouse, 06 Jan 2005
!! ------------------------------------------------------------------------------
!  implicit none
!  include 'commons.h'
!
!! local
!  character(len=max_string_len_rout), save :: here
!  integer bas_i
!  integer orb_i
!
!! header
!  if (header_exe) then
!
!   here = 'coef_cusp_bld'
!
!   call object_create ('coef_cusp')
!
!   call object_needed ('nbasis')
!   call object_needed ('orb_tot_nb')
!   call object_needed ('coef')
!   call object_needed ('orb_occ_in_wf_lab_nb')
!   call object_needed ('orb_occ_in_wf_lab')
!
!   return
!
!  endif
!
!! begin
!!  write(6,*) trim(here),': entering'
!
!! preconditions
!  if (ncent /= 1 ) then
!   write(6,*) trim(here),': ncent=',ncent,' /= 1'
!   write(6,*) trim(here),': implemented only for one center'
!   call die (here)
!  endif
!
!  if (inum_orb /= 0) then
!   write(6,*) trim(here),': inum_orb=',inum_orb,' /= 0'
!   write(6,*) trim(here),': not implemented for numerical orbitals'
!   call die (here)
!  endif
!
!  if (ijas /= 4) then
!   write(6,*) trim(here),': ijas=',ijas,' /= 4'
!   write(6,*) trim(here),': implemented only for jastrow 4'
!   call die (here)
!  endif
!
!! allocations
!  call object_alloc ('coef_cusp', coef_cusp, nbasis, orb_tot_nb)
!
!  aa1 = a4(1,1,1)
!
!  do orb_i = 1, orb_occ_in_wf_lab_nb
!
!   orb_occ = orb_occ_in_wf_lab (orb_i)
!
!   if (.not. is_orb_s (orb_occ)) cycle
!
!   do bas_i = 1, nbasis
!    do orb_i = 1, orb_tot_nb
!      coef_sav (bas_i, orb_i) = coef (bas_i, orb_i, 1)
!    enddo
!
!  enddo ! orb_i
!
!  end subroutine coef_sav_bld

! ==============================================================================
  subroutine cusp_en_orb
! ------------------------------------------------------------------------------
! Description   : check or impose e-N cusp conditions orbitals
!
! Created       : J. Toulouse, 06 Jan 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'cusp_en_orb'
  real(dp), allocatable :: diff (:)

! begin
!  if (icusp < 0) return

  call object_provide ('ncent')
  call object_provide ('orb_tot_nb')
  call object_provide ('nloc')
  call object_provide ('numr')

  call alloc('diff', diff, ncent*orb_tot_nb)

  if((nloc.eq.0. .or. nloc.eq.5) .and. numr.le.0) then
     if (l_cusp_en) then
       icusp = 1
       if(l_cusp_en_occ) then
        write(6,'(a)') 'imposing e-N cusp conditions on occupied orbitals:'
       else
        write(6,'(a)') 'imposing e-N cusp conditions on all orbitals:'
       endif
     else
        icusp = -1
        write(6,'(a)') 'checking e-N cusp conditions on orbitals:'
     endif
     call equiv_bas
     call cuspco(diff,1)
  endif

  call object_modified ('coef')

  end subroutine cusp_en_orb

! ==============================================================================
  subroutine write_orbitals_pw_real
! ------------------------------------------------------------------------------
! Description : write in a file orbitals on plane waves that have already been
! Description : converted to be real.
!
! Created     : J. Toulouse, 26 Mar 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'write_orbitals_pw_real'
  integer file_unit
  integer ikv, k, iband, jorb, igv

! begin
  call object_provide ('nkvec')
  call object_provide ('ngvec_orb')
  call object_provide ('nband')
  call object_provide ('rkvec')
  call object_provide ('ndim')
  call object_provide ('c_rp')
  call object_provide ('c_rm')
  call object_provide ('c_ip')
  call object_provide ('c_im')

! open file
  call open_file_or_die (file_orbitals_pw_out, file_unit)

! write in file
  write (file_unit,'(i3,i6,a)') nkvec, ngvec_orb, ' nkvec, ngvec_orb'
  do ikv = 1, nkvec
    write (file_unit,'(i3,i3,3f12.8,a)') ikv, nband(ikv), (rkvec (k,ikv), k = 1, ndim), ' ikvec, nband, rkvec'
    do iband = 1, nband(ikv)
      jorb = jorb + 1
      write (file_unit,'(i3,f12.7,a)') iband, 0.d0, ' band, eigenvalue'
      do igv = 1, ngvec_orb
        write (file_unit,' (4e16.8,a)') c_rp(igv,jorb),c_rm(igv,iband),c_ip(igv,jorb),c_im(igv,jorb), ' c_rp,c_rm,c_ip,c_im'
      enddo ! igv
    enddo ! iband
  enddo ! ikv

  close (file_unit)

  end subroutine write_orbitals_pw_real

end module orbitals_mod

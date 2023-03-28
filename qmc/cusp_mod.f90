module cusp_mod

  use all_tools_mod
  use basis_mod
  use orbitals_mod

! Declaration of global variables and default values
  logical                             :: l_check_cusp_en = .false.
  logical                             :: l_impose_cusp_en = .false.
  logical                             :: l_impose_cusp_en_occ = .false.
  logical                             :: l_impose_cusp_en_opt = .false.

  contains

!===========================================================================
  subroutine orb_cusp_menu
!---------------------------------------------------------------------------
! Description : menu for e-N cusp on orbitals
!
! Created     : J. Toulouse, 13 Jan 2006
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout),save :: lhere = 'orb_cusp_menu'
  integer i, it

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for cusp menu:'
   write(6,'(a)') ' cusp'
   write(6,'(a)') '  check_cusp_en = [bool] : check e-n cusp conditions? (default=false)'
   write(6,'(a)') '  impose_cusp_en = [bool] : impose e-n cusp conditions on all orbitals (default=false)'
   write(6,'(a)') '  impose_cusp_en_opt = [bool] : impose e-n cusp conditions during orbital optimization (default=false)'
   write(6,'(a)') '  impose_cusp_en_occ = [bool] : impose e-n cusp conditions on occupied orbitals only during orbital optimization (default=false)'
   write(6,'(a)') ' end'

  case ('check_cusp_en')
   call get_next_value (l_check_cusp_en)

  case ('impose_cusp_en')
   call get_next_value (l_impose_cusp_en)

  case ('impose_cusp_en_occ')
   call get_next_value (l_impose_cusp_en_occ)
   if (l_impose_cusp_en_occ) l_impose_cusp_en = .true.

  case ('impose_cusp_en_opt')
   call get_next_value (l_impose_cusp_en_opt)
   if (l_impose_cusp_en_opt) l_impose_cusp_en = .true.

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  if (l_impose_cusp_en_opt) then
   call die (lhere, 'option impose_cusp_en_opt=true not fully checked yet.')
  endif
  if (l_impose_cusp_en_occ) then
   call die (lhere, 'option impose_cusp_en_occ=true not fully checked yet.')
  endif

  if (l_check_cusp_en .or. l_impose_cusp_en) then

!   initialization for cusp
    call object_provide ('nbasis_ctype')
    call object_provide ('iwctype')
    call object_provide ('ncent')
    call alloc ('imnbas', imnbas, ncent)
    imnbas(1)=1
    do i=1,ncent-1
       it=iwctype(i)
       imnbas(i+1)=imnbas(i)+nbasis_ctype(it)
    enddo
    call ie

!   imposing or checking e-n cusp conditions
    call cusp_en_orb
  endif

  end subroutine orb_cusp_menu

! ==============================================================================
  subroutine cusp_en_orb
! ------------------------------------------------------------------------------
! Description   : check or impose e-n cusp conditions on orbital coefficients
!
! Created       : J. Toulouse, 06 Jan 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  real(dp), allocatable :: diff (:)
  integer bas_i, orb_i
! begin
!  if (icusp < 0) return

  call object_provide ('ncent')
  call object_provide ('orb_tot_nb')
  call object_provide ('nloc')
  call object_provide ('numr')
  call object_provide ('lo')
  call object_provide ('nbasis')
  call object_provide ('coef')
  call object_provide ('norb')
  call alloc('diff', diff, ncent*orb_tot_nb)

! With the current implementation checking or imposing the e-n cusp will only work with a purely analytical basis,
! not a mixed analytic-numerical basis
! if((nloc.eq.0. .or. nloc.eq.5) .and. numr.le.0) then
  if((nloc.eq.0. .or. nloc.eq.5) .and. minval(zex(:,1)).ne.0.d0) then

    call coef_orb_on_norm_basis_from_coef(1)
    call object_provide ('coef_orb_on_norm_basis')
   
    if (l_check_cusp_en .and. .not. l_impose_cusp_en) then
      icusp = -1
      write(6,'(a)') 'checking e-n cusp conditions on orbitals:'
    endif

    if (l_impose_cusp_en) then
      write(6,'(a)') '   Orbitals before imposition of e-n cusp conditions:'
      do orb_i=1,norb
        write(6,'(100f10.6)') (coef_orb_on_norm_basis(bas_i,orb_i,1),bas_i=1,nbasis)
      enddo
      write(6,'(a)') '   ----------------------------------------'   
      do orb_i=1,norb
        write(6,'(100f10.6)') (coef(bas_i,orb_i,1),bas_i=1,nbasis)
      enddo
      write(6,'(a)') '   ----------------------------------------'

      icusp = 1
      if (l_impose_cusp_en_occ) then
        write(6,'(a)') 'imposing e-n cusp conditions on occupied orbitals:'
      else
        write(6,'(a)') 'imposing e-n cusp conditions on orbitals:'
      endif
    endif

    call equiv_bas
    call cuspco(diff,1)
    call object_modified ('coef')
    write(6,'(a)') '   Orbitals after imposition of e-n cusps conditions:'
    call coef_orb_on_norm_basis_from_coef (1)
    do orb_i=1,norb
      write(6,'(100f10.6)') (coef_orb_on_norm_basis(bas_i,orb_i,1),bas_i=1,nbasis)
    enddo
    write(6,'(a)') '   ----------------------------------------'
    do orb_i=1,norb
      write(6,'(100f10.6)') (coef(bas_i,orb_i,1),bas_i=1,nbasis)
    enddo
    write(6,'(a)') '   ----------------------------------------'
  endif
  
  end subroutine cusp_en_orb

end module cusp_mod

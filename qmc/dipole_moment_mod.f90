 module dipole_moment_mod

  use all_tools_mod
  use electrons_mod
  use nuclei_mod

! Declaration of global variables and default values
  real(dp), allocatable                       :: dipole_moment_origin (:)
  real(dp), allocatable                       :: dipole_moment_nucl (:)
  real(dp), allocatable                       :: dipole_moment (:,:)
  real(dp), allocatable                       :: dipole_moment_av (:)
  real(dp), allocatable                       :: dipole_moment_av_err (:)
  real(dp), parameter                         :: atomic_unit_to_debye = 1.d0/0.393456d0

  contains

!===========================================================================
  subroutine dipole_moment_menu
!---------------------------------------------------------------------------
! Description : menu for calculation of dipole moment
!
! Created     : J. Toulouse, 02 Feb 2008
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'dipole_moment_menu'
  integer dipole_moment_origin_dim

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))
   case ('help')
    write(6,'(a)') 'HELP for menu dipole_moment'
    write(6,'(a)') ': dipole_moment'
    write(6,'(a)') ':  origin real  real  real end : origin with respect to which dipole moment is calculated (default = mass center)'
    write(6,'(a)') ': end'

   case ('origin')
    call get_next_value_list ('dipole_moment_origin', dipole_moment_origin, dipole_moment_origin_dim)
    call require ('dipole_moment_origin_dim == ndim', dipole_moment_origin_dim == ndim)

   case ('end')
    exit

   case default
    call die (lhere, 'unknown keyword >'+trim(word)+'<.')

  end select

  enddo ! end loop over menu lines

  call object_average_request ('dipole_moment_av')
  call object_error_request ('dipole_moment_av_err')

  end subroutine dipole_moment_menu

! ==============================================================================
  subroutine dipole_moment_origin_bld
! ------------------------------------------------------------------------------
! Description   : default origin for dipole moment
!
! Created       : J. Toulouse, 02 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_origin')

   call object_needed ('ndim')
   call object_needed ('mass_nucl_center')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_origin', dipole_moment_origin, ndim)
  dipole_moment_origin (:) = mass_nucl_center (:)

  end subroutine dipole_moment_origin_bld

! ==============================================================================
  subroutine dipole_moment_nucl_bld
! ------------------------------------------------------------------------------
! Description   : nuclear contribution to dipole moment
!
! Created       : J. Toulouse, 02 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, cent_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_nucl')

   call object_needed ('ndim')
   call object_needed ('znuc')
   call object_needed ('iwctype')
   call object_needed ('ncent')
   call object_needed ('cent')
   call object_needed ('dipole_moment_origin')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_nucl', dipole_moment_nucl, ndim)

  do dim_i = 1, ndim
   dipole_moment_nucl (dim_i) = 0.d0
   do cent_i = 1, ncent
    dipole_moment_nucl (dim_i) = dipole_moment_nucl (dim_i) + atomic_unit_to_debye * znuc(iwctype(cent_i)) * (cent (dim_i, cent_i) - dipole_moment_origin (dim_i))
   enddo ! cent_i
  enddo ! dim_i

  end subroutine dipole_moment_nucl_bld

! ==============================================================================
  subroutine dipole_moment_bld
! ------------------------------------------------------------------------------
! Description   : total dipole moment
!
! Created       : J. Toulouse, 02 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer elec_i, dim_i, walk_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment')
   call object_average_walk_define ('dipole_moment', 'dipole_moment_av')
   call object_error_define ('dipole_moment_av', 'dipole_moment_av_err')

   call object_needed ('nwalk')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('dipole_moment_nucl')
   call object_needed ('coord_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('dipole_moment', dipole_moment, ndim, nwalk)
  call object_alloc ('dipole_moment_av', dipole_moment_av, ndim)
  call object_alloc ('dipole_moment_av_err', dipole_moment_av_err, ndim)

  do walk_i = 1, nwalk
   do dim_i = 1, ndim
    dipole_moment (dim_i, walk_i) = dipole_moment_nucl (dim_i)
    do elec_i = 1, nelec
     dipole_moment (dim_i, walk_i) = dipole_moment (dim_i, walk_i) - atomic_unit_to_debye * (coord_elec_wlk (dim_i, elec_i, walk_i) - dipole_moment_origin (dim_i))
    enddo ! cent_i
   enddo ! dim_i
  enddo ! walk_i

  end subroutine dipole_moment_bld

end module dipole_moment_mod

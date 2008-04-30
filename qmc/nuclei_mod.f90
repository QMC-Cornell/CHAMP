module nuclei_mod

  use all_tools_mod

! Declaration of global variables and default values
  logical                   :: l_convert_from_angstrom_to_bohr = .false.
  real(dp), allocatable     :: dist_nn  (:,:)
  real(dp), allocatable     :: mass_nucl  (:)
  real(dp), allocatable     :: mass_nucl_center (:)
  real(dp)                  :: mass_nucl_total

  contains

!===========================================================================
  subroutine nuclei_menu
!---------------------------------------------------------------------------
! Description : menu for nuclei
!
! Created     : J. Toulouse, 30 Apr 2008
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'nuclei_menu'
  integer cent_i

! initialization
  l_convert_from_angstrom_to_bohr = .false.

  write(6,*)
  write(6,'(a)') 'Beginning of nuclei menu ---------------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') ' HELP for nuclei menu:'
   write(6,'(a)') '  nuclei'
   write(6,'(a)') '   convert_from_angstrom_to_bohr = [logical] : convert nuclear coordinates from Angstrom to Bohr units (default=false)'
   write(6,'(a)') '  end'
   write(6,*)

  case ('convert_from_angstrom_to_bohr')
   call get_next_value (l_convert_from_angstrom_to_bohr)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

! convert unit of nuclear coordinates
  if (l_convert_from_angstrom_to_bohr) then
   write(6,'(a)') ' Converting nuclear coordinates from Angstrom to Bohr:'
   call object_provide ('ncent')
   call object_provide ('ndim')
   call object_provide ('cent')
   cent(1:ndim,1:ncent) = cent(1:ndim,1:ncent) * angstrom_to_bohr
   call object_modified ('cent')
   do cent_i = 1,ncent
    write(6,'(a,i4,a,3f8.5)') ' center # ',cent_i,' : ',cent(1:ndim,cent_i)
   enddo
!  recalculate nuclear potential energy
   call object_provide ('znuc')
   call object_provide ('iwctype')
   call pot_nn(cent,znuc,iwctype,ncent,pecent)
   write(6,'(a,f14.7)') ' recalculting nuclear potential energy: pecent=',pecent
   call object_modified ('pecent')
  endif

  write(6,'(a)') 'End of nuclei menu ---------------------------------------------------------------------------------------'

  end subroutine nuclei_menu

! ==============================================================================
  subroutine mass_nucl_bld
! ------------------------------------------------------------------------------
! Description   : atomic mass of each nuclear center
!
! Created       : J. Toulouse, 02 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer cent_i
  integer znuc_int
  real(dp) :: znuc_to_mass_nucl (10) = (/1, 4, 7, 9, 11, 12, 14, 16, 19, 20/)

! begin

! header
  if (header_exe) then

   call object_create ('mass_nucl')

   call object_needed ('ncent')
   call object_needed ('iwctype')
   call object_needed ('znuc')

   return

  endif

! allocations
  call object_alloc ('mass_nucl', mass_nucl, ncent)

  do cent_i = 1, ncent
   znuc_int = int(znuc(iwctype(cent_i)))
   if (znuc_int >= 1 .and. znuc_int <= 10) then
    mass_nucl (cent_i) = znuc_to_mass_nucl (znuc_int)
   else
    call die (here, 'atomic mass associated to nuclear charge ='+ znuc_int+' is unknown.')
   endif
  enddo ! cent_i


  end subroutine mass_nucl_bld

! ==============================================================================
  subroutine mass_nucl_total_bld
! ------------------------------------------------------------------------------
! Description   : total nuclear mass
!
! Created       : J. Toulouse, 03 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer cent_i

! begin

! header
  if (header_exe) then

   call object_create ('mass_nucl_total')

   call object_needed ('ncent')
   call object_needed ('mass_nucl')

   return

  endif

! allocations
  mass_nucl_total = 0.d0

  do cent_i = 1, ncent
    mass_nucl_total = mass_nucl_total + mass_nucl (cent_i)
  enddo ! cent_i

  end subroutine mass_nucl_total_bld

! ==============================================================================
  subroutine mass_nucl_center_bld
! ------------------------------------------------------------------------------
! Description   : calculate nuclear center of mass
!
! Created       : J. Toulouse, 03 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, cent_i

! begin

! header
  if (header_exe) then

   call object_create ('mass_nucl_center')

   call object_needed ('ndim')
   call object_needed ('ncent')
   call object_needed ('cent')
   call object_needed ('mass_nucl_total')

   return

  endif

! allocations
  call object_alloc ('mass_nucl_center', mass_nucl_center, ndim)

  do dim_i = 1, ndim
   mass_nucl_center (dim_i) = 0.d0
   do cent_i = 1, ncent
    mass_nucl_center (dim_i) = mass_nucl_center (dim_i) + mass_nucl (cent_i) * cent (dim_i, cent_i)
   enddo ! cent_i
    mass_nucl_center (dim_i) = mass_nucl_center (dim_i) / mass_nucl_total
  enddo ! dim_i

!  write (6,*) "center of mass:", mass_nucl_center(:)

  end subroutine mass_nucl_center_bld

! ==============================================================================
  subroutine dist_nn_bld
! ------------------------------------------------------------------------------
! Description   : distance |Ri - Rj| between two nuclei
!
! Created       : J. Toulouse, 26 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer  cent_i, cent_j, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('dist_nn')

   call object_needed ('ncent')
   call object_needed ('ndim')
   call object_needed ('cent')

   return

  endif

! allocations
  call object_alloc ('dist_nn', dist_nn, ncent, ncent)
  dist_nn (:,:) = 0.d0

  do cent_i = 1, ncent
    do cent_j = cent_i+1, ncent
      do dim_i = 1, ndim
        dist_nn (cent_i, cent_j) = dist_nn (cent_i, cent_j) + (cent(dim_i,cent_i) - cent(dim_i,cent_j))**2
      enddo ! dim_i
      dist_nn (cent_i, cent_j) = dsqrt (dist_nn (cent_i, cent_j))
      dist_nn (cent_j, cent_i) = dist_nn (cent_i, cent_j)
    enddo ! cent_j
  enddo ! cent_i

  end subroutine dist_nn_bld

end module nuclei_mod

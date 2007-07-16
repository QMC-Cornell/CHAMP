module restart_mod

  use all_tools_mod

! Declaration of global variables and default values

  contains

! ==============================================================================
  subroutine dumper_dmc
! ------------------------------------------------------------------------------
! Description   : Choice of dmc dumper routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character (len=max_string_len_rout), save :: lhere = 'dumper_dmc'

! begin
# if defined (MPI)
    if (l_mode_dmc_mov1_mpi1) then
     call dumper_dmc_mov1_mpi1
    else
     call dumper_dmc_mov1_mpi2
    endif
# else
    if (l_mode_dmc_mov1) then
     call dumper_dmc_mov1
    else
     call dumper_dmc_movall
    endif
# endif

  end subroutine dumper_dmc

! ==============================================================================
  subroutine startr_dmc
! ------------------------------------------------------------------------------
! Description   : Choice of dmc startr routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character (len=max_string_len_rout), save :: lhere = 'startr_dmc'

! begin
# if defined (MPI)
    if (l_mode_dmc_mov1_mpi1) then
     call startr_dmc_mov1_mpi1
    else
     call startr_dmc_mov1_mpi2
    endif
# else
    if (l_mode_dmc_mov1) then
     call startr_dmc_mov1
    else
     call startr_dmc_movall
    endif
# endif

  end subroutine startr_dmc

end module restart_mod

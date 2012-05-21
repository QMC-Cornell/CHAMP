module nodes_mod

  use basic_tools_mod
  use constants_mod
  use strings_tools_mod


! Declaration of global variables and default values
  integer*8, external            :: address

  type type_node
   character(len=max_string_len) :: routine_name
   integer*8                     :: routine_address
   integer, allocatable          :: objects_needed_index(:)
   integer, allocatable          :: objects_create_index(:)
   logical                       :: valid
   logical                       :: debug = .false.
   logical                       :: entered = .false.
   integer                       :: calls_nb = 0
   real(dp)                      :: cpu_duration = 0.d0
  end type type_node

  integer                        :: nodes_nb = 0
  integer                        :: node_current_index = 0
  type (type_node), save         :: nodes (max_nodes_nb)

  logical                        :: header_exe = .false.
  logical                        :: l_trace = .false.

  integer, allocatable           :: nodes_indexes_sort (:)
  integer, allocatable           :: nodes_calls_nb_sort (:)
  real(dp), allocatable          :: nodes_cpu_duration_sort (:)

  contains

! ==============================================================================
  subroutine catalog_one_node (routine_name, routine, node_ind)
! ------------------------------------------------------------------------------
! Description   : catalog a new node
! Description   : warning: check the addresses in fortran and c
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent (in) :: routine_name
  external routine

! input/output
  integer, optional, intent(inout) :: node_ind

! local
  character(len=max_string_len_rout) :: lhere = 'catalog_one_node'
  integer node_i, nod_ind
  integer*8 rout_address

! begin

! make sure that node_ind = 0
  if (present(node_ind)) then
   if (node_ind /= 0) then
    write(6,*) trim(lhere), ': node_name=',trim(routine_name)
    write(6,*) trim(lhere), ': calling with node_ind=',node_ind, '/=0'
    call die (lhere)
   endif
  endif

! index of node
  nod_ind = node_index (routine_name)

! if routine already catalogued, die
  if (nod_ind /= 0) then
    call die (lhere, 'node >'+trim(routine_name)+'< catalogued more than once.')
  endif

! check maximun number of nodes
  if (nodes_nb == max_nodes_nb) then
   write(6,*) trim(lhere),': maximum number of nodes max_nodes_nb=',max_nodes_nb,' reached'
   call die (lhere)
  endif

! get routine's adress
  rout_address = address(routine)

! if a previously catalogued node has the same address, die
  do node_i = 1, nodes_nb
   if (rout_address == nodes(node_i)%routine_address) then
     write(6,'(6a)') trim(lhere),': nodes ', trim(routine_name),' and ', trim(nodes(node_i)%routine_name),' have the same address'
     write(6,'(6a)') trim(lhere),': check the declaration of these nodes in the file catalog_routines_mod.f90'
     call die (lhere)
   endif
  enddo ! node_i

! catalog node
  nodes_nb = nodes_nb + 1
  nod_ind = nodes_nb

  nodes(nod_ind)%routine_name = routine_name
  nodes(nod_ind)%routine_address = rout_address
  nodes(nod_ind)%valid = .false.

! save node index
  if (present(node_ind)) then
   node_ind = nod_ind
  endif

 end subroutine catalog_one_node

!===============================================================================
  subroutine execute_node_headers
! ------------------------------------------------------------------------------
! Description   : execute node headers
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  integer node_i

! begin
  header_exe = .true.

  do node_i = 1, nodes_nb
   node_current_index = node_i
   call exe_by_address_0 (nodes(node_i)%routine_address)
  enddo

  header_exe  = .false.
  node_current_index = 0

 end subroutine execute_node_headers

! ==============================================================================
  function node_index (node_name) result(result)
! ------------------------------------------------------------------------------
! Description   : return the index of a node, 0 if node not catalogued
!
! Created       : J. Toulouse, 15 Dec 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: node_name

! output
  integer                      :: result

! local
  integer nod_i
  integer length, half
  character(len=max_string_len) current_node_name

! begin

  length = len_trim(node_name)
  half   = max(1,length/2)

! check if node is catalogued
   do nod_i = 1, nodes_nb

    current_node_name = nodes(nod_i)%routine_name

!    if (length /= len_trim(current_node_name)) cycle  ! too slow
    if (node_name(1:1) /= current_node_name(1:1) ) cycle
    if (node_name(length:length) /= current_node_name(length:length) ) cycle
    if (node_name(half:half) /= current_node_name(half:half) ) cycle

    if (node_name == current_node_name ) then
      result = nod_i
      return
    endif

   enddo

! if object not found
  result = 0

  return

 end function node_index

! ==============================================================================
  subroutine node_add (node_name)
! ------------------------------------------------------------------------------
! Description   : add a new node to the catalog
!
! Created       : J. Toulouse, 18 Dec 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: node_name

! local
  character(len=max_string_len_rout), save :: lhere = 'node_add'

! begin

  if (nodes_nb == max_nodes_nb) then
   write(6,*) trim(lhere),': maximum number of nodes max_nodes_nb=',max_nodes_nb,' reached'
   call die (lhere)
  endif

  nodes_nb = nodes_nb + 1
  nodes(nodes_nb)%routine_name = node_name

 end subroutine node_add

!===============================================================================
  subroutine node_enter (node_name)
! ------------------------------------------------------------------------------
! Description   :
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'node_enter'
  character(len=*) :: node_name
  integer node_ind

! begin
  if (l_trace) then
    write(6,'(a,a)') 'l_trace: entering ',trim(node_name)
  endif

! index of node
  node_ind = node_index (node_name)

! if node not found, die
   if (node_ind == 0) then
      call die(lhere, 'entering node >'+trim(node_name)+'< that not catalogued.')
   endif

! die if we were already entered in this node but not exited
   if (nodes(node_ind)%entered) then
      call die(lhere, 'entering a second time node >'+trim(node_name)+'< but first exiting has not been notified.')
   endif

   nodes(node_ind)%entered = .true.

 end subroutine node_enter

!===============================================================================
  subroutine node_exit (node_name)
! ------------------------------------------------------------------------------
! Description   :
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout) :: here
  character(len=*) :: node_name
  integer node_ind

! begin
  here ='node_exit'

  if (l_trace) then
    write(6,'(a,a)') 'l_trace: exiting  ',trim(node_name)
  endif

! index of node
  node_ind = node_index (node_name)


! if node not found, die
   if (node_ind == 0) then
      write(6,*) trim(here),': exiting node ', trim(node_name),' that not catalogued'
      call die(here)
   endif

! if we have not entered this node previously
   if (.not. nodes(node_ind)%entered) then
      write(6,*) trim(here),': exiting node ', trim(node_name),' but entering has not been notified.'
      call die(here)
   endif

!   call cpu_time (nodes(node_ind)%cpu_end)
!
!   write(6,*) 'cpu_end=',nodes(node_ind)%cpu_end
!
   nodes(node_ind)%entered = .false.
!   nodes(node_ind)%cpu_duration = nodes(node_ind)%cpu_duration + (nodes(node_ind)%cpu_end - nodes(node_ind)%cpu_start)
!
 end subroutine node_exit


!===============================================================================
  subroutine nodes_statistics
! ------------------------------------------------------------------------------
! Description   : print statisitics for nodes
!
! Created       : J. Toulouse, 16 Dec 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout) :: here
  character(len=30) :: node_name
  integer nod_i, nod_j, nod_srt_i
  real(dp) temp
  real(dp) cpu_duration_total

! begin
  here ='nodes_statistics'

  cpu_duration_total = 0.d0

! sort by decreasing cpu times
  call alloc ('nodes_indexes_sort', nodes_indexes_sort, nodes_nb)
  call alloc ('nodes_cpu_duration_sort', nodes_cpu_duration_sort, nodes_nb)
  call alloc ('nodes_calls_nb_sort', nodes_calls_nb_sort, nodes_nb)

  do nod_i = 1, nodes_nb
   nodes_indexes_sort (nod_i) = nod_i
   nodes_cpu_duration_sort (nod_i) = nodes(nod_i)%cpu_duration
   nodes_calls_nb_sort (nod_i) = nodes(nod_i)%calls_nb
  enddo

  do nod_i = 1, nodes_nb
   do nod_j = nod_i+1, nodes_nb

   if (nodes_cpu_duration_sort (nod_j) > nodes_cpu_duration_sort (nod_i)) then
    temp = nodes_cpu_duration_sort (nod_i)
    nodes_cpu_duration_sort (nod_i) = nodes_cpu_duration_sort (nod_j)
    nodes_cpu_duration_sort (nod_j) = temp
    temp = nodes_indexes_sort (nod_i)
    nodes_indexes_sort (nod_i) = nodes_indexes_sort (nod_j)
    nodes_indexes_sort (nod_j) = temp
    temp = nodes_calls_nb_sort (nod_i)
    nodes_calls_nb_sort (nod_i) = nodes_calls_nb_sort (nod_j)
    nodes_calls_nb_sort (nod_j) = temp
   endif

   if (nodes_cpu_duration_sort (nod_j) == nodes_cpu_duration_sort (nod_i)) then
   if (nodes_calls_nb_sort (nod_j) > nodes_calls_nb_sort (nod_i)) then
    temp = nodes_cpu_duration_sort (nod_i)
    nodes_cpu_duration_sort (nod_i) = nodes_cpu_duration_sort (nod_j)
    nodes_cpu_duration_sort (nod_j) = temp
    temp = nodes_indexes_sort (nod_i)
    nodes_indexes_sort (nod_i) = nodes_indexes_sort (nod_j)
    nodes_indexes_sort (nod_j) = temp
    temp = nodes_calls_nb_sort (nod_i)
    nodes_calls_nb_sort (nod_i) = nodes_calls_nb_sort (nod_j)
    nodes_calls_nb_sort (nod_j) = temp
   endif
   endif

  enddo
  enddo

    write(6,*)
    write(6,'(a)') '----------------------------------------------------------------------'
    write(6,'(a)') '                        NODE STATISTICS'
    write(6,'(a)') ''
    write(6,'(a)') '    node                           calls           cpu time'
  do nod_i = 1, nodes_nb
   nod_srt_i = nodes_indexes_sort (nod_i)
   if (nodes(nod_srt_i)%calls_nb >= 1) then
   node_name = trim(nodes(nod_srt_i)%routine_name)
   write(6,'(a,i10,es15.8)') node_name, nodes(nod_srt_i)%calls_nb, nodes(nod_srt_i)%cpu_duration
    cpu_duration_total = cpu_duration_total +  nodes(nod_srt_i)%cpu_duration
   endif
  enddo
   write(6,'(a)') '------------------------------------------------------------------'
   write(6,'(a,es15.8)') 'TOTAL                                   ',cpu_duration_total
   write(6,*)

 end subroutine nodes_statistics

end module nodes_mod

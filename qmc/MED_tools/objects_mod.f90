module objects_mod

! size -> mysize for pathscale compiler

  use basic_tools_mod
  use constants_mod
  use nodes_mod
  use routines_mod
  use types_mod
  use strings_tools_mod

 type type_object
   character(len=max_string_len_obj) :: name
   logical                           :: associated
   character(len=max_string_len_type):: type
   integer, allocatable              :: dimensions (:)
   logical                           :: walkers = .false.
   logical                           :: unweighted = .false.
   logical                           :: valid
   logical                           :: freezed = .false.
   logical                           :: object_depend_valid = .false.
   integer                           :: node_create_index = 0
   integer, allocatable              :: nodes_depend_index (:)

!  pointers to objects
   integer,  pointer                 :: pointer_integer_0
   integer,  pointer                 :: pointer_integer_1 (:)
   integer,  pointer                 :: pointer_integer_2 (:,:)
   integer,  pointer                 :: pointer_integer_3 (:,:,:)
   integer,  pointer                 :: pointer_integer_4 (:,:,:,:)
   real(dp), pointer                 :: pointer_double_0
   real(dp), pointer                 :: pointer_double_1 (:)
   real(dp), pointer                 :: pointer_double_2 (:,:)
   real(dp), pointer                 :: pointer_double_3 (:,:,:)
   real(dp), pointer                 :: pointer_double_4 (:,:,:,:)
   real(dp), pointer                 :: pointer_double_5 (:,:,:,:,:)
   complex(dpc), pointer             :: pointer_complex_1 (:)
   logical, pointer                  :: pointer_logical_0
   logical, pointer                  :: pointer_logical_1 (:)
   logical, pointer                  :: pointer_logical_2 (:,:)
   logical, pointer                  :: pointer_logical_3 (:,:,:)
   character(len=max_string_len), pointer  :: pointer_string_1 (:)

!  for saving and restoring the value of object
   logical                           :: saved    = .false.
   real(dp)                          :: save_integer_0
   real(dp)                          :: save_double_0
   real(dp), allocatable             :: save_double_1 (:)
   real(dp), allocatable             :: save_double_2 (:,:)
   real(dp), allocatable             :: save_double_3 (:,:,:)

!  for average and statistical error
   real(dp)                          :: sum_double_0
   real(dp)                          :: sum_blk_double_0
   real(dp)                          :: sum_blk_square_double_0
   real(dp)                          :: previous_av1_double_0
   real(dp)                          :: previous_av2_double_0
   real(dp), allocatable             :: sum_double_1 (:)
   real(dp), allocatable             :: sum_blk_double_1 (:)
   real(dp), allocatable             :: sum_blk_square_double_1 (:)
   real(dp), allocatable             :: previous_av1_double_1 (:)
   real(dp), allocatable             :: previous_av2_double_1 (:)
   real(dp), allocatable             :: sum_double_2 (:,:)
   real(dp), allocatable             :: sum_blk_double_2 (:,:)
   real(dp), allocatable             :: sum_blk_square_double_2 (:,:)
   real(dp), allocatable             :: previous_av1_double_2 (:,:)
   real(dp), allocatable             :: previous_av2_double_2 (:,:)
   real(dp)                          :: variance_double_0
   complex(dpc), allocatable         :: sum_complex_1 (:)
   complex(dpc), allocatable         :: sum_blk_complex_1 (:)
   complex(dpc), allocatable         :: sum_blk_square_complex_1 (:)
   complex(dpc), allocatable         :: previous_av1_complex_1 (:)

  end type type_object


  integer                            :: objects_nb = 0
  type (type_object), save           :: objects (max_objects_nb)

  logical                            :: l_print_all_objects = .false.

  integer                            :: objects_sort_index (max_objects_nb)

  character(len=max_string_len_rout) :: current_node_exe = 'undefined'

! Interfaces

!===============================================================
  interface object_provide
!---------------------------------------------------------------
   module procedure object_provide_default,  &
                    object_provide_in_node

  end interface object_provide

!===============================================================
  interface object_provide_by_index
!---------------------------------------------------------------
   module procedure object_provide_default_by_index,  &
                    object_provide_in_node_by_index

  end interface object_provide_by_index

!===============================================================
  interface object_alloc
!---------------------------------------------------------------
   module procedure object_alloc_integer_1,  &
                    object_alloc_integer_2,  &
                    object_alloc_integer_3,  &
                    object_alloc_integer_row_1,  &
                    object_alloc_double_1,   &
                    object_alloc_double_2,   &
                    object_alloc_double_3,   &
                    object_alloc_double_4,   &
                    object_alloc_double_5,   &
                    object_alloc_double_row_1,   &
                    object_alloc_complex_1,   &
                    object_alloc_logical_1,  &
                    object_alloc_logical_2,  &
                    object_alloc_logical_3,  &
                    object_alloc_string_1,   &
                    object_alloc_string_row_1

  end interface object_alloc

!===============================================================
  interface object_associate
!---------------------------------------------------------------
   module procedure object_associate_integer_0, &
                    object_associate_integer_1, &
                    object_associate_integer_2, &
                    object_associate_integer_3, &
                    object_associate_integer_row_1, &
                    object_associate_double_0,  &
                    object_associate_double_1,  &
                    object_associate_double_2,  &
                    object_associate_double_3,  &
                    object_associate_double_4,  &
                    object_associate_double_5,  &
                    object_associate_double_row_1, &
                    object_associate_complex_1,  &
                    object_associate_logical_0, & !fp
                    object_associate_logical_1, &
                    object_associate_logical_2, &
                    object_associate_logical_3, &
                    object_associate_string_1,  &
                    object_associate_string_row_1

  end interface object_associate

!===============================================================
  interface object_release
!---------------------------------------------------------------
   module procedure object_release_integer_1,  &
                    object_release_integer_2,  &
                    object_release_integer_3,  &
                    object_release_integer_row_1, &
                    object_release_double_1,   &
                    object_release_double_2,   &
                    object_release_double_3,   &
                    object_release_double_4,   &
                    object_release_double_5,   &
                    object_release_double_row_1,   &
                    object_release_complex_1,  &
                    object_release_logical_1,  &
                    object_release_logical_2,  &
                    object_release_logical_3, &
                    object_release_string_1,  &
                    object_release_string_row_1

  end interface object_release

!===============================================================
  interface object_write
!---------------------------------------------------------------
   module procedure object_write, &
                    object_write_no_routine_name

  end interface object_write

!===============================================================
  interface object_write_2
!---------------------------------------------------------------
   module procedure object_write_2_no_routine_name

  end interface object_write_2

  contains

! ==============================================================================
  function object_index (object_name) result(result)
! ------------------------------------------------------------------------------
! Description   : return the index of an object, 0 if object not catalogued
! Description   : lookup is speeded up by ordering of objects
!
! Created       : J. Toulouse, 08 Jan 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! output
  integer                      :: result

! local
  integer obj_i, obj_1, obj_2
  integer obj_ind_i, obj_ind_1, obj_ind_2
  character(len=max_string_len_obj) current_object_name

! begin
  if (objects_nb == 0) then
    result = 0
    return
  endif

  obj_1 = 1
  obj_ind_1 = objects_sort_index (obj_1)
  if (object_name == objects(obj_ind_1)%name) then
    result = obj_ind_1
    return
  endif

  if (objects_nb == 1) then
   result = 0
   return
  endif

  obj_2 = objects_nb
  obj_ind_2 = objects_sort_index (obj_2)
  if (object_name == objects(obj_ind_2)%name) then
    result = obj_ind_2
    return
  endif

! loop over object indexes
    do
    obj_i = int((obj_1 + obj_2)/2.d0)

    if (obj_i == obj_1 .or.  obj_i == obj_2 ) then
     result = 0
     return
    endif

    obj_ind_i = objects_sort_index (obj_i)
    current_object_name = objects(obj_ind_i)%name
    if (object_name == current_object_name) then
      result = obj_ind_i
      return
    endif

    if (object_name < current_object_name) then
     obj_2 = obj_i
    else
     obj_1 = obj_i
    endif

    enddo ! end loop

  return
  end function object_index

! ==============================================================================
  function object_index_or_die (object_name) result(object_ind)
! ------------------------------------------------------------------------------
! Description   : return the index of an object, die if object not catalogued
!
! Created       : J. Toulouse, 20 Apr 2007
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! output
  integer object_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_index_or_die'

! begin

! index of object
  object_ind = object_index (object_name)

! if object not found, die
  if (object_ind == 0) then
      write(6,*)
      write(6,'(4a)') trim(lhere),': object >',trim(object_name),'< not catalogued. You should check that:'
      write(6,'(4a)') trim(lhere),': - this object has been declared as created object of a node, i.e. call object_create (',trim(object_name),')'
      write(6,'(2a)') trim(lhere),': - the creating node of this object has been catalogued in catalog_routines_mod.f90, i.e. call catalog_one_node (...)'
      write(6,'(2a)') trim(lhere),': or that'
      write(6,'(4a)') trim(lhere),': - this object has at least been declared as an object using call object_modified (',trim(object_name),')'
      call die (lhere)
  endif

  return
  end function object_index_or_die

! ==============================================================================
  subroutine object_add (object_name)
! ------------------------------------------------------------------------------
! Description   : catalogue a new object and order object list
!
! Created       : J. Toulouse, 18 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_add'
  integer obj_i, object_sort_pos

! begin
  if (objects_nb == max_objects_nb) then
   call die (lhere, 'maximum number of objects max_objects_nb='+max_objects_nb+' reached. Increase max_objects_nb in objects_mod.f90')
  endif

  objects_nb = objects_nb + 1
  objects(objects_nb)%name = object_name
  objects(objects_nb)%associated = .false.
  objects(objects_nb)%type = 'undefined'
  objects(objects_nb)%valid = .false.
  objects(objects_nb)%node_create_index = 0

! indexes of sorted objects
  object_sort_pos = objects_nb
  do obj_i = 1, objects_nb-1
    if (object_name < objects(objects_sort_index(obj_i))%name) then
      object_sort_pos = obj_i
      exit
    endif
  enddo
  do obj_i = objects_nb, object_sort_pos + 1, -1
      objects_sort_index (obj_i) = objects_sort_index (obj_i-1)
  enddo
  objects_sort_index (object_sort_pos) = objects_nb

  end subroutine object_add

! ==============================================================================
  subroutine object_alloc_by_index_double_0 (object_ind)
! ------------------------------------------------------------------------------
! Description   : allocate pointer of an object
!
! Created       : J. Toulouse, 06 Sep 2007
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_by_index_double_0'
  integer all_err

! begin
  allocate (objects(object_ind)%pointer_double_0, stat = all_err)
  if (all_err /= 0) then
   call die (lhere,'allocation of object >'+trim(objects(object_ind)%name)+'< failed.')
  endif
  objects(object_ind)%pointer_double_0 = 0.d0

  end subroutine object_alloc_by_index_double_0

! ==============================================================================
  subroutine object_add_and_index (object_ind)
! ------------------------------------------------------------------------------
! Description   : catalogue a new object (with no name) and return its index
!
! Created       : J. Toulouse, 06 Sep 2007
! ------------------------------------------------------------------------------
  implicit none

! output
  integer, intent(out) :: object_ind

! local
  character(len=max_string_len_obj) object_name
  integer object_number

! begin

! choose a default name
  object_number = objects_nb+1
  object_name = 'object'+object_number

! add object and return index
  call object_add_once_and_index (object_name, object_ind)

  end subroutine object_add_and_index

! ==============================================================================
  subroutine object_add_once_and_index (object_name, object_ind)
! ------------------------------------------------------------------------------
! Description   : catalogue a new object if not already catalogued
!
! Created       : J. Toulouse, 11 Nov 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! output
  integer, intent(out) :: object_ind

! local

! begin

! index of object
  object_ind = object_index (object_name)

! add object if not already catalogued
  if (object_ind == 0 ) then
    call object_add (object_name)
    object_ind = objects_nb
  endif

  end subroutine object_add_once_and_index

! ==============================================================================
  subroutine object_create (object_name, object_ind)
! ------------------------------------------------------------------------------
! Description   : catalog a new object created
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! output
  integer, intent(out), optional :: object_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_create'
  integer object_ind_local
  integer obj_i

! begin

! this routine must be called from a header of a node, i.e. for a defined node_current_index
  if (node_current_index == 0) then
      write(6,'(4a)') trim(lhere),': create node for object >', trim(object_name),'< not defined'
      write(6,'(4a)') trim(lhere),': this probably happened because the subroutine object_create was called from outside a header'
      call die (lhere)
  endif

! if object already catalogued as a created object, die
  object_ind_local = object_index (object_name)
  if (object_ind_local /= 0 ) then
     if (objects(object_ind_local)%node_create_index /= 0) then
       write(6,'(4a)') trim(lhere),': object >', trim(object_name),'< created at more than one place:'
       write(6,'(6a)') trim(lhere),': first at node >',trim(nodes(objects(object_ind_local)%node_create_index)%routine_name),'< and then again at node >',trim(nodes(node_current_index)%routine_name),'<'
       call die (lhere)
     endif
  endif

! check if object_ind not already set
  if (present(object_ind)) then
   if (object_ind /= 0 ) then
     write(6,'(4a,i3)') trim(lhere),': on entering the subroutine object_create, object ', trim(object_name),' has already a defined index different from 0'
     write(6,'(2a,i3)') trim(lhere),': object_ind= ',object_ind
     call die (lhere)
   endif
  endif

! add object if not already catalogued
  call object_add_once_and_index (object_name, object_ind_local)

! check if the current created object is not as the same time a needed object of the same node!
  do obj_i = 1, mysize(nodes(node_current_index)%objects_needed_index)
    if (object_ind_local == nodes(node_current_index)%objects_needed_index (obj_i) ) then
      write(6,'(5a)') trim(lhere),': object ', trim(object_name),' is needed and created by the same node ', trim(nodes(node_current_index)%routine_name)
      call die (lhere)
    endif
  enddo

! the current node creates the current object
  objects(object_ind_local)%node_create_index = node_current_index

! the current create object is created by the current node
  call append_once(nodes(node_current_index)%objects_create_index, object_ind_local)

! output object_ind
  if (present(object_ind)) then
   object_ind = object_ind_local
  endif

!  write(6,*) trim(lhere),': object ',trim(object_name)

  end subroutine object_create

! ==============================================================================
  subroutine object_needed (object_name)
! ------------------------------------------------------------------------------
! Description   : catalog a new object needed
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_needed'
  integer object_ind, obj_i

! begin

! this routine must be called from a header of a node, i.e. for a defined node_current_index
  if (node_current_index == 0) then
      write(6,'(4a)') trim(lhere),': needed node for object ', trim(object_name),' not defined'
      write(6,'(2a)') trim(lhere),': this probably happened because the subroutine object_needed was called from outside a header'
      call die (lhere)
  endif

! add object if not already catalogued
  call object_add_once_and_index (object_name, object_ind)

! check if the current needed object is not as the same time a created object of the same node!
  do obj_i = 1, mysize(nodes(node_current_index)%objects_create_index)
    if (object_ind == nodes(node_current_index)%objects_create_index (obj_i) ) then
      call die (lhere, 'object >'+trim(object_name)+'< is needed and created by the same node >'+trim(nodes(node_current_index)%routine_name)+'<.')
    endif
  enddo

! the current needed object is needed by the current node
  call append_once(nodes(node_current_index)%objects_needed_index, object_ind)

! the current node depends on the current needed object
  call append_once(objects(object_ind)%nodes_depend_index, node_current_index)

  end subroutine object_needed

! ==============================================================================
   function object_valid (object_name)
! ------------------------------------------------------------------------------
! Description   : object valid?
!
! Created       : J. Toulouse, 25 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! output
  logical object_valid

! local
  integer object_ind

! begin

! index of object
  call object_add_once_and_index (object_name, object_ind)

  object_valid = object_valid_by_index (object_ind)

  end function object_valid

! ==============================================================================
  function object_valid_by_index (object_index)
! ------------------------------------------------------------------------------
! Description   : object valid?
!
! Created       : J. Toulouse, 25 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_index

! output
  logical object_valid_by_index

! begin
  object_valid_by_index = objects(object_index)%valid

  end function object_valid_by_index

! ==============================================================================
  subroutine object_valid_or_die (object_name)
! ------------------------------------------------------------------------------
! Description   : check if object is valid, die if not
!
! Created       : J. Toulouse, 08 Nov 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_valid_or_die'

! begin
  if (.not. object_valid (object_name)) then
   call die (lhere, 'object >'+trim(object_name)+'< is not marked as valid.')
  endif

  end subroutine object_valid_or_die

! ==============================================================================
  subroutine object_associated_or_die_by_index (object_ind)
! ------------------------------------------------------------------------------
! Description   : check if object is associated, die if not
!
! Created       : J. Toulouse, 05 Feb 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) ::  object_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_associated_or_die_by_index'

! begin
  if (.not. objects(object_ind)%associated) then
   call die (lhere, 'object >'+trim(objects(object_ind)%name)+'< is not associated.')
  endif

  end subroutine object_associated_or_die_by_index

! ==============================================================================
  subroutine object_provide_default (object_name)
! ------------------------------------------------------------------------------
! Description   : provide an object
!
! Created       : J. Toulouse, 15 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  integer object_ind

! begin
!  write(6,'(3a)') trim(lhere),': calling object_provide for object ',trim(object_name)

! index of object
  object_ind = object_index_or_die (object_name)

  call object_provide_by_index (object_ind)

  end subroutine object_provide_default

! ==============================================================================
  subroutine object_provide_in_node (node_name, object_name)
! ------------------------------------------------------------------------------
! Description   : provide an object and add the denpendence of node on object
!
! Created       : J. Toulouse, 15 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: node_name, object_name

! local
  character(len=max_string_len_rout), save :: lhere ='object_provide_in_node'
  integer object_ind, node_ind

! begin
!  write(6,*) trim(lhere),': node: ',trim(node_name),', object: ',trim(object_name)

! index of object
  object_ind = object_index_or_die (object_name)

  call object_provide_by_index (object_ind)

! add dependency:

! index of node
  node_ind = node_index (node_name)

! if node not found, die
  if (node_ind == 0) then
     call die(lhere, 'node >'+trim(node_name)+'< not catalogued.')
  endif

! the current node depends on the current object
  call append_once(objects(object_ind)%nodes_depend_index, node_ind)

! add also the reverse dependency for object_depend_valid to work properly, but it can be problemic...
! the current object is needed by the current node
  call append_once(nodes(node_ind)%objects_needed_index, object_ind)

  end subroutine object_provide_in_node

! ==============================================================================
  subroutine object_provide_default_by_index (object_index)
! ------------------------------------------------------------------------------
! Description   : provide an object by its index
!
! Created       : J. Toulouse, 15 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_index

! local
  character(len=max_string_len_rout), save :: lhere = 'object_provide_default_by_index'

! begin

! if object is valid, do nothing
  if (objects(object_index)%valid) return

! if object has no create node, die
  if (objects(object_index)%node_create_index == 0) then
    write(6,*)
    write(6,'(4a)') trim(lhere),': object ', trim(objects(object_index)%name),' is marked as not valid and cannot be provided since it does not have a catalogued creating node.'
    if (trim(here) /= 'undefined') then
     write(6,'(4a)') trim(lhere),': This object is asked in routine ',trim(here),'.'
    endif
    write(6,'(2a)') trim(lhere),': Check that this object is declared as created in a node and that this node is catalogued in catalog_routines_mod.f90.'
    call die (lhere)
  endif

! execute create node of object
  call node_exe_by_index (objects(object_index)%node_create_index)

  end subroutine object_provide_default_by_index

! ==============================================================================
  subroutine object_provide_in_node_by_index (node_ind, object_ind)
! ------------------------------------------------------------------------------
! Description   : provide an object ny its index and add the denpendence of node on object
!
! Created       : J. Toulouse, 29 Mar 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: node_ind, object_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_provide_in_node_by_index'

! begin

! if object not catalogued, die
   if (object_ind == 0) then
     write(6,'(4a,i3)') trim(lhere),': in routine ', trim(nodes(node_ind)%routine_name),', an object is asked but has an undefined index = ', object_ind
     write(6,'(2a)') trim(lhere),': most likely the index of this object has not been given in its creation routine: call object_create (object_name, object_index)'
     call die (lhere)
   endif

  call object_provide_by_index (object_ind)

! add dependency:

! if node not found, die
  if (node_ind == 0) then
      call die(lhere, 'node index ='+node_ind)
  endif

! the current node depends on the current object
  call append_once(objects(object_ind)%nodes_depend_index, node_ind)

! add also the reverse dependency for object_depend_valid to work properly, but it can be problemic...
! the current object is needed by the current node
  call append_once(nodes(node_ind)%objects_needed_index, object_ind)

  end subroutine object_provide_in_node_by_index

! ==============================================================================
  subroutine object_provide_from_node_by_index (node_name, object_index)
! ------------------------------------------------------------------------------
! Description   : provide an object by its index, passing which node is asking it
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_index
  character(len=*), intent(in)   :: node_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_provide_from_node_by_index'

! begin

! if object is valid, do nothing
  if (objects(object_index)%valid) return

! if object has no create node, die
  if (objects(object_index)%node_create_index == 0) then
    write(6,*)
    write(6,'(4a)') trim(lhere),': object >',trim(objects(object_index)%name),'< is marked as not valid and cannot be provided since it does not have a catalogued creating node.'
    write(6,'(4a)') trim(lhere),': This object is asked by node >',trim(node_name),'<.'
    write(6,'(2a)') trim(lhere),': Check that this object is declared as created in a node and that this node is catalogued in catalog_routines_mod.f90'
    call die (lhere)
  endif

! execute create node of object
  call node_exe_by_index (objects(object_index)%node_create_index)

  end subroutine object_provide_from_node_by_index

!===============================================================================
  subroutine node_exe (node_name)
! ------------------------------------------------------------------------------
! Description   : execute a node by its name
!
! Created       : J. Toulouse, 15 Dec 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: node_name

! local
  character(len=max_string_len_rout), save  :: lhere =  'node_exe'
  integer node_ind

! index of node
  node_ind = node_index (node_name)

! if node not found, die
  if (node_ind == 0) then
    call die (lhere, 'node >'+trim(node_name)+'< not catalogued.')
  endif

  call node_exe_by_index (node_ind)

  end subroutine node_exe

!===============================================================================
  recursive subroutine node_exe_by_index (node_index)
! ------------------------------------------------------------------------------
! Description   : execute a node by its index
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: node_index

! local
  integer obj_i
  real(dp) cpu1, cpu2

! begin
! if the current node is valid, do nothing
  if (nodes(node_index)%valid) return

! provide all the needed objects of the current nodes (we cannot simply execute the needed nodes
! since some objects may not have a creating node)
  do obj_i = 1, mysize(nodes(node_index)%objects_needed_index)
   call object_provide_from_node_by_index (nodes(node_index)%routine_name, nodes(node_index)%objects_needed_index(obj_i))
  enddo

# if defined (DEBUG)
    if (l_trace) then
      write(6,'(a,i3,2a)') 'trace: execute node # ',node_index,': ', trim(nodes(node_index)%routine_name)
    endif
# endif

! set name of routine
  here = nodes(node_index)%routine_name

! execute the current node
  call cpu_time(cpu1)
  call exe_by_address_0 (nodes(node_index)%routine_address)
  call cpu_time(cpu2)

! reset name of routine
  here = 'undefined'

  nodes(node_index)%calls_nb = nodes(node_index)%calls_nb + 1
  nodes(node_index)%cpu_duration =  nodes(node_index)%cpu_duration + cpu2 - cpu1

! validate the current node
  nodes(node_index)%valid = .true.

! validate objects created by the current node
  do obj_i = 1, mysize(nodes(node_index)%objects_create_index)
   objects(nodes(node_index)%objects_create_index(obj_i))%valid = .true.
  enddo

! print all objects
# if defined (DEBUG)
  if (l_print_all_objects) then
   do obj_i = 1, mysize(nodes(node_index)%objects_create_index)
    call object_write_by_index (nodes(node_index)%objects_create_index(obj_i))
   enddo
  endif
# endif

  end subroutine node_exe_by_index

! ==============================================================================
  subroutine object_modified (object_name)
! ------------------------------------------------------------------------------
! Description   : when an object has been modified, validate it and
! Description   : invalidate all objects depending on it
! Description   : also catalog object if necessary
!
! Created       : J. Toulouse, 15 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)             :: object_name

! local
  integer                                  :: object_ind

! begin
!  write(6,*) trim(lhere),': modify object by name ', trim(object_name)

! index of object, and catalog it if necessary
  call object_add_once_and_index (object_name, object_ind)

  call object_modified_by_index (object_ind)

  end subroutine object_modified

! ==============================================================================
  subroutine object_modified_by_index (object_index)
! ------------------------------------------------------------------------------
! Description   : when an object has been modified, invalidate all objects depending on it
! Description   : using the index of the object
!
! Created       : J. Toulouse, 15 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_index

! local
  integer node_i, obj_i

! begin

! temporary
!  call object_modified2_by_index (object_index) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  return   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! validate current object
  objects(object_index)%valid = .true.

! improve: the following part takes time if a lot of calls to object_modified_by_index

! invalidate objects created by the nodes depending of the current object
! (if there is no nodes depending on this object then nodes_depend_nb = 0)
  do node_i = 1, mysize(objects(object_index)%nodes_depend_index)
   do obj_i = 1, mysize(nodes(objects(object_index)%nodes_depend_index(node_i))%objects_create_index)
     call object_invalidate_by_index (nodes(objects(object_index)%nodes_depend_index(node_i))%objects_create_index(obj_i))
   enddo
  enddo

  end subroutine object_modified_by_index

! ==============================================================================
  subroutine object_modified2_by_index (object_index)
! ------------------------------------------------------------------------------
! Description   : when an object has been modified, invalidate all objects depending on it
! Description   : using the index of the object
! Description   : and mark that the current object has been validated for the objects
! Description   : on which the current object depends
! Description   : This routine is more general and safe than object_modified_by_index
! Description   : but it is (a bit?) more costly, so it is not yet used everywhere
!
! Created       : J. Toulouse, 13 Apr 2007
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_index

! local
  integer node_i, obj_i

! begin

! validate current object
  objects(object_index)%valid = .true.

! improve: the following part takes time if a lot of calls to object_modified_by_index

! invalidate objects created by the nodes depending of the current object
! (if there is no nodes depending on this object then nodes_depend_nb = 0)
  do node_i = 1, mysize(objects(object_index)%nodes_depend_index)
   do obj_i = 1, mysize(nodes(objects(object_index)%nodes_depend_index(node_i))%objects_create_index)
     call object_invalidate_by_index (nodes(objects(object_index)%nodes_depend_index(node_i))%objects_create_index(obj_i))
   enddo
  enddo

! mark that the current object has been validated for the objects on which the current object depends
! (in order to invalidate properly this object later on)
  if (objects(object_index)%node_create_index /= 0) then
   do obj_i = 1, mysize(nodes(objects(object_index)%node_create_index)%objects_needed_index)
     call object_depend_valid_by_index (nodes(objects(object_index)%node_create_index)%objects_needed_index (obj_i))
   enddo
  endif

  end subroutine object_modified2_by_index

! ==============================================================================
  recursive subroutine object_depend_valid_by_index (object_index)
! ------------------------------------------------------------------------------
! Description   : record that an object depending on the current object
! Description   : has been valididated, in order to properley invalidate it later on
!
! Created       : J. Toulouse, 13 Apr 2007
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_index

! local
  integer obj_i

! begin
!  write(6,'(4a)') trim(lhere),': >',trim(objects(object_index)%name),'<.'
  objects(object_index)%object_depend_valid = .true.

! recusrive call to objects on which the current object depends
  if (objects(object_index)%node_create_index /= 0) then
   do obj_i = 1, mysize(nodes(objects(object_index)%node_create_index)%objects_needed_index)
     call object_depend_valid_by_index (nodes(objects(object_index)%node_create_index)%objects_needed_index (obj_i))
   enddo
  endif

  end subroutine object_depend_valid_by_index

! ==============================================================================
  subroutine object_invalidate (object_name)
! ------------------------------------------------------------------------------
! Description   : invalidate an object and catalogue it if necessary
!
! Created       : J. Toulouse, 15 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  integer object_ind

! begin
!  write (6,'(3a)') 'invalidate object >',trim(object_name),'<.'

! index of object, and catalog it if necessary
  call object_add_once_and_index (object_name, object_ind)

  call object_invalidate_by_index (object_ind)

  end subroutine object_invalidate

! ==============================================================================
  recursive subroutine object_invalidate_by_index (object_index)
! ------------------------------------------------------------------------------
! Description   : invalidate an object by its index and all objects depending on it
!
! Created       : J. Toulouse, 15 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_index

! local
  integer node_i, obj_i


! old: if current object is already invalid, do nothing
!  if (.not. objects(object_index)%valid) return

! new: return only if current object not valid AND not valid object depending on the current object
  if (.not. (objects(object_index)%valid .or. objects(object_index)%object_depend_valid)) then
!   write (6,'(3a)') 'object >',trim(objects(object_index)%name),'< is already invalid.'
   return
  endif

!  write (6,'(3a)') 'invalidate object >',trim(objects(object_index)%name),'<.'

! if object has been asked to be freezed, leave the object valid and return
  if (objects(object_index)%freezed) return

! invalidate current object unless the object has been asked to be freezed
  objects(object_index)%valid = .false.
  objects(object_index)%object_depend_valid = .false.

! invalidate create node of current object if it exists
  if (objects(object_index)%node_create_index /= 0 ) then
   nodes(objects(object_index)%node_create_index)%valid = .false.
  endif

! invalidate objects created by the nodes depending of the current object
! (if there are no nodes depending on this object then nodes_depend_nb = 0)
  do node_i = 1, mysize(objects(object_index)%nodes_depend_index)
   do obj_i = 1, mysize(nodes(objects(object_index)%nodes_depend_index(node_i))%objects_create_index)
     call object_invalidate_by_index(nodes(objects(object_index)%nodes_depend_index(node_i))%objects_create_index(obj_i))
   enddo
  enddo

  end subroutine object_invalidate_by_index

!===========================================================================
  subroutine object_associate_integer_0 (object_name, object)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, target, intent(in)         :: object

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'integer_0'

! associate pointer
  objects(object_ind)%pointer_integer_0 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_integer_0

!===========================================================================
  subroutine object_associate_integer_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, target, intent(in)         :: object (:)
  integer, intent(in)                 :: dim1

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type  = 'integer_1'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)

! associate pointer
  objects(object_ind)%pointer_integer_1 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_integer_1

!===========================================================================
  subroutine object_associate_integer_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)       :: object_name
  integer, target, intent(in)        :: object (:,:)
  integer, intent(in)                :: dim1, dim2

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'integer_2'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)
  call append(objects(object_ind)%dimensions, dim2)

! associate pointer
  objects(object_ind)%pointer_integer_2 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_integer_2

!===========================================================================
  subroutine object_associate_integer_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, target, intent(in)         :: object (:,:,:)
  integer, intent(in)                 :: dim1, dim2, dim3

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'integer_3'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)
  call append(objects(object_ind)%dimensions, dim2)
  call append(objects(object_ind)%dimensions, dim3)

! associate pointer
  objects(object_ind)%pointer_integer_3 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_integer_3

!===========================================================================
  subroutine object_associate_integer_row_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                   :: object_name
  type (type_integer_row), intent(in), target    :: object (:)
  integer, intent(in)                            :: dim1

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'integer_row_1'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)

!  associate pointer
!  objects(object_ind)%pointer_integer_row_1 => object
!  objects(object_ind)%associated = .true.

  end subroutine object_associate_integer_row_1


!===========================================================================
  subroutine object_associate_double_0 (object_name, object)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name
  real(dp), target, intent(in)            :: object

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'double_0'

! associate pointer
  objects(object_ind)%pointer_double_0 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_double_0

!===========================================================================
  subroutine object_associate_double_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  real(dp), target, intent(in)        :: object (:)
  integer, intent(in)                 :: dim1

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'double_1'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)

! associate pointer
  objects(object_ind)%pointer_double_1 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_double_1

!===========================================================================
  subroutine object_associate_double_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  real(dp), target, intent(in)        :: object (:,:)
  integer, intent(in)                 :: dim1, dim2

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'double_2'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)
  call append(objects(object_ind)%dimensions, dim2)

! associate pointer
  objects(object_ind)%pointer_double_2 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_double_2

!===========================================================================
  subroutine object_associate_double_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  real(dp), target, intent(in)        :: object (:,:,:)
  integer, intent(in)                 :: dim1, dim2, dim3

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'double_3'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)
  call append(objects(object_ind)%dimensions, dim2)
  call append(objects(object_ind)%dimensions, dim3)

! associate pointer
  objects(object_ind)%pointer_double_3 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_double_3

!===========================================================================
  subroutine object_associate_double_4 (object_name, object, dim1, dim2, dim3, dim4)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  real(dp), target, intent(in)        :: object (:,:,:,:)
  integer, intent(in)                 :: dim1, dim2, dim3, dim4

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'double_4'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)
  call append(objects(object_ind)%dimensions, dim2)
  call append(objects(object_ind)%dimensions, dim3)
  call append(objects(object_ind)%dimensions, dim4)

! associate pointer
  objects(object_ind)%pointer_double_4 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_double_4

!===========================================================================
  subroutine object_associate_double_5 (object_name, object, dim1, dim2, dim3, dim4, dim5)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 26 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  real(dp), target, intent(in)        :: object (:,:,:,:,:)
  integer, intent(in)                 :: dim1, dim2, dim3, dim4, dim5

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'double_5'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)
  call append(objects(object_ind)%dimensions, dim2)
  call append(objects(object_ind)%dimensions, dim3)
  call append(objects(object_ind)%dimensions, dim4)
  call append(objects(object_ind)%dimensions, dim5)

! associate pointer
  objects(object_ind)%pointer_double_5 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_double_5

!===========================================================================
  subroutine object_associate_double_row_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)               :: object_name
  type (type_real_row), target, intent(in)   :: object (:)
  integer, intent(in)                        :: dim1

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'double_row_1'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)

!  associate pointer
!  objects(object_ind)%pointer_integer_row_1 => object
!  objects(object_ind)%associated = .true.

  end subroutine object_associate_double_row_1

!===========================================================================
  subroutine object_associate_complex_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 28 Sep 2013
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  complex(dpc), target, intent(in)    :: object (:)
  integer, intent(in)                 :: dim1

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'complex_1'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)

! associate pointer
  objects(object_ind)%pointer_complex_1 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_complex_1

!===========================================================================
  subroutine object_associate_logical_0 (object_name, object) !fp
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : F. Petruzielo, 21 Jul 2008
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  logical, target, intent(in)         :: object

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'logical_0'

! associate pointer
  objects(object_ind)%pointer_logical_0 => object
  objects(object_ind)%associated = .true.

 end subroutine object_associate_logical_0

!===========================================================================
  subroutine object_associate_logical_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  logical, target, intent(in)         :: object (:)
  integer, intent(in)                 :: dim1

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'logical_1'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)

! associate pointer
  objects(object_ind)%pointer_logical_1 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_logical_1

!===========================================================================
  subroutine object_associate_logical_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)   :: object_name
  logical, target, intent(in)    :: object (:,:)
  integer, intent(in)            :: dim1, dim2

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'logical_2'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)
  call append(objects(object_ind)%dimensions, dim2)

! associate pointer
  objects(object_ind)%pointer_logical_2 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_logical_2

!===========================================================================
  subroutine object_associate_logical_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)           :: object_name
  logical, target, intent(in)            :: object (:,:,:)
  integer, intent(in)                    :: dim1, dim2, dim3

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'logical_3'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)
  call append(objects(object_ind)%dimensions, dim2)
  call append(objects(object_ind)%dimensions, dim3)

! associate pointer
  objects(object_ind)%pointer_logical_3 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_logical_3

!===========================================================================
  subroutine object_associate_string_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 31 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                         :: object_name
  character(len=max_string_len), target, intent(in)    :: object (:)
  integer, intent(in)                                  :: dim1

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type    = 'string_1'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)

! associate pointer
  objects(object_ind)%pointer_string_1 => object
  objects(object_ind)%associated = .true.

  end subroutine object_associate_string_1

!===========================================================================
  subroutine object_associate_string_row_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer to object
!
! Created     : J. Toulouse, 11 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)               :: object_name
  type (type_string_row), target, intent(in)   :: object (:)
  integer, intent(in)                        :: dim1

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'string_row_1'

! store dimensions
  call append(objects(object_ind)%dimensions, dim1)

!  associate pointer
!  objects(object_ind)%pointer_integer_row_1 => object
!  objects(object_ind)%associated = .true.

  end subroutine object_associate_string_row_1

!===========================================================================
  subroutine object_associate_by_index_double_0 (object_ind)
!---------------------------------------------------------------------------
! Description : associate pointer of object by its index
!
! Created     : J. Toulouse, 06 Sep 2007
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in)  :: object_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_associate_by_index_double_0'
  integer all_err

! begin

! if object already associated, return
  if (objects(object_ind)%associated) return

! store type
  objects(object_ind)%type = 'double_0'

! associate (allocate) pointer
  allocate (objects(object_ind)%pointer_double_0, stat = all_err)
  if (all_err /= 0) then
   call die (lhere,'allocation of object >'+trim(objects(object_ind)%name)+'< failed.')
  endif
  objects(object_ind)%pointer_double_0 = 0.d0
  objects(object_ind)%associated = .true.

  end subroutine object_associate_by_index_double_0

!===========================================================================
  subroutine object_associate_by_index_double_1 (object_ind, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer of object by its index
!
! Created     : J. Toulouse, 07 Sep 2007
! Modified    : J. Toulouse, 14 Mar 2010, resize
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in)  :: object_ind, dim1

! local
  character(len=max_string_len_rout), save :: lhere = 'object_associate_by_index_double_1'
  integer all_err, object_dim, dim_min
  real(dp), allocatable                :: object_temp (:)

! begin

! associate object if not already associated
  if (.not. objects(object_ind)%associated) then
!   store type
    objects(object_ind)%type = 'double_1'
!   store dimensions
    call append(objects(object_ind)%dimensions, dim1)
!   associate (allocate) pointer
    allocate (objects(object_ind)%pointer_double_1(dim1), stat = all_err)
    if (all_err /= 0) then
     call die (lhere,'allocation of object >'+trim(objects(object_ind)%name)+'< failed.')
    endif
    objects(object_ind)%pointer_double_1 = 0.d0
    objects(object_ind)%associated = .true.

! resize object if already associated with different dimension
  else
    object_dim = size(objects(object_ind)%pointer_double_1)
    if (object_dim /= dim1) then
     dim_min =  min(object_dim, dim1)
     call alloc ('object_temp', object_temp, dim_min)
     object_temp(:) = objects(object_ind)%pointer_double_1(1:dim_min)
     call object_deassociate (objects(object_ind)%name)

     objects(object_ind)%type = 'double_1'
     call append(objects(object_ind)%dimensions, dim1)
     allocate (objects(object_ind)%pointer_double_1(dim1), stat = all_err)
     if (all_err /= 0) then
      call die (lhere,'allocation of object >'+trim(objects(object_ind)%name)+'< failed.')
     endif
     objects(object_ind)%pointer_double_1 = 0.d0
     objects(object_ind)%associated = .true.

     objects(object_ind)%pointer_double_1(1:dim_min) = object_temp(:)
     call release ('object_temp', object_temp)
!    resize also the corresponding objects for averages and errors
!     call alloc ('objects(object_ind)%sum_double_1',objects(object_ind)%sum_double_1, dim1)
!     call alloc ('objects(object_ind)%sum_blk_double_1',objects(object_ind)%sum_blk_double_1, dim1)
!!     call alloc ('objects(object_ind)%previous_double_1',objects(object_ind)%previous_double_1, dim1)
    endif
  endif

  end subroutine object_associate_by_index_double_1

!===========================================================================
  subroutine object_associate_by_index_double_2 (object_ind, dim1, dim2)
!---------------------------------------------------------------------------
! Description : associate pointer of object by its index
!
! Created     : J. Toulouse, 07 Sep 2007
! Modified    : J. Toulouse, 14 Mar 2010, resize
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in)  :: object_ind, dim1, dim2

! local
  character(len=max_string_len_rout), save :: lhere = 'object_associate_by_index_double_2'
  integer all_err, object_dim1, object_dim2, dim_min1, dim_min2
  real(dp), allocatable                :: object_temp (:,:)

! begin

! associate object if not already associated
  if (.not. objects(object_ind)%associated) then
!   store type
    objects(object_ind)%type = 'double_2'
!   store dimensions
    call append(objects(object_ind)%dimensions, dim1)
    call append(objects(object_ind)%dimensions, dim2)
!   associate (allocate) pointer
    allocate (objects(object_ind)%pointer_double_2(dim1, dim2), stat = all_err)
    if (all_err /= 0) then
     call die (lhere,'allocation of object >'+trim(objects(object_ind)%name)+'< failed.')
    endif
    objects(object_ind)%pointer_double_2 = 0.d0
    objects(object_ind)%associated = .true.

! resize object if already associated with different dimensions
  else
    object_dim1 = size(objects(object_ind)%pointer_double_2,1)
    object_dim2 = size(objects(object_ind)%pointer_double_2,2)
    if (object_dim1 /= dim1 .or. object_dim2 /= dim2) then
     dim_min1 =  min(object_dim1, dim1)
     dim_min2 =  min(object_dim2, dim2)
     call alloc ('object_temp', object_temp, dim_min1, dim_min2)
     object_temp(:,:) = objects(object_ind)%pointer_double_2(1:dim_min1,1:dim_min2)
     call object_deassociate (objects(object_ind)%name)

     objects(object_ind)%type = 'double_2'
     call append(objects(object_ind)%dimensions, dim1)
     call append(objects(object_ind)%dimensions, dim2)
     allocate (objects(object_ind)%pointer_double_2(dim1, dim2), stat = all_err)
     if (all_err /= 0) then
      call die (lhere,'allocation of object >'+trim(objects(object_ind)%name)+'< failed.')
     endif
     objects(object_ind)%pointer_double_2 = 0.d0
     objects(object_ind)%associated = .true.

     objects(object_ind)%pointer_double_2(1:dim_min1,1:dim_min2) = object_temp(:,:)
     call release ('object_temp', object_temp)
!    resize also the corresponding objects for averages and errors
!   resize also the corresponding objects for averages and errors
!    object_ind = object_index (object_name)
!    call alloc ('objects(object_ind)%sum_double_1',objects(object_ind)%sum_double_1, dim1)
!    call alloc ('objects(object_ind)%sum_blk_double_1',objects(object_ind)%sum_blk_double_1, dim1)
!!    call alloc ('objects(object_ind)%previous_double_1',objects(object_ind)%previous_double_1, dim1)
!    call alloc ('objects(object_ind)%sum_double_2',objects(object_ind)%sum_double_2, dim1, dim2)
!    call alloc ('objects(object_ind)%sum_blk_double_2',objects(object_ind)%sum_blk_double_2, dim1, dim2)
!!    call alloc ('objects(object_ind)%previous_double_2',objects(object_ind)%previous_double_2, dim1, dim2)
    endif
  endif

  end subroutine object_associate_by_index_double_2

!===========================================================================
  subroutine object_associate_by_index_complex_1 (object_ind, dim1)
!---------------------------------------------------------------------------
! Description : associate pointer of object by its index
!
! Created     : J. Toulouse, 28 Sep 2013
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in)  :: object_ind, dim1

! local
  character(len=max_string_len_rout), save :: lhere = 'object_associate_by_index_complex_1'
  integer all_err, object_dim, dim_min
  complex(dpc), allocatable                :: object_temp (:)

! begin

! associate object if not already associated
  if (.not. objects(object_ind)%associated) then
!   store type
    objects(object_ind)%type = 'complex_1'
!   store dimensions
    call append(objects(object_ind)%dimensions, dim1)
!   associate (allocate) pointer
    allocate (objects(object_ind)%pointer_complex_1(dim1), stat = all_err)
    if (all_err /= 0) then
     call die (lhere,'allocation of object >'+trim(objects(object_ind)%name)+'< failed.')
    endif
    objects(object_ind)%pointer_complex_1 = 0.d0
    objects(object_ind)%associated = .true.

! resize object if already associated with different dimension
  else
    object_dim = size(objects(object_ind)%pointer_complex_1)
    if (object_dim /= dim1) then
     dim_min =  min(object_dim, dim1)
     call alloc ('object_temp', object_temp, dim_min)
     object_temp(:) = objects(object_ind)%pointer_complex_1(1:dim_min)
     call object_deassociate (objects(object_ind)%name)

     objects(object_ind)%type = 'complex_1'
     call append(objects(object_ind)%dimensions, dim1)
     allocate (objects(object_ind)%pointer_complex_1(dim1), stat = all_err)
     if (all_err /= 0) then
      call die (lhere,'allocation of object >'+trim(objects(object_ind)%name)+'< failed.')
     endif
     objects(object_ind)%pointer_complex_1 = 0.d0
     objects(object_ind)%associated = .true.

     objects(object_ind)%pointer_complex_1(1:dim_min) = object_temp(:)
     call release ('object_temp', object_temp)
!    resize also the corresponding objects for averages and errors
!     call alloc ('objects(object_ind)%sum_complex_1',objects(object_ind)%sum_complex_1, dim1)
!     call alloc ('objects(object_ind)%sum_blk_complex_1',objects(object_ind)%sum_blk_complex_1, dim1)
!!     call alloc ('objects(object_ind)%previous_complex_1',objects(object_ind)%previous_complex_1, dim1)
    endif
  endif

  end subroutine object_associate_by_index_complex_1

!===========================================================================
  subroutine object_deassociate (object_name)
!---------------------------------------------------------------------------
! Description : deassociate object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                       :: object_name

! local
  integer object_ind

! begin

! index of object, catalog if necessary
  call object_add_once_and_index (object_name, object_ind)

  objects(object_ind)%type = 'undefined'
  objects(object_ind)%associated = .false.
  deallocate(objects(object_ind)%dimensions)

  end subroutine object_deassociate

!===========================================================================
!  recursive subroutine object_alloc_integer_1 (object_name, object, dim1) !old
  subroutine object_alloc_integer_1 (object_name, object, dim1) !new
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1

! input/output
  integer, allocatable, intent(inout) :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_integer_1'
  integer all_err
  integer object_dim, dim_min
  integer, allocatable                :: object_temp (:)

! begin

! allocate object if not already allocated
  if (.not. allocated(object)) then
   allocate (object(dim1), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,i8)') trim(lhere),': dimension of object is ', dim1
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:) = 0
   call object_associate (object_name, object, dim1)

!  resize object if already allocated with different dimension
   else
   object_dim = size(object)
   if (object_dim /= dim1) then
    dim_min =  min(object_dim, dim1)
    call alloc ('object_temp', object_temp, dim_min)
    object_temp(:) = object(1:dim_min)
    call object_release (object_name, object)
!    call object_alloc (object_name, object, dim1) ! old
    allocate (object(dim1)) !new
    object(:) = 0 ! new
    call object_associate (object_name, object, dim1) ! new
    object(1:dim_min) = object_temp(:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_integer_1

!===========================================================================
  subroutine object_alloc_integer_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1, dim2

! input/output
  integer, allocatable, intent(inout) :: object (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_integer_2'
  integer all_err
  integer object_dim1, dim_min1
  integer object_dim2, dim_min2
  integer, allocatable                :: object_temp (:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1,dim2), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,2i8)') trim(lhere),': dimensions are ', dim1, dim2
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:,:) = 0
   call object_associate (object_name, object, dim1, dim2)

   else
!  resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   if (object_dim1 /= dim1 .or. object_dim2 /= dim2) then
    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    call alloc ('object_temp', object_temp, dim_min1, dim_min2)
    object_temp(:,:) = object(1:dim_min1,1:dim_min2)
    call object_release (object_name, object)
    allocate (object(dim1,dim2))
    object(:,:) = 0
    call object_associate (object_name, object, dim1, dim2)
    object(1:dim_min1,1:dim_min2) = object_temp(:,:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_integer_2

!===========================================================================
  subroutine object_alloc_integer_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1, dim2, dim3

! input/output
  integer, allocatable, intent(inout) :: object (:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_integer_3'
  integer all_err
  integer object_dim1, dim_min1
  integer object_dim2, dim_min2
  integer object_dim3, dim_min3
  integer, allocatable                 :: object_temp (:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1,dim2,dim3), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,3i8)') trim(lhere),': dimensions are ', dim1, dim2, dim3
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:,:,:) = 0
   call object_associate (object_name, object, dim1, dim2, dim3)

   else
!  resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3) then
    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)
    call alloc ('object_temp', object_temp, dim_min1, dim_min2, dim_min3)
    object_temp(:,:,:) = object(1:dim_min1,1:dim_min2,1:dim_min3)
    call object_release (object_name, object)
    allocate (object(dim1,dim2,dim3))
    object(:,:,:) = 0
    call object_associate (object_name, object, dim1, dim2, dim3)
    object(1:dim_min1,1:dim_min2,1:dim_min3) = object_temp(:,:,:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_integer_3

!===========================================================================
  subroutine object_alloc_integer_row_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1

! input/output
  type (type_integer_row), allocatable, intent(inout)  :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_integer_row_1'
  integer all_err
  integer object_dim, dim_min
  type (type_integer_row), allocatable          :: object_temp (:)

! begin
# if defined (DEBUG)
   call routine_enter (lhere)
# endif

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,i8)') trim(lhere),': dimension is ', dim1
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   call object_associate (object_name, object, dim1)

   else

!  resize object if already allocated with different dimension
   object_dim = size(object)
   if(object_dim /= dim1) then
    dim_min =  min(object_dim, dim1)
    call alloc ('object_temp', object_temp, dim_min)
    object_temp(:) = object(1:dim_min)
    call object_release (object_name, object)
    allocate (object(dim1))
    call object_associate (object_name, object, dim1)
    object(1:dim_min) = object_temp(:)
    call release ('object_temp', object_temp)
   endif

  endif

# if defined (DEBUG)
   call routine_exit (lhere)
# endif

  end subroutine object_alloc_integer_row_1

!===========================================================================
  subroutine object_alloc_double_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)         :: object_name
  integer, intent(in)                  :: dim1

! input/output
  real(dp), allocatable, intent(inout) :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_double_1'
  integer all_err, object_dim, dim_min, object_ind
  real(dp), allocatable                :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,i8)') trim(lhere),': dimension is ', dim1
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:) = 0.d0
   call object_associate (object_name, object, dim1)

   else
!  resize object if already allocated with different dimension
   object_dim = size(object)
   if (object_dim /= dim1) then
    dim_min =  min(object_dim, dim1)
    call alloc ('object_temp', object_temp, dim_min)
    object_temp(:) = object(1:dim_min)
    call object_release (object_name, object)
    allocate (object(dim1))
    object(:) = 0.d0
    call object_associate (object_name, object, dim1)
    object(1:dim_min) = object_temp(:)
    call release ('object_temp', object_temp)
!   resize also the corresponding objects for averages and errors
    object_ind = object_index (object_name)
    call alloc ('objects(object_ind)%sum_double_1',objects(object_ind)%sum_double_1, dim1)
    call alloc ('objects(object_ind)%sum_blk_double_1',objects(object_ind)%sum_blk_double_1, dim1)
!    call alloc ('objects(object_ind)%previous_double_1',objects(object_ind)%previous_double_1, dim1)
   endif

  endif

  end subroutine object_alloc_double_1

!===========================================================================
  subroutine object_alloc_double_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)         :: object_name
  integer, intent(in)                  :: dim1, dim2

! input/output
  real(dp), allocatable, intent(inout) :: object (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_double_2'
  integer all_err, object_ind
  integer object_dim1, dim_min1
  integer object_dim2, dim_min2
  real(dp), allocatable                :: object_temp (:,:)

! begin

! allocate object if not already allocated
  if (.not. allocated(object)) then

   allocate (object(dim1,dim2), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,2i8)') trim(lhere),': dimensions are ', dim1, dim2
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:,:) = 0.d0
   call object_associate (object_name, object, dim1, dim2)

   else
!  resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   if (object_dim1 /= dim1 .or. object_dim2 /= dim2) then
    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    call alloc ('object_temp', object_temp, dim_min1, dim_min2)
    object_temp(:,:) = object(1:dim_min1,1:dim_min2)
    call object_release (object_name, object)
    allocate (object(dim1,dim2))
    object(:,:) = 0.d0
    call object_associate (object_name, object, dim1, dim2)
    object(1:dim_min1,1:dim_min2) = object_temp(:,:)
    call release ('object_temp', object_temp)
!   resize also the corresponding objects for averages and errors
    object_ind = object_index (object_name)
    call alloc ('objects(object_ind)%sum_double_1',objects(object_ind)%sum_double_1, dim1)
    call alloc ('objects(object_ind)%sum_blk_double_1',objects(object_ind)%sum_blk_double_1, dim1)
!    call alloc ('objects(object_ind)%previous_double_1',objects(object_ind)%previous_double_1, dim1)
    call alloc ('objects(object_ind)%sum_double_2',objects(object_ind)%sum_double_2, dim1, dim2)
    call alloc ('objects(object_ind)%sum_blk_double_2',objects(object_ind)%sum_blk_double_2, dim1, dim2)
!    call alloc ('objects(object_ind)%previous_double_2',objects(object_ind)%previous_double_2, dim1, dim2)
   endif

  endif

  end subroutine object_alloc_double_2

!===========================================================================
  subroutine object_alloc_double_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)         :: object_name
  integer, intent(in)                  :: dim1, dim2, dim3

! input/output
  real(dp), allocatable, intent(inout) :: object (:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_double_3'
  integer all_err, object_ind
  integer object_dim1, dim_min1
  integer object_dim2, dim_min2
  integer object_dim3, dim_min3
  real(dp), allocatable                :: object_temp (:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1,dim2,dim3), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,3i8)') trim(lhere),': dimensions are ', dim1, dim2, dim3
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:,:,:) = 0.d0
   call object_associate (object_name, object, dim1, dim2, dim3)

   else

!  resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)
   if (object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3) then
    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)
    call alloc ('object_temp', object_temp, dim_min1, dim_min2, dim_min3)
    object_temp(:,:,:) = object(1:dim_min1,1:dim_min2,1:dim_min3)
    call object_release (object_name, object)
    allocate (object(dim1,dim2,dim3))
    object(:,:,:) = 0.d0
    call object_associate (object_name, object, dim1, dim2, dim3)
    object(1:dim_min1,1:dim_min2,1:dim_min3) = object_temp(:,:,:)
    call release ('object_temp', object_temp)
!   resize also the corresponding objects for averages and errors
    object_ind = object_index (object_name)
    call alloc ('objects(object_ind)%sum_double_2',objects(object_ind)%sum_double_2, dim1, dim2)
    call alloc ('objects(object_ind)%sum_blk_double_2',objects(object_ind)%sum_blk_double_2, dim1, dim2)
!    call alloc ('objects(object_ind)%previous_double_2',objects(object_ind)%previous_double_2, dim1, dim2)
   endif

  endif

  end subroutine object_alloc_double_3

!===========================================================================
  subroutine object_alloc_double_4 (object_name, object, dim1, dim2, dim3, dim4)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1, dim2, dim3, dim4

! input/output
  real(dp), allocatable               :: object (:,:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_double_4'
  integer all_err
  integer object_dim1, dim_min1
  integer object_dim2, dim_min2
  integer object_dim3, dim_min3
  integer object_dim4, dim_min4
  real(dp), allocatable                :: object_temp (:,:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1,dim2,dim3,dim4), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,4i8)') trim(lhere),': dimensions are ', dim1, dim2, dim3, dim4
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:,:,:,:) = 0.d0
   call object_associate (object_name, object, dim1, dim2, dim3, dim4)

   else
!  resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)
   object_dim4 = size(object,4)
   if (object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3 .or. object_dim4 /= dim4) then
    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)
    dim_min4 =  min(object_dim4, dim4)
    call alloc ('object_temp', object_temp, dim_min1, dim_min2, dim_min3, dim_min4)
    object_temp(:,:,:,:) = object(1:dim_min1,1:dim_min2,1:dim_min3,1:dim_min4)
    call object_release (object_name, object)
    allocate (object(dim1,dim2,dim3,dim4))
    object(:,:,:,:) = 0.d0
    call object_associate (object_name, object, dim1, dim2, dim3, dim4)
    object(1:dim_min1,1:dim_min2,1:dim_min3,1:dim_min4) = object_temp(:,:,:,:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_double_4

!===========================================================================
  subroutine object_alloc_double_5 (object_name, object, dim1, dim2, dim3, dim4, dim5)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 26 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)         :: object_name
  integer, intent(in)                  :: dim1, dim2, dim3, dim4, dim5

! input/output
  real(dp), allocatable, intent(inout) :: object (:,:,:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_double_5'
  integer all_err
  integer object_dim1, dim_min1
  integer object_dim2, dim_min2
  integer object_dim3, dim_min3
  integer object_dim4, dim_min4
  integer object_dim5, dim_min5
  real(dp), allocatable                :: object_temp (:,:,:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1,dim2,dim3,dim4,dim5), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,5i8)') trim(lhere),': dimensions are ', dim1, dim2, dim3, dim4, dim5
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:,:,:,:,:) = 0.d0
   call object_associate (object_name, object, dim1, dim2, dim3, dim4, dim5)

!  resize object if already allocated with different dimension
   else
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)
   object_dim4 = size(object,4)
   object_dim5 = size(object,5)
   if (object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3 .or. object_dim4 /= dim4 .or. object_dim5 /= dim5) then
    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)
    dim_min4 =  min(object_dim4, dim4)
    dim_min5 =  min(object_dim5, dim5)
    call alloc ('object_temp', object_temp, dim_min1, dim_min2, dim_min3, dim_min4, dim_min5)
    object_temp(:,:,:,:,:) = object(1:dim_min1,1:dim_min2,1:dim_min3,1:dim_min4,1:dim_min5)
    call object_release (object_name, object)
    allocate (object(dim1,dim2,dim3,dim4,dim5))
    object(:,:,:,:,:) = 0.d0
    call object_associate (object_name, object, dim1, dim2, dim3, dim4, dim5)
    object(1:dim_min1,1:dim_min2,1:dim_min3,1:dim_min4,1:dim_min5) = object_temp(:,:,:,:,:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_double_5

!===========================================================================
  subroutine object_alloc_double_row_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1

! input/output
  type (type_real_row), allocatable, intent(inout)  :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_double_row_1'
  integer all_err
  integer object_dim, dim_min
  type (type_real_row), allocatable    :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,i8)') trim(lhere),': dimensions are ', dim1
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   call object_associate (object_name, object, dim1)

! resize object if already allocated with different dimension
   else
   object_dim = size(object)
   if (object_dim /= dim1) then
    dim_min =  min(object_dim, dim1)
    call alloc ('object_temp', object_temp, dim_min)
    object_temp(:) = object(1:dim_min)
    call object_release (object_name, object)
    allocate (object(dim1))
    call object_associate (object_name, object, dim1)
    object(1:dim_min) = object_temp(:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_double_row_1

!===========================================================================
  subroutine object_alloc_complex_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 28 Sep 2013
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)         :: object_name
  integer, intent(in)                  :: dim1

! input/output
  complex(dpc), allocatable, intent(inout) :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_complex_1'
  integer all_err, object_dim, dim_min, object_ind
  complex(dpc), allocatable                :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,i8)') trim(lhere),': dimension is ', dim1
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:) = 0.d0
   call object_associate (object_name, object, dim1)

   else
!  resize object if already allocated with different dimension
   object_dim = size(object)
   if (object_dim /= dim1) then
    dim_min =  min(object_dim, dim1)
    call alloc ('object_temp', object_temp, dim_min)
    object_temp(:) = object(1:dim_min)
    call object_release (object_name, object)
    allocate (object(dim1))
    object(:) = 0.d0
    call object_associate (object_name, object, dim1)
    object(1:dim_min) = object_temp(:)
    call release ('object_temp', object_temp)
!   resize also the corresponding objects for averages and errors
    object_ind = object_index (object_name)
    call alloc ('objects(object_ind)%sum_complex_1',objects(object_ind)%sum_complex_1, dim1)
    call alloc ('objects(object_ind)%sum_blk_complex_1',objects(object_ind)%sum_blk_complex_1, dim1)
!    call alloc ('objects(object_ind)%previous_complex_1',objects(object_ind)%previous_complex_1, dim1)
   endif

  endif

  end subroutine object_alloc_complex_1

!===========================================================================
  subroutine object_alloc_logical_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1

! input/output
  logical, allocatable, intent(inout) :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_logical_1'
  integer all_err
  integer object_dim, dim_min
  logical , allocatable          :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,i8)') trim(lhere),': dimensions are ', dim1
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:) = .false.
   call object_associate (object_name, object, dim1)

!  resize object if already allocated with different dimension
   else
   object_dim = size(object)
   if (object_dim /= dim1) then
    dim_min =  min(object_dim, dim1)
    call alloc ('object_temp', object_temp, dim_min)
    object_temp(:) = object(1:dim_min)
    call object_release (object_name, object)
    allocate (object(dim1))
    object(:) = .false.
    call object_associate (object_name, object, dim1)
    object(1:dim_min) = object_temp(:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_logical_1

!===========================================================================
  subroutine object_alloc_logical_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1, dim2

! input/output
  logical, allocatable, intent(inout) :: object (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_logical_2'
  integer all_err
  integer object_dim1, dim_min1
  integer object_dim2, dim_min2
  logical, allocatable                :: object_temp (:,:)

! begin

! allocate object if not already allocated
  if (.not. allocated(object)) then
   allocate (object(dim1,dim2), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,2i8)') trim(lhere),': dimensions are ', dim1, dim2
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:,:) = .false.
   call object_associate (object_name, object, dim1, dim2)

!  resize object if already allocated with different dimension
   else
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   if (object_dim1 /= dim1 .or. object_dim2 /= dim2) then
    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    call alloc ('object_temp', object_temp, dim_min1, dim_min2)
    object_temp(:,:) = object(1:dim_min1,1:dim_min2)
    call object_release (object_name, object)
    allocate (object(dim1,dim2))
    object(:,:) = .false.
    call object_associate (object_name, object, dim1, dim2)
    object(1:dim_min1,1:dim_min2) = object_temp(:,:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_logical_2

!===========================================================================
  subroutine object_alloc_logical_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1, dim2, dim3

! input/output
  logical, allocatable, intent(inout) :: object (:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_logical_3'
  integer all_err
  integer object_dim1, dim_min1
  integer object_dim2, dim_min2
  integer object_dim3, dim_min3
  logical, allocatable                :: object_temp (:,:,:)

! begin

! allocate object if not already allocated
  if (.not. allocated(object)) then
   allocate (object(dim1,dim2,dim3), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,3i8)') trim(lhere),': dimensions are ', dim1, dim2, dim3
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   object(:,:,:) = .false.
   call object_associate (object_name, object, dim1, dim2, dim3)

!  resize object if already allocated with different dimension
   else
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)
   if (object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3) then
    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)
    call alloc ('object_temp', object_temp, dim_min1, dim_min2, dim_min3)
    object_temp(:,:,:) = object(1:dim_min1,1:dim_min2,1:dim_min3)
    call object_release (object_name, object)
    allocate (object(dim1,dim2,dim3))
    call object_associate (object_name, object, dim1, dim2, dim3)
    object(1:dim_min1,1:dim_min2,1:dim_min3) = object_temp(:,:,:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_logical_3

!===========================================================================
  subroutine object_alloc_string_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 31 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1

! input/output
  character(len=max_string_len), allocatable, intent(inout)  :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_string_1'
  integer all_err
  integer object_dim, dim_min
  character(len=max_string_len), allocatable       :: object_temp (:)

! begin

! allocate object if not already allocated
  if (.not. allocated(object)) then
   allocate (object(dim1), stat = all_err)
   if (all_err /= 0) then
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   call object_associate (object_name, object, dim1)

!  resize object if already allocated with different dimension
   else
   object_dim = size(object)
   if (object_dim /= dim1) then
    dim_min =  min(object_dim, dim1)
    call alloc ('object_temp', object_temp, dim_min)
    object_temp(:) = object(1:dim_min)
    call object_release (object_name, object)
    allocate (object(dim1))
    call object_associate (object_name, object, dim1)
    object(1:dim_min) = object_temp(:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_string_1

!===========================================================================
  subroutine object_alloc_string_row_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate an object and associate its name with its address
! Description : or resize it if already allocated
!
! Created     : J. Toulouse, 11 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1

! input/output
  type (type_string_row), allocatable, intent(inout)  :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'object_alloc_string_row_1'
  integer object_dim, dim_min, all_err
  type (type_string_row), allocatable    :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1), stat = all_err)
   if (all_err /= 0) then
    write(6,'(2a,i8)') trim(lhere),': dimensions are ', dim1
    call die (lhere,'allocation of object >'+trim(object_name)+'< failed.')
   endif
   call object_associate (object_name, object, dim1)

! resize object if already allocated with different dimension
   else
   object_dim = size(object)
   if (object_dim /= dim1) then
    dim_min =  min(object_dim, dim1)
    call alloc ('object_temp', object_temp, dim_min)
    object_temp(:) = object(1:dim_min)
    call object_release (object_name, object)
    allocate (object(dim1))
    call object_associate (object_name, object, dim1)
    object(1:dim_min) = object_temp(:)
    call release ('object_temp', object_temp)
   endif

  endif

  end subroutine object_alloc_string_row_1

!===========================================================================
  subroutine object_release_integer_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  integer, allocatable, intent(inout)     :: object(:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_integer_1

!===========================================================================
  subroutine object_release_integer_2 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  integer, allocatable, intent(inout)     :: object(:,:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_integer_2

!===========================================================================
  subroutine object_release_integer_3 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  integer, allocatable, intent(inout)     :: object(:,:,:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_integer_3

!===========================================================================
  subroutine object_release_integer_row_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                            :: object_name

! input/output
  type (type_integer_row), allocatable, intent(inout)     :: object(:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_integer_row_1

!===========================================================================
  subroutine object_release_double_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  real(dp), allocatable, intent(inout)    :: object(:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_double_1

!===========================================================================
  subroutine object_release_double_2 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  real(dp), allocatable, intent(inout)    :: object(:,:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_double_2

!===========================================================================
  subroutine object_release_double_3 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  real(dp), allocatable, intent(inout)    :: object(:,:,:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_double_3

!===========================================================================
  subroutine object_release_double_4 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  real(dp), allocatable, intent(inout)    :: object(:,:,:,:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_double_4

!===========================================================================
  subroutine object_release_double_5 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 26 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  real(dp), allocatable, intent(inout)    :: object(:,:,:,:,:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_double_5

!===========================================================================
  subroutine object_release_double_row_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                         :: object_name

! input/output
  type (type_real_row), allocatable, intent(inout)     :: object(:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_double_row_1

!===========================================================================
  subroutine object_release_complex_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 28 Sep 2013
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  complex(dpc), allocatable, intent(inout)    :: object(:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_complex_1

!===========================================================================
  subroutine object_release_logical_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  logical, allocatable, intent(inout)     :: object(:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_logical_1

!===========================================================================
  subroutine object_release_logical_2 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/ouput
  logical, allocatable, intent(inout)     :: object(:,:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_logical_2

!===========================================================================
  subroutine object_release_logical_3 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 01 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)            :: object_name

! input/output
  logical, allocatable, intent(inout)     :: object(:,:,:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_logical_3

!===========================================================================
  subroutine object_release_string_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 31 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                                 :: object_name

! input/output
  character(len=max_string_len), allocatable, intent(inout)    :: object(:)

! local
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_string_1

!===========================================================================
  subroutine object_release_string_row_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate and deassociate an object
!
! Created     : J. Toulouse, 11 Dec 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                         :: object_name

! input/output
  type (type_string_row), allocatable, intent(inout)     :: object(:)

! begin
  call release (object_name, object)
  call object_deassociate (object_name)

  end subroutine object_release_string_row_1

! ===================================================================================
  subroutine object_write (routine_name, object_name)
! -----------------------------------------------------------------------------------
! Description   : write an object
!
! Created       : J. Toulouse, 18 Oct 2005
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*) routine_name, object_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_write'
  character(len=max_string_len_type) :: object_type
  integer object_ind
  integer dim1, dim2
  integer i, j

! begin

! index of object
  object_ind = object_index_or_die (object_name)

! type of object
  object_type = objects(object_ind)%type
  if (object_type == '') then
    write(6,'(4a)') trim(lhere),': type of object ',trim(object_name),' is unknown'
    write(6,'(2a)') trim(lhere),': object has probably not been allocated with object_alloc'
   call die (lhere)
  endif

  select case (trim(object_type))
  case ('double_0')
     write(6,'(4a,es15.8)') trim(routine_name),': ',trim(object_name),'=' ,objects(object_ind)%pointer_double_0

  case ('double_1')
   dim1 = objects(object_ind)%dimensions(1)
   do i = 1, dim1
     write(6,'(4a,i3,a,es15.8)') trim(routine_name),': ',trim(object_name),'(',i,')=' ,objects(object_ind)%pointer_double_1(i)
   enddo

  case('double_2')
   dim1 = objects(object_ind)%dimensions(1)
   dim2 = objects(object_ind)%dimensions(2)
   do i = 1, dim1
    do j = 1, dim2
      write(6,'(4a,i3,a,i3,a,es15.8)') trim(routine_name),': ',trim(object_name),'(',i,',',j,')=' ,objects(object_ind)%pointer_double_2(i,j)
    enddo
   enddo

   case default
     call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
  end select

  end subroutine object_write

! ===================================================================================
  subroutine object_write_no_routine_name (object_name)
! -----------------------------------------------------------------------------------
! Description   : write an object
!
! Created       : J. Toulouse, 18 Oct 2005
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  integer object_ind

! begin

! index of object
  object_ind = object_index_or_die (object_name)

  call object_write_by_index (object_ind)

  end subroutine object_write_no_routine_name

! ===================================================================================
  subroutine object_write_2_no_routine_name (object_name1, object_name2)
! -----------------------------------------------------------------------------------
! Description   : write two objects
!
! Created       : J. Toulouse, 06 Dec 2005
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*) object_name1, object_name2

! local
  character(len=max_string_len_rout), save :: lhere = 'object_write_2_no_routine_name'
  character(len=max_string_len_type) :: object_type1, object_type2
  integer object_ind1, object_ind2
  integer dim11, dim21, dim31
  integer dim12, dim22, dim32
  integer i, j, k

! begin

! index of object
  object_ind1 = object_index (object_name1)
  if (object_ind1 == 0) then
   call die (lhere, 'object '+trim(object_name1)+' is not catalogued.')
  endif

! index of object
  object_ind2 = object_index (object_name2)
  if (object_ind2 == 0) then
   call die (lhere, 'object '+trim(object_name2)+' is not catalogued.')
  endif

! types of objects
  object_type1 = objects(object_ind1)%type
  object_type2 = objects(object_ind2)%type

  if (object_type1 == '') then
    write(6,'(4a)') trim(lhere),': type of object ',trim(object_name1),' is unknown'
    write(6,'(2a)') trim(lhere),': object has probably not been allocated with object_alloc'
   call die (lhere)
  endif
  if (object_type2 == '') then
    write(6,'(4a)') trim(lhere),': type of object ',trim(object_name2),' is unknown'
    write(6,'(2a)') trim(lhere),': object has probably not been allocated with object_alloc'
   call die (lhere)
  endif
  if (object_type1 /= object_type2) then
    write(6,'(5a)') trim(lhere),': type of object ',trim(object_name1),' is ', trim(object_type1)
    write(6,'(5a)') trim(lhere),': type of object ',trim(object_name2),' is ', trim(object_type2)
    write(6,'(2a)') trim(lhere),': they should be the same.'
   call die (lhere)
  endif

  select case (trim(object_type1))
   case ('double_0')
   write(6,'(2a,es15.8,3a,es15.8)') trim(object_name1),'=' ,objects(object_ind1)%pointer_double_0, ' ', trim(object_name2),'=' ,objects(object_ind2)%pointer_double_0

   case ('double_1')
   dim11 = objects(object_ind1)%dimensions(1)
   dim12 = objects(object_ind2)%dimensions(1)
   if (dim11 /= dim12) then
    write(6,'(4a,i8)') trim(lhere),': dimensions of object ',trim(object_name1),' is', dim11
    write(6,'(4a,i8)') trim(lhere),': dimensions of object ',trim(object_name2),' is', dim12
    write(6,'(2a)') trim(lhere),': they should be the same.'
    call die (lhere)
   endif
   do i = 1, dim11
     write(6,'(2a,i4,a,es15.8,3a,i4,a,es15.8)') trim(object_name1),'(',i,')=' ,objects(object_ind1)%pointer_double_1(i) , ' ',trim(object_name2),'(',i,')=' ,objects(object_ind2)%pointer_double_1(i)
   enddo

   case ('double_2')
   dim11 = objects(object_ind1)%dimensions(1)
   dim21 = objects(object_ind1)%dimensions(2)
   dim12 = objects(object_ind2)%dimensions(1)
   dim22 = objects(object_ind2)%dimensions(2)
   if (dim11 /= dim12 .or. dim21 /= dim22) then
    write(6,'(4a,2i8)') trim(lhere),': dimensions of object ',trim(object_name1),' are', dim11, dim21
    write(6,'(4a,2i8)') trim(lhere),': dimensions of object ',trim(object_name2),' are', dim12, dim22
    write(6,'(2a)') trim(lhere),': they should be the same.'
    call die (lhere)
   endif
   do i = 1, dim11
    do j = 1, dim21
      write(6,'(2a,i4,a,i4,a,es15.8,3a,i4,a,i4,a,es15.8)') trim(object_name1),'(',i,',',j,')=' ,objects(object_ind1)%pointer_double_2(i,j), ' ',trim(object_name2),'(',i,',',j,')=' ,objects(object_ind2)%pointer_double_1(i)
    enddo
   enddo

   case ('double_3')
   dim11 = objects(object_ind1)%dimensions(1)
   dim21 = objects(object_ind1)%dimensions(2)
   dim31 = objects(object_ind1)%dimensions(3)
   dim12 = objects(object_ind2)%dimensions(1)
   dim22 = objects(object_ind2)%dimensions(2)
   dim32 = objects(object_ind2)%dimensions(3)
   if (dim11 /= dim12 .or. dim21 /= dim22 .or. dim31 /= dim32) then
    write(6,'(4a,3i8)') trim(lhere),': dimensions of object ',trim(object_name1),' are', dim11, dim21, dim31
    write(6,'(4a,3i8)') trim(lhere),': dimensions of object ',trim(object_name2),' are', dim12, dim22, dim32
    write(6,'(2a)') trim(lhere),': they should be the same.'
    call die (lhere)
   endif
   do i = 1, dim11
    do j = 1, dim21
     do k = 1, dim31
      write(6,'(2a,i4,a,i4,a,i4,a,es15.8,3a,i4,a,i4,a,i4,a,es15.8)') trim(object_name1),'(',i,',',j,',',k,')=' ,objects(object_ind1)%pointer_double_3(i,j,k), ' ',trim(object_name2),'(',i,',',j,',',k,')=' ,objects(object_ind2)%pointer_double_3(i,j,k)
     enddo
    enddo
   enddo

   case ('integer_2')
   dim11 = objects(object_ind1)%dimensions(1)
   dim21 = objects(object_ind1)%dimensions(2)
   dim12 = objects(object_ind2)%dimensions(1)
   dim22 = objects(object_ind2)%dimensions(2)
   if (dim11 /= dim12 .or. dim21 /= dim22) then
    write(6,'(4a,2i8)') trim(lhere),': dimensions of object ',trim(object_name1),' are', dim11, dim21
    write(6,'(4a,2i8)') trim(lhere),': dimensions of object ',trim(object_name2),' are', dim12, dim22
    write(6,'(2a)') trim(lhere),': they should be the same.'
    call die (lhere)
   endif
   do i = 1, dim11
    do j = 1, dim21
      write(6,'(2a,i3,a,i3,a,i8,2a,i3,a,i3,a,i8)') trim(object_name1),'(',i,',',j,')=' ,objects(object_ind1)%pointer_integer_2(i,j), trim(object_name2),'(',i,',',j,')=' ,objects(object_ind2)%pointer_integer_2(i,j)
    enddo
   enddo

   case default
     call die (lhere, 'object type >'+trim(object_type1)+'< not handled.')
   end select

  end subroutine object_write_2_no_routine_name

! ===================================================================================
  subroutine object_write_by_index (object_ind)
! -----------------------------------------------------------------------------------
! Description   : write an object
!
! Created       : J. Toulouse, 18 Oct 2005
! -----------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) ::  object_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_write_by_index'
  character(len=max_string_len_obj) :: object_name
  character(len=max_string_len_type) :: object_type
  integer i, j, k, l, m

! begin

! object name
  object_name = objects(object_ind)%name

! object type
  object_type = objects(object_ind)%type

  select case (trim(object_type))
  case ('double_0')
     write(6,'(2a,es15.8)') trim(object_name),'=' ,objects(object_ind)%pointer_double_0

  case ('double_1')
   do i = 1, objects(object_ind)%dimensions(1)
     write(6,'(2a,i3,a,es15.8)') trim(object_name),'(',i,')=' ,objects(object_ind)%pointer_double_1(i)
   enddo

  case ('double_2')
   do j = 1, objects(object_ind)%dimensions(2)
    do i = 1, objects(object_ind)%dimensions(1)
      write(6,'(2a,i3,a,i3,a,es15.8)') trim(object_name),'(',i,',',j,')=' ,objects(object_ind)%pointer_double_2(i,j)
    enddo
   enddo

  case ('double_3')
   do k = 1, objects(object_ind)%dimensions(3)
    do j = 1, objects(object_ind)%dimensions(2)
     do i = 1, objects(object_ind)%dimensions(1)
      write(6,'(2a,i3,a,i3,a,i3,a,es15.8)') trim(object_name),'(',i,',',j,',',k,')=' ,objects(object_ind)%pointer_double_3(i,j,k)
     enddo
    enddo
   enddo

  case ('double_4')
   do l = 1, objects(object_ind)%dimensions(4)
    do k = 1, objects(object_ind)%dimensions(3)
     do j = 1, objects(object_ind)%dimensions(2)
      do i = 1, objects(object_ind)%dimensions(1)
      write(6,'(2a,i3,a,i3,a,i3,a,i3,a,es15.8)') trim(object_name),'(',i,',',j,',',k,',',l,')=' ,objects(object_ind)%pointer_double_4(i,j,k,l)
      enddo
     enddo
    enddo
   enddo

  case ('double_5')
   do m = 1, objects(object_ind)%dimensions(5)
    do l = 1, objects(object_ind)%dimensions(4)
     do k = 1, objects(object_ind)%dimensions(3)
      do j = 1, objects(object_ind)%dimensions(2)
       do i = 1, objects(object_ind)%dimensions(1)
       write(6,'(2a,i3,a,i3,a,i3,a,i3,a,i3,a,es15.8)') trim(object_name),'(',i,',',j,',',k,',',l,',',m,')=' ,objects(object_ind)%pointer_double_5(i,j,k,l,m)
       enddo
      enddo
     enddo
    enddo
   enddo

  case ('integer_0')
   write(6,'(2a,i3)') trim(object_name),'=' ,objects(object_ind)%pointer_integer_0

  case ('integer_1')
   do i = 1, objects(object_ind)%dimensions(1)
     write(6,'(2a,i3,a,i3)') trim(object_name),'(',i,')=' ,objects(object_ind)%pointer_integer_1(i)
   enddo

  case ('integer_2')
   do j = 1, objects(object_ind)%dimensions(2)
    do i = 1, objects(object_ind)%dimensions(1)
      write(6,'(2a,i3,a,i3,a,i8)') trim(object_name),'(',i,',',j,')=' ,objects(object_ind)%pointer_integer_2(i,j)
    enddo
   enddo

  case ('integer_3')
   do k = 1, objects(object_ind)%dimensions(3)
    do j = 1, objects(object_ind)%dimensions(2)
     do i = 1, objects(object_ind)%dimensions(1)
      write(6,'(2a,i3,a,i3,a,i3,a,i8)') trim(object_name),'(',i,',',j,',',k,')=' ,objects(object_ind)%pointer_integer_3(i,j,k)
     enddo
    enddo
   enddo

  case ('logical_2')
   do j = 1, objects(object_ind)%dimensions(2)
    do i = 1, objects(object_ind)%dimensions(1)
      write(6,'(2a,i3,a,i3,a,l)') trim(object_name),'(',i,',',j,')=' ,objects(object_ind)%pointer_logical_2(i,j)
    enddo
   enddo

  case ('logical_0') !fp
        write(6,'(2a,l)') trim(object_name),'=' ,objects(object_ind)%pointer_logical_0 !fp
     
! skip writing of objects of remaining types for now
  case ('string_1')
     write (6,'(5a)') 'Warning: skip writing object >',trim(object_name),'< of type >',trim(object_type),'<'
  case ('logical_1')
     write (6,'(5a)') 'Warning: skip writing object >',trim(object_name),'< of type >',trim(object_type),'<'
  case ('integer_row_1')
     write (6,'(5a)') 'Warning: skip writing object >',trim(object_name),'< of type >',trim(object_type),'<'
  case ('double_row_1')
     write (6,'(5a)') 'Warning: skip writing object >',trim(object_name),'< of type >',trim(object_type),'<'
  case ('undefined')
     write (6,'(5a)') 'Warning: skip writing object >',trim(object_name),'< of type >',trim(object_type),'<'

  case default
     call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
  end select

  end subroutine object_write_by_index

! ==============================================================================
  subroutine object_save (object_name)
! ------------------------------------------------------------------------------
! Description   : save current value of an object by its name
!
! Created       : J. Toulouse, 17 Feb 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  integer object_ind

! begin

! index of object
  object_ind = object_index_or_die (object_name)

  call object_save_by_index (object_ind)

  end subroutine object_save

! ==============================================================================
  recursive subroutine object_save_by_index (object_index)
! ------------------------------------------------------------------------------
! Description   : save current value of an object by its index
!
! Created       : J. Toulouse, 17 Feb 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_index

! local
  character(len=max_string_len_rout), save :: lhere = 'object_save_by_index'
  character(len=max_string_len_type) :: object_type
  integer dim1, dim2, dim3
!  integer obj_i, obj_j

! begin

! provide object if necessary
  call object_provide_by_index (object_index)

! check is object is associated
  call object_associated_or_die_by_index (object_index)

! type of object
  object_type = objects(object_index)%type

  select case (trim(object_type))
  case ('integer_0')
    objects(object_index)%save_integer_0 = objects(object_index)%pointer_integer_0

  case ('double_0')
    objects(object_index)%save_double_0 = objects(object_index)%pointer_double_0

  case ('double_1')
    dim1 = objects(object_index)%dimensions(1)
    call alloc ('objects(object_index)%save_double_1', objects(object_index)%save_double_1, dim1)
    objects(object_index)%save_double_1 (:) = objects(object_index)%pointer_double_1 (:)

  case ('double_2')
    dim1 = objects(object_index)%dimensions(1)
    dim2 = objects(object_index)%dimensions(2)
    call alloc ('objects(object_index)%save_double_2', objects(object_index)%save_double_2, dim1, dim2)
    objects(object_index)%save_double_2 (:,:) = objects(object_index)%pointer_double_2 (:,:)

  case ('double_3')
    dim1 = objects(object_index)%dimensions(1)
    dim2 = objects(object_index)%dimensions(2)
    dim3 = objects(object_index)%dimensions(3)
    call alloc ('objects(object_index)%save_double_3', objects(object_index)%save_double_3, dim1, dim2, dim3)
    objects(object_index)%save_double_3 (:,:,:) = objects(object_index)%pointer_double_3 (:,:,:)

  case default
     call die (lhere, 'type >'+trim(object_type)+'< of object >'+trim(objects(object_index)%name)+'< is unknown.')
  end select

! object has been saved
  objects(object_index)%saved = .true.

! save all the created objects by needed nodes
!  if (objects(object_index)%node_create_index > 0) then
!  do obj_i = 1, size(nodes(objects(object_index)%node_create_index)%objects_needed_index)
!   if (objects(nodes(objects(object_index)%node_create_index)%objects_needed_index(obj_i))%node_create_index > 0 ) then
!   do obj_j = 1, size(nodes(objects(nodes(objects(object_index)%node_create_index)%objects_needed_index(obj_i))%node_create_index)%objects_create_index)
!     call  object_save_by_index (nodes(objects(nodes(objects(object_index)%node_create_index)%objects_needed_index(obj_i))%node_create_index)%objects_create_index(obj_j))
!   enddo
!   endif
!  enddo
!  endif

  end subroutine object_save_by_index

! ==============================================================================
  subroutine object_restore (object_name)
! ------------------------------------------------------------------------------
! Description   : restore previously-saved value of an object by its index
!
! Created       : J. Toulouse, 17 Feb 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  integer object_ind

! begin

! index of object
  object_ind = object_index_or_die (object_name)

  call object_restore_by_index (object_ind)

  end subroutine object_restore

! ==============================================================================
  subroutine object_restore_by_index (object_index)
! ------------------------------------------------------------------------------
! Description   : restore previously-saved value of an object by its index
! Description   : and all the needed objects
!
! Created       : J. Toulouse, 17 Feb 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_index

! local
  character(len=max_string_len_rout), save :: lhere = 'object_restore_by_index'
  character(len=max_string_len_type) :: object_type
!  integer dim1, dim2, obj_i, obj_j

! begin
!  write (6,'(3a)') 'restoring object >',trim(objects(object_index)%name),'<.'

! begin by restore all the created objects by needed nodes
!  if (objects(object_index)%node_create_index > 0) then
!  do obj_i = 1, size(nodes(objects(object_index)%node_create_index)%objects_needed_index)
!   if (objects(nodes(objects(object_index)%node_create_index)%objects_needed_index(obj_i))%node_create_index > 0) then
!   do obj_j = 1, size(nodes(objects(nodes(objects(object_index)%node_create_index)%objects_needed_index(obj_i))%node_create_index)%objects_create_index)
!     call  object_restore_by_index (nodes(objects(nodes(objects(object_index)%node_create_index)%objects_needed_index(obj_i))%node_create_index)%objects_create_index(obj_j))
!    enddo
!    endif
!  enddo
!  endif

! check is object has been saved
  if (.not. objects(object_index)%saved) then
     call die (lhere, 'object >'+trim(objects(object_index)%name)+'< has never been saved.')
  endif

! type of object
  object_type = objects(object_index)%type

  select case (trim(object_type))
  case ('integer_0')
    objects(object_index)%pointer_integer_0 = objects(object_index)%save_integer_0

  case ('double_0')
    objects(object_index)%pointer_double_0 = objects(object_index)%save_double_0

  case ('double_1')
    objects(object_index)%pointer_double_1 (:) = objects(object_index)%save_double_1 (:)

  case ('double_2')
    objects(object_index)%pointer_double_2 (:,:) = objects(object_index)%save_double_2 (:,:)

  case ('double_3')
    objects(object_index)%pointer_double_3 (:,:,:) = objects(object_index)%save_double_3 (:,:,:)

  case default
    call die (lhere, 'type >'+trim(object_type)+'< of object >'+trim(objects(object_index)%name)+'< is unknown.')
  end select

! object has been modified
  call object_modified2_by_index (object_index)

  end subroutine object_restore_by_index

! ==============================================================================
  subroutine object_freeze (object_name)
! ------------------------------------------------------------------------------
! Description   : freeze an object (keep it always valid)
!
! Created       : J. Toulouse, 19 Apr 2007
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  integer object_ind

! begin

! index of object
  object_ind = object_index_or_die (object_name)

  objects(object_ind)%freezed = .true.

  end subroutine object_freeze

! ==============================================================================
  subroutine object_zero (object_name)
! ------------------------------------------------------------------------------
! Description   : zero an object
!
! Created       : J. Toulouse, 19 Apr 2007
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_zero'
  character(len=max_string_len_type) :: object_type
  integer object_ind

! begin

! index of object
  object_ind = object_index_or_die (object_name)

! die if object not associated
  call object_associated_or_die_by_index (object_ind)

! type of object
  object_type = objects(object_ind)%type

  select case (trim(object_type))
  case ('integer_0')
    objects(object_ind)%pointer_integer_0 = 0

  case ('double_0')
    objects(object_ind)%pointer_double_0 = 0.d0

  case ('double_1')
    objects(object_ind)%pointer_double_1 (:) = 0.d0

  case ('double_2')
    objects(object_ind)%pointer_double_2 (:,:) = 0.d0

  case ('double_3')
    objects(object_ind)%pointer_double_3 (:,:,:) = 0.d0

  case default
    call die (lhere, 'type >'+trim(object_type)+'< of object >'+trim(objects(object_ind)%name)+'< not handled.')
  end select

  call object_modified_by_index (object_ind)

  end subroutine object_zero

!! ==============================================================================
! subroutine object_not_NaN_or_die_by_index (object_index)
!! ------------------------------------------------------------------------------
!! Description   : check if object is not a NaN, die otherwise
!!
!! Created       : J. Toulouse, 19 Apr 2008
!! ------------------------------------------------------------------------------
!  implicit none
!
!! input
!  integer, intent(in) :: object_index
!
!! local
!  character(len=max_string_len_rout), save :: lhere = 'object_not_NaN_or_die_by_index'
!  character(len=max_string_len_type) :: object_type
!  logical l
!
! begin
!
!! type of object
!  object_type = objects(object_index)%type
!
!  select case (trim(object_type))
!  case ('integer_0')
!    if (objects(object_index)%pointer_integer_0 == 'NaN') then
!     write (6,'(a)') 'debug: object >'+objects(object_index)%name+'< is NaN.'
!    endif
!
!  case ('integer_1')
!    where (objects(object_index)%pointer_integer_1 == 'NaN')
!!     write (6,'(a)') 'debug: object >'+objects(object_index)%name+'< is NaN.'
!     l=.true.
!    endwhere
!
!  case ('integer_2')
!    where (objects(object_index)%pointer_integer_2 == 'NaN')
!     write (6,'(a)') 'debug: object >'+objects(object_index)%name+'< is NaN.'
!    endwhere
!
!  case ('double_0')
!    if (objects(object_index)%pointer_double_0 == 'NaN') then
!     write (6,'(a)') 'debug: object >'+objects(object_index)%name+'< is NaN.'
!    endif
!
!  case ('double_1')
!    where (objects(object_index)%pointer_double_1 == 'NaN')
!     write (6,'(a)') 'debug: object >'+objects(object_index)%name+'< is NaN.'
!    endwhere
!
!  case ('double_2')
!    where (objects(object_index)%pointer_double_2 == 'NaN')
!     write (6,'(a)') 'debug: object >'+objects(object_index)%name+'< is NaN.'
!    endwhere
!
!  case ('double_3')
!    where (objects(object_index)%pointer_double_3 == 'NaN')
!     write (6,'(a)') 'debug: object >'+objects(object_index)%name+'< is NaN.'
!    endwhere
!
!  case default
!  end select
!
!  end subroutine object_not_NaN_or_die_by_index

end module objects_mod

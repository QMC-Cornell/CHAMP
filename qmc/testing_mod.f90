module testing_mod

  use all_tools_mod
  real(dp), allocatable         :: aa(:)
  real(dp), pointer :: pa(:)

  contains

  subroutine testing

  implicit none
  character*(max_string_len_rout), save :: lhere= 'test'

  real(dp), allocatable         :: mat_test(:,:)
  real(dp), allocatable         :: mat_test_inv(:,:)
  integer i, j

  call alloc ('mat_test', mat_test, 3, 3)
  call alloc ('mat_test_inv', mat_test_inv, 3, 3)
  mat_test (1,1) = 1.d0 
  mat_test (1,2) = 2.d0 
  mat_test (1,3) = 3.d0 
  mat_test (2,1) = -3.d0 
  mat_test (2,2) = -5.d0 
  mat_test (2,3) = 0.5d0 
  mat_test (3,1) = 0.5d0 
  mat_test (3,2) = 0.5d0 
  mat_test (3,3) = 2.d0 
  
  do i =1,3 
  write(6,*) "mat_test=",(mat_test(i,j),j=1,3)
  enddo

  call inverse_by_svd (mat_test, mat_test_inv, 3, 1.d-10)
  do i =1,3 
  write(6,*) "mat_test_inv=",(mat_test_inv(i,j),j=1,3)
  enddo

  return
  call alloc_test ('aa',aa,3)


  aa(1)=3.d0
  aa(2)=1.d0
  aa(3)=2.d0
  write(6,*) trim(lhere),': aa=',aa

  write(6,*) trim(lhere),': reshape(aa,(/2/))=', reshape(aa,(/2/))
  write(6,*) trim(lhere),': reshape(aa,(/4/))=', reshape(aa,(/4/))
  stop

  aa(1:0) =  -55555555555555555555555553.d0
  write(6,*) trim(lhere),': aa(1:0)=',aa(1:0)
  write(6,*) trim(lhere),': aa(1:3)=',aa(1:3)


  stop

  write(6,*) trim(lhere),': aa=',aa(:)


!  pa => aa

  write(6,*) trim(lhere),': a=',aa(:)
  write(6,*) trim(lhere),': pa=',pa(:)

  pa(2)=-6.d0

  write(6,*) trim(lhere),': a=',aa(:)
  write(6,*) trim(lhere),': pa=',pa(:)

  write(6,*) trim(lhere),': end of test'
  stop

end subroutine testing

!===========================================================================
  subroutine alloc_test (object_name, object, dim1)
!---------------------------------------------------------------------------
!  Description : Allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1

! output
  real(dp) , allocatable,target          :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_double_1'
  integer all_err
  integer i, object_dim, dim_min
  real(dp) , allocatable          :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1), stat = all_err)

   pa => object

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die(lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim = size(object)

   if(object_dim /= dim1) then

    dim_min =  min(object_dim, dim1)

    call alloc ('object_temp', object_temp, dim_min)


    do i = 1, dim_min
     object_temp(i) = object(i)
    enddo

    call release (object_name, object)

    call alloc (object_name, object, dim1)

    do i = 1, dim_min
     object(i) = object_temp(i)
    enddo

    call release ('object_temp', object_temp)

   endif

  endif

  end subroutine alloc_test

end module testing_mod

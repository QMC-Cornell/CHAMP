module matrix_tools_mod

  use basic_tools_mod
  use strings_tools_mod

  contains

! ==============================================================================
  subroutine inverse_by_svd (matrix, matrix_inv, dim, threshold)
! ------------------------------------------------------------------------------
! Description   : Calculate inverse of square matrix by SVD with a threshold on singular values
!
! Created       : J. Toulouse, 04 Nov 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent (in)      :: matrix (:,:)
  integer,          intent (in)      :: dim
  real(dp), intent (in)      :: threshold

! output
  real(dp), intent (out)      :: matrix_inv (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'inverse_by_svd'
  integer i, j, k
  real(dp), allocatable :: mat_u(:,:)
  real(dp), allocatable :: mat_v(:,:)
  real(dp), allocatable :: mat_w(:)
  real(dp), allocatable :: mat_w_inv(:)
  integer w_kept_nb
  real(dp), allocatable :: mat_a(:,:)
  real(dp), allocatable :: mat_vt(:,:)
  real(dp), allocatable :: work (:)
  integer lwork
  integer info

! begin
!  lhere = 'inverse_by_svd'
!  write(6,*) trim(lhere),': entering'
  
! temporary arrays for SVD
  call alloc ('mat_u', mat_u, dim, dim)
  call alloc ('mat_v', mat_v, dim, dim)
  call alloc ('mat_w', mat_w, dim)
  call alloc ('mat_w_inv', mat_w_inv, dim)

  lwork = 10 * dim
  call alloc ('mat_a', mat_a, dim, dim)
  call alloc ('mat_vt', mat_vt, dim, dim)
  call alloc ('mat_w', mat_w, dim)
  call alloc ('work', work, lwork)

! SVD from numerical recipes
!  mat_u = matrix ! matrix to be inverted
!  call svdcmp (mat_u, dim, dim, dim, dim, mat_w, mat_v)
!    
!  write(6,*) trim(lhere),': SVD from numerical recipes:'
!  write(6,*) trim(lhere),': mat_u=',mat_u
!  write(6,*) trim(lhere),': mat_v=',mat_v
!  write(6,*) trim(lhere),': mat_w=',mat_w

! SVD from Lapack
  mat_a = matrix ! matrix to be inverted
!  write(6,*) trim(lhere),': mat_a=',mat_a
!  write(6,*) trim(lhere),': before dgesvd'
  call dgesvd( 'A', 'A', dim, dim, mat_a, dim, mat_w, mat_u, dim, mat_vt, dim, work, lwork, info)
!  write(6,*) trim(lhere),': after dgesvd'
  if (info /= 0) then
   call die (lhere, 'problem in dgesvd')
  endif

 mat_v = transpose(mat_vt)

!  write(6,*) trim(lhere),': SVD from Lapack:'
!  write(6,*) trim(lhere),': mat_u=',mat_u
!  write(6,*) trim(lhere),': mat_v=',mat_v
!  write(6,*) trim(lhere),': mat_w=',mat_w

! Singular values
!JT  do i = 1, dim
!JT   write(6,*) trim(lhere), ': i=',i,' mat_w=',mat_w(i)
!JT  enddo 

! Inverse singular values  (drop small singular values)
!JT  write(6,*) trim(lhere), ': threshold on singular values: ',threshold
  w_kept_nb = dim
  do i = 1, dim
   if (mat_w (i) < threshold) then
     mat_w_inv (i) = 0.d0
     w_kept_nb = w_kept_nb - 1
   else
     mat_w_inv (i) = 1.d0 /  mat_w (i)
   endif
  enddo
!JT  write(6,*) trim(lhere), ': number of singular values dropped:', dim - w_kept_nb

! calculate inverse matrix
  matrix_inv = 0.d0
  
   do i = 1, dim
     do j = 1, dim
       do k = 1, dim
         matrix_inv (i, j) = matrix_inv (i, j) + mat_v (i, k) * mat_w_inv (k) * mat_u (j, k)
       enddo
     enddo
   enddo
 
!  write(6,*) trim(lhere), ': matrix_inv=',matrix_inv
  
! release arrays for SVD
  call release ('mat_u', mat_u)
  call release ('mat_v', mat_v)
  call release ('mat_w', mat_w)
  call release ('mat_w_inv', mat_w_inv)
  call release ('mat_a', mat_a)
  call release ('mat_vt', mat_vt)
  call release ('work', work)

!  write(6,*) trim(lhere), ': exiting'
  end subroutine inverse_by_svd

! ==============================================================================
  subroutine eigensystem (matrix, eigenvectors, eigenvalues, dim)
! ------------------------------------------------------------------------------
! Description   : compute eigenvectors and eigenvalues of a real symmetrix matrix
! Description   : using Lapack routine
!
! Created       : J. Toulouse, 11 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)   :: matrix (:,:)
  integer,  intent(in)   :: dim

! output
  real(dp), intent(out)  :: eigenvectors (:,:)
  real(dp), intent(out)  :: eigenvalues (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'eigensystem'
  real(dp), allocatable :: mat_a(:,:)
  real(dp), allocatable :: work (:)
  real(dp) :: matrix_check
  integer lwork, info, i, j, k

! begin
  lwork = 10 * dim
  call alloc ('mat_a', mat_a, dim, dim)
  call alloc ('work', work, lwork)

  mat_a (:,:) = matrix (:,:)

  call dsyev('V','U',dim, mat_a, dim, eigenvalues, work, lwork, info)
  if (info /= 0) then
   call die (lhere, 'exiting dsyev with info='+info+' /= 0.')
  endif

  eigenvectors (:,:) = mat_a (:,:)

! checkings
  call is_a_number_or_die ('eigenvalues', eigenvalues)

! check recovery of original matrix after diagonalization
  do i = 1, dim
   do j = 1, dim
     matrix_check = 0.d0
     do k = 1, dim
      matrix_check = matrix_check + eigenvectors (i, k) * eigenvalues (k) * eigenvectors (j, k)
     enddo ! k
     if (abs(matrix_check-matrix(i,j)) > 1.d-8) then
        call die (lhere, 'low accuracy in diagonalization; the error on a matrix element is '+abs(matrix_check-matrix(i,j)))
     endif
   enddo ! j
  enddo ! i

! release arrays
  call release ('mat_a', mat_a)
  call release ('work', work)

  end subroutine eigensystem

! ==============================================================================
  real(dp) function trace (matrix)
! ------------------------------------------------------------------------------
! Description   : returns trace of a square matrix 
!
! Created       : J. Toulouse, 20 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  real(dp)             , intent (in)   :: matrix (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'trace'
  integer dim1, dim2
  integer i

! begin
  dim1 = size(matrix,1)
  dim2 = size(matrix,2)
  
  if (dim1 /= dim2) then
   write(6,*) trim(lhere),': in routine', trim(here),' trace requested of a non-square matrix!'
   write(6,*) trim(lhere),': dim1=',dim1,' /= dim2=',dim2
   call die (lhere)
  endif

  trace = 0.d0
  do i = 1, dim1
   trace = trace + matrix (i,i) 
  enddo

  end function trace

! ==============================================================================
  integer function kronecker_delta (i, j)
! ------------------------------------------------------------------------------
! Description   : kronecker_delta
!
! Created       : J. Toulouse, 19 Mar 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent (in)   :: i, j

! begin
  if (i == j) then
   kronecker_delta = 1
  else
   kronecker_delta = 0
  endif 

  end function kronecker_delta

end module matrix_tools_mod

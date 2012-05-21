module matrix_tools_mod

  use basic_tools_mod
  use strings_tools_mod
  use variables_mod


!===============================================================
  interface inverse_by_svd
!---------------------------------------------------------------
   module procedure inverse_by_svd_double , &
                    inverse_by_svd_complex
  end interface inverse_by_svd

  contains

! ==============================================================================
  subroutine inverse_by_svd_double (matrix, matrix_inv, dim, threshold)
! ------------------------------------------------------------------------------
! Description   : Calculate inverse of square matrix by SVD with a threshold on singular values
!
! Created       : J. Toulouse, 04 Nov 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent (in)      :: matrix (:,:)
  integer,  intent (in)      :: dim
  real(dp), intent (in)      :: threshold

! output
  real(dp), intent (out)     :: matrix_inv (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'inverse_by_svd_double'
  integer i, j, k, w_kept_nb, lwork, info
  real(dp), allocatable :: mat_a(:,:), mat_u(:,:), mat_vt(:,:) 
  real(dp), allocatable :: mat_w(:), work (:)

  if (dim == 0) return

! temporary arrays for SVD
  call alloc ('mat_a', mat_a, dim, dim)
  call alloc ('mat_u', mat_u, dim, dim)
  call alloc ('mat_vt', mat_vt, dim, dim)
  call alloc ('mat_w', mat_w, dim)
  mat_a = matrix

! calculate optimal value of lwork
  call alloc ('work', work, 1)
  call dgesvd ('A', 'A', dim, dim, mat_a, dim, mat_w, mat_u, dim, mat_vt, dim, work, -1, info)
  if(info /= 0) then
   call die (lhere, 'problem in dgesvd (while calculating optimal value of lwork): info='+info+' /= 0')
  endif
  lwork =  work(1)
  call alloc ('work', work, lwork)

! SVD decomposition: A = U . W . V^T
  call dgesvd ('A', 'A', dim, dim, mat_a, dim, mat_w, mat_u, dim, mat_vt, dim, work, lwork, info)
  call release ('mat_a', mat_a)
  call release ('work', work)
  if (info /= 0) then
   call die (lhere, 'problem in dgesvd: info='+info+' /= 0')
  endif

! calculate inverse matrix: A^-1 = V . W^-1 . U^T
  matrix_inv = 0.d0
  w_kept_nb = 0
  do k = 1, dim
    if (mat_w (k) >= threshold) then
      w_kept_nb = w_kept_nb + 1
      do i = 1, dim
        do j = 1, dim
          matrix_inv (i, j) = matrix_inv (i, j) + mat_vt (k, i) * (1.d0 / mat_w (k)) * mat_u (j, k)
        enddo
      enddo
   endif
  enddo
!  write(6,*) trim(lhere), ': number of singular values dropped:', dim - w_kept_nb

!  write(6,*) trim(lhere), ': matrix_inv=',matrix_inv

  call release ('mat_u', mat_u)
  call release ('mat_vt', mat_vt)
  call release ('mat_w', mat_w)

  end subroutine inverse_by_svd_double

! ==============================================================================
  subroutine inverse_by_svd_complex (matrix, matrix_inv, dim, threshold)
! ------------------------------------------------------------------------------
! Description   : Calculate inverse of square matrix by SVD with a threshold on singular values
!
! Created       : J. Toulouse, 02 May 2010
! ------------------------------------------------------------------------------
  implicit none

! input
  double complex, intent (in)  :: matrix (:,:)
  integer,  intent (in)        :: dim
  real(dp), intent (in)        :: threshold

! output
  double complex, intent (out) :: matrix_inv (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'inverse_by_svd_complex'
  integer i, j, k, w_kept_nb, lwork, info
  double complex, allocatable :: mat_a(:,:), mat_u(:,:), mat_vh(:,:), work(:)
  real(dp), allocatable :: mat_w(:), rwork (:)

  if (dim == 0) return

! temporary arrays for SVD
  call alloc ('mat_a', mat_a, dim, dim)
  call alloc ('mat_u', mat_u, dim, dim)
  call alloc ('mat_vh', mat_vh, dim, dim)
  call alloc ('mat_w', mat_w, dim)
  call alloc ('rwork', rwork, 5*dim)
  mat_a = matrix

! calculate optimal value of lwork
  call alloc ('work', work, 1)
  call die (lhere, 'need zgesvd')
!  call zgesvd ('A', 'A', dim, dim, mat_a, dim, mat_w, mat_u, dim, mat_vh, dim, work, -1, rwork, info)
  if(info /= 0) then
   call die (lhere, 'problem in zgesvd (while calculating optimal value of lwork): info='+info+' /= 0')
  endif
  lwork =  work(1)
  call alloc ('work', work, lwork)

! SVD decomposition: A = U . W . V^H
  call die (lhere, 'need zgesvd')
!  call zgesvd ('A', 'A', dim, dim, mat_a, dim, mat_w, mat_u, dim, mat_vh, dim, work, lwork, rwork, info)
  call release ('mat_a', mat_a)
  call release ('work', work)
  call release ('rwork', rwork)
  if (info /= 0) then
   call die (lhere, 'problem in zgesvd: info='+info+' /= 0')
  endif

! calculate inverse matrix: A^-1 = V . W^-1 . U^H
  matrix_inv = 0.d0
  w_kept_nb = 0
  do k = 1, dim
    if (mat_w (k) >= threshold) then
      w_kept_nb = w_kept_nb + 1
      do i = 1, dim
        do j = 1, dim
          matrix_inv (i, j) = matrix_inv (i, j) + dconjg(mat_vh (k, i)) * (1.d0 / mat_w (k)) * dconjg(mat_u (j, k))
        enddo
      enddo
   endif
  enddo
!  write(6,*) trim(lhere), ': number of singular values dropped:', dim - w_kept_nb

!  write(6,*) trim(lhere), ': matrix_inv=',matrix_inv

  call release ('mat_u', mat_u)
  call release ('mat_vh', mat_vh)
  call release ('mat_w', mat_w)

  end subroutine inverse_by_svd_complex

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
  integer, save :: warnings_nb = 0
  integer       :: warnings_nb_max = 20

! begin
  if (dim == 0) return

  lwork = 10 * dim
  call alloc ('mat_a', mat_a, dim, dim)
  call alloc ('work', work, lwork)

  mat_a (:,:) = matrix (:,:)

  call dsyev ('V','U',dim, mat_a, dim, eigenvalues, work, lwork, info)
  call release ('work', work)
  if (info /= 0) then
   call die (lhere, 'exiting dsyev with info='+info+' /= 0.')
  endif

  eigenvectors (:,:) = mat_a (:,:)
  call release ('mat_a', mat_a)

! checkings
!  call is_a_number_or_die ('eigenvalues', eigenvalues)

! check recovery of original matrix after diagonalization
  do i = 1, dim
   do j = 1, dim
     matrix_check = 0.d0
     do k = 1, dim
      matrix_check = matrix_check + eigenvectors (i, k) * eigenvalues (k) * eigenvectors (j, k)
     enddo ! k
     if(warnings_nb < warnings_nb_max .and. abs(matrix_check-matrix(i,j)) > 1.d-7) then
       write(6,'(''Warning: low accuracy in diagonalization; the error on a matrix element'',2i4,'' is'',d12.4)') &
         i,j,matrix_check-matrix(i,j)
       l_warning = .true.
       warnings_nb = warnings_nb + 1
       if (warnings_nb == warnings_nb_max) then
        write(6,'(a)') 'all further similar warnings will be suppressed'
       endif
     endif
! JT: Warning: comment out this stop for now
!     if(abs(matrix_check-matrix(i,j)) > 1.d-2) then
!       call die (lhere, 'low accuracy in diagonalization; the error on a matrix element is '+abs(matrix_check-matrix(i,j)))
!     endif
   enddo ! j
  enddo ! i

  end subroutine eigensystem

! ==============================================================================
  subroutine to_the_power (matrix, dim, power_n, matrix_out)
! ------------------------------------------------------------------------------
! Description   : returns the matrix to the power n
! Description   : valid for a real symmetric matrix (uses eigensystem)
!
! Created       : B. Mussard, 09 Mar 2010
! ------------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)   :: matrix (:,:)
  integer,  intent(in)   :: dim
  real(dp), intent(in)   :: power_n

! output
  real(dp), intent(out)  :: matrix_out (:,:)

! local
!  character(len=max_string_len_rout), save :: lhere = 'to_the_power'
  integer eigen_i, bas_i, bas_j, bas_k
  real(dp), allocatable  :: eigenvectors (:,:)
  real(dp), allocatable  :: eigenvalues (:)

! begin
  call alloc('eigenvectors',eigenvectors,dim,dim)
  call alloc('eigenvalues',eigenvalues,dim)

  call eigensystem (matrix, eigenvectors, eigenvalues, dim)
 
  if (power_n < 1) then
     do eigen_i = 1, dim
        if (eigenvalues (eigen_i) < -1.d-5) then
            write(6,'(a)') 'one eigenvalue is too negative to be equalled to zero'
        else if (eigenvalues (eigen_i) < 0) then
            eigenvalues (eigen_i) = 0
        endif
     enddo
  endif
            

  do bas_i = 1, dim
    do bas_j = 1, dim
       matrix_out(bas_i, bas_j) = 0.d0
       do bas_k = 1, dim
          matrix_out(bas_i, bas_j) = matrix_out(bas_i, bas_j) + eigenvectors (bas_i, bas_k) * eigenvalues(bas_k)**power_n * eigenvectors (bas_j, bas_k)
       enddo ! bas_k
    enddo ! bas_j
  enddo ! bas_i

  end subroutine to_the_power

! ==============================================================================
  real(dp) function matrix_determinant (matrix, n)
! ------------------------------------------------------------------------------
! Description   : returns determinant of a matrix
! WARNING: not checked!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Function to find the determinant of a square matrix
! Author : Louisda16th a.k.a Ashwith J. Rego
! Description: The subroutine is based on two key points:
! 1] A determinant is unaltered when row operations are performed: Hence, using this principle,
! row operations (column operations would work as well) are used
! to convert the matrix into upper traingular form
! 2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
! Created       : J. Toulouse, 29 Oct 2009
! ------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(dp), DIMENSION(n,n) :: matrix
  INTEGER, INTENT(IN) :: n
  REAL(dp) :: m, temp
  INTEGER :: i, j, k, l
  LOGICAL :: DetExists = .TRUE.

  l = 1
  !Convert to upper triangular form
  DO k = 1, n-1
          IF (matrix(k,k) == 0) THEN
                  DetExists = .FALSE.
                  DO i = k+1, n
                          IF (matrix(i,k) /= 0) THEN
                                  DO j = 1, n
                                          temp = matrix(i,j)
                                          matrix(i,j)= matrix(k,j)
                                          matrix(k,j) = temp
                                  END DO
                                  DetExists = .TRUE.
                                  l=-l
                                  EXIT
                          ENDIF
                  END DO
                  IF (DetExists .EQV. .FALSE.) THEN
                          matrix_determinant = 0
                          return
                  END IF
          ENDIF
          DO j = k+1, n
                  m = matrix(j,k)/matrix(k,k)
                  DO i = k+1, n
                          matrix(j,i) = matrix(j,i) - m*matrix(k,i)
                  END DO
          END DO
  END DO
  
  !Calculate determinant by finding product of diagonal elements
  matrix_determinant = l
  DO i = 1, n
          matrix_determinant = matrix_determinant * matrix(i,i)
  END DO
             
  end function matrix_determinant

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

module basis_mod

  use all_tools_mod
  use control_mod

! Declaration of global variables and default values
  real(dp), allocatable           :: phin_norm (:,:)
  real(dp), allocatable           :: phin_ortho (:,:)
  integer, allocatable            :: basis_fns_cent (:)
  integer, allocatable            :: basis_fns_by_center_type_nb (:)
  type (type_string_row), allocatable      :: basis_fns_type_by_center_type (:)
  type (type_real_row), allocatable        :: basis_fns_expo_by_center_type (:)
  type (type_integer_row), allocatable     :: basis_fns_rad_by_center_type (:)
  character(len=max_string_len), allocatable :: basis_fns_type (:)
  character(len=max_string_len), allocatable :: basis_fns_name (:)
  real(dp), allocatable           :: basis_ovlp (:,:)
  real(dp), allocatable           :: basis_ovlp_eigvec (:,:)
  real(dp), allocatable           :: basis_ovlp_eigval (:)
  real(dp), allocatable           :: basis_ovlp_12 (:,:)
  real(dp), allocatable           :: basis_ovlp_m12 (:,:)
  real(dp), allocatable           :: zex_sav (:)
  real(dp), allocatable           :: zex2_sav (:,:)
  real(dp), allocatable           :: zex_best (:)
  real(dp), allocatable           :: norm_basis (:)
  integer                         :: exp_opt_lab_read_nb
  integer, allocatable            :: exp_opt_lab_read (:)
  character(len=max_string_len)   :: basis_functions_varied = 'normalized'
  logical                         :: l_optimize_log_exp = .false.
  character(len=max_string_len)   :: which_analytical_basis = 'none' !fp
  logical                         :: l_purely_analytical_basis = .true.
  integer, parameter              :: max_n_gauss_slat = 5 !fp
  real(dp), dimension(max_n_gauss_slat), parameter   :: norm_gauss_slat_exp_1 = (/ 1.12646742161049d0, 0.5766099503612371d0, 0.1965811411216233d0, 0.05027565586963441d0, 0.01028077216468808d0 /)  !fp

  real(dp), external              :: gamma1

  contains

!===========================================================================
  subroutine basis_menu
!---------------------------------------------------------------------------
! Description : menu for basis
!
! Created     : J. Toulouse, 18 Apr 2007
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'basis_menu'
  real(dp), allocatable  :: exponents (:)
  integer exponents_nb, bas_i, bas_c_i, cent_i, cent_type_i, ict, m

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of basis menu ----------------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for basis menu:'
   write(6,'(a)') 'basis'
   write(6,'(a)') '  ibasis = [integer]: type of basis (default: 1)'
   write(6,'(a)') '  which_analytical_basis = [none|slater|gaussian|gauss-slater] (default: slater for all-electron, none otherwise)'
   write(6,'(a)') '  numr = [integer] options for numerical basis (default: -3)'
   write(6,'(a)') '  basis_functions : function types and exponents'
   write(6,'(a)') '   1  ! atom type'
   write(6,'(a)') '   1S   9.5799650 '
   write(6,'(a)') '   2PX  4.6785620'
   write(6,'(a)') '   2PY  4.6785620'
   write(6,'(a)') '   2PZ  4.6785620'
   write(6,'(a)') '   2  ! atom type'
   write(6,'(a)') '   1S   2.4731060 '
   write(6,'(a)') '  end'
   write(6,'(a)') '  exponents 0.2 0.1 end (not used anymore)'
   write(6,'(a)') '  basis_functions_varied = [unnormalized|normalized|orthonormalized] : choice of basis functions for exponent optimization (default=normalized)'
   write(6,'(a)') '  optimize_log_exp = [bool] : optimize logarithm of exponents (default=false)'
   write(6,'(a)') '  optimized_exponents 1 2 3 4 end : list of labels of exponents to optimize (default: all exponents)'
   write(6,'(a)') 'end'
   write(6,*)

  case ('ibasis')
   call get_next_value (ibasis)

  case ('numr')
   call get_next_value (numr)
   call object_modified ('numr')

  case ('basis_functions')
   if (numr <=0) then
    call basis_functions_slater
   else
    call basis_functions_num
   endif

  case ('exponents') ! not used anymore
   call get_next_value_list ('exponents', exponents, exponents_nb)
   call object_provide ('nbasis')
   write(6,*) 'nbasis=',nbasis
   write(6,*) 'exponents_nb=',exponents_nb
   call require (lhere, 'exponents_nb == nbasis', exponents_nb == nbasis)
   write(6,*)
   write(6,'(a,500f12.6)')  ' exponents: ', exponents (:)
   zex (1:nbasis,1) = exponents (:)
   call object_modified ('zex')
   call copy_zex_zex2
!   write(6,'(a,500f12.6)')  'zex2: ', zex2


  case ('basis_functions_varied')
   call get_next_value (basis_functions_varied)

  case ('which_analytical_basis')  !fp
   call get_next_value (which_analytical_basis) !fp
   call object_modified('which_analytical_basis') !fp

  case ('optimize_log_exp')
   call get_next_value (l_optimize_log_exp)

  case ('optimized_exponents')
   call get_next_value_list_object ('exp_opt_lab_read', exp_opt_lab_read, exp_opt_lab_read_nb)
   call object_provide ('nbasis')
   if (exp_opt_lab_read_nb > nbasis) then
    call die (lhere, 'number of optimized exponents read = '+exp_opt_lab_read_nb+' > number of basis functions nbasis ='+nbasis)
   endif
   call object_modified ('exp_opt_lab_read_nb')

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  if (use_parser) then

  if(ibasis == 1) then
    write(6,'(a)') ' use localized basis functions (ibasis=1)'
    notype=0
   else
    call die (lhere, 'case ibasis='+ibasis+'not yet implemented for new input')
  endif

! default for all-electron calculations: analytical Slater functions
  call object_provide ('nloc')
  if (nloc ==0) then
   which_analytical_basis = 'slater'
  endif

  call object_provide ('zex')
  if(minval(zex(:,1)).ne.0.d0) then
    l_purely_analytical_basis = .true.
  else
    l_purely_analytical_basis = .false.
  endif

  select case(trim(which_analytical_basis))
  case ('none')
    write(6,'(a)') ' the basis functions are all numerical'
  case ('slater')
    write(6,'(a)') ' the analytical basis functions are Slater functions'
  case ('gaussian')
    write(6,'(a)') ' the analytical basis functions are Gaussian functions'
  case ('gauss-slater')
    write(6,'(a)') ' the analytical basis functions are Gauss-Slater functions'
  case default
   call die (lhere, 'unknown keyword which_analytical_basis='+trim(which_analytical_basis)+'')
  end select
  
  write(6,'(a,i5)') ' numr=',numr

!   if (numr>0) then
!      if (trim(which_analytical_basis) .eq. 'none') then !fp
!         write(6,'(a)') ' numerical basis functions used'
!      else
!         write(6,*) 'numr>0 but which_analytical_basis .ne. none'
!         stop 'numr>0 but which_analytical_basis .ne. none'
!      end if
!   end if

  call object_provide ('nbasis')
  write(6,'(/a,i5)') ' number of basis functions = ',nbasis

! distinct_radial_bas must be called to update n_bas2, etc...
  call distinct_radial_bas

! if(ibasis.eq.1.and.numr.gt.0.and.inum_orb.eq.0) call read_bas_num(1)
! if(ibasis.eq.1 .and. minval(zex(:,1)).eq.0.d0 .and. inum_orb.eq.0) call read_bas_num(1)
  if(ibasis.eq.1 .and. inum_orb.eq.0) call read_bas_num(1)

  call object_provide ('nctype')
  call object_provide ('basis_fns_by_center_type_nb')
  call object_provide ('basis_fns_type_by_center_type')

! analytic (Slater) basis
  write(6,*) 
  if (numr<=0) then
    call object_provide ('basis_fns_expo_by_center_type')
    do cent_type_i = 1, nctype
      write (6, '(a,i5)') ' center type # ',cent_type_i
      do bas_c_i = 1, basis_fns_by_center_type_nb (cent_type_i)
        write (6, '(a,i5,a,a6,a,f12.6)') ' function # ', bas_c_i, ': ',trim(basis_fns_type_by_center_type (cent_type_i)%row(bas_c_i)), '       exponent = ', basis_fns_expo_by_center_type (cent_type_i)%row(bas_c_i)
      enddo
    enddo

! numerical basis
  else
    call object_provide ('basis_fns_rad_by_center_type')
    do cent_type_i = 1, nctype
      write (6, '(a,i5)') ' center type # ',cent_type_i
      do bas_c_i = 1, basis_fns_by_center_type_nb (cent_type_i)
        write (6, '(a,i5,a,a6,a,i5)') ' function # ', bas_c_i, ': ',trim(basis_fns_type_by_center_type (cent_type_i)%row(bas_c_i)), '       radial function # ', basis_fns_rad_by_center_type (cent_type_i)%row(bas_c_i)
        if (basis_fns_rad_by_center_type (cent_type_i)%row(bas_c_i) > nrbas(cent_type_i)) then
          call die (lhere, ' for center type # '+cent_type_i+', radial function # '+basis_fns_rad_by_center_type (cent_type_i)%row(bas_c_i)+' > nrbas='+nrbas(cent_type_i))
        endif
      enddo
    enddo

  endif


  if(ibasis.eq.1) then
! irecursion_ylm=0 use Cyrus' spherical harmonics (upto g functions)
! irecursion_ylm=1 use Ryo' spherical harmonics (any L)
! Note that at present it always calculates upto lmax (set in basis_fns.f) and so it takes long if lmax is large.
! Change it to calculate upto largest l actually used.
        irecursion_ylm=0
!       irecursion_ylm=1
!       write(6,'(''Warning temporarily set irecursion_ylm=1'')')
!        call read_orb_loc
    if(irecursion_ylm.eq.0)then
      write(6,'(a)') ' not using recursion for spherical harmonics'
     elseif(irecursion_ylm.eq.1) then
      write(6,'(a)') ' using recursion for spherical harmonics'
      call setup_spherical_harmonics
      call setup_coefficients_ylm
     else
      stop 'irecursion_ylm must be 0 or 1'
    endif
  endif

! Check that irecursion_ylm=1 if l of basis function >=4
  do bas_i=1,nbasis
    if(l_bas(bas_i).ge.5 .and. irecursion_ylm.eq.0) then
      call die (lhere, ' basis functions with l>=5 are used, set irecursion_ylm=1')
    endif
  enddo  

  endif ! if use_parser

  write(6,'(a)') 'End of basis menu ----------------------------------------------------------------------------------------'

  end subroutine basis_menu

!===========================================================================
  subroutine basis_functions_slater
!---------------------------------------------------------------------------
! Description : read and set up Slater basis functions
!
! Created     : J. Toulouse, 03 Mar 2009
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'basis_functions_slater'
  character(len=500) line
  character(len=max_string_len) basis_fns_type_temp
  real(dp) basis_fns_expo_temp
  integer cent_type_i, bas_i, iostat,i, cent_i, bas_c_i

! begin
  call object_provide ('nctype')

  call object_alloc ('basis_fns_by_center_type_nb', basis_fns_by_center_type_nb, nctype)
  call object_alloc ('basis_fns_type_by_center_type', basis_fns_type_by_center_type, nctype)
  call object_alloc ('basis_fns_expo_by_center_type', basis_fns_expo_by_center_type, nctype)
  basis_fns_by_center_type_nb (:) = 0

  do
   read(unit_input,'(a)',iostat=iostat) line
   if(iostat < 0) then
     call die (lhere, 'error while reading basis functions')
   endif
!   write(6,*) 'line >',trim(line),'<'

!  convert to lower case
   call upplow (line)

!  skip if empty line
   if (trim(line) == '') cycle

!  exit when 'end' is read
   if (index(line,'end') /= 0) exit

!  new center type
   if (is_string_integer (line)) then
 
     cent_type_i = string_to_integer (line)

     if (cent_type_i > nctype) then
       call die (lhere, ' center type index = '+cent_type_i+' > number of center types ='+nctype)
     endif

     cycle
   endif

!  read basis function type and exponent
   basis_fns_by_center_type_nb (cent_type_i) = basis_fns_by_center_type_nb (cent_type_i) + 1
   
   read(line,*,iostat=iostat) basis_fns_type_temp , basis_fns_expo_temp
   if(iostat < 0) then
     call die (lhere, 'error while reading basis functions')
   endif
   
   
   call append (basis_fns_type_by_center_type (cent_type_i)%row, basis_fns_type_temp)
   call append (basis_fns_expo_by_center_type (cent_type_i)%row, basis_fns_expo_temp)
 
  enddo

  call alloc ('nbasis_ctype', nbasis_ctype, nctype)
  nbasis_ctype (1:nctype) = basis_fns_by_center_type_nb (1:nctype)
  mbasis_ctype = maxval (nbasis_ctype)
  call object_modified ('nbasis_ctype')
  call object_modified ('mbasis_ctype')
  call object_modified ('basis_fns_by_center_type_nb')
  call object_modified ('basis_fns_type_by_center_type')
  call object_modified ('basis_fns_expo_by_center_type')

! initialization
  bas_i = 0

  call object_provide ('ncent')
  call object_provide ('nwf')

  do cent_i = 1, ncent
    cent_type_i = iwctype (cent_i)

    do bas_c_i = 1, basis_fns_by_center_type_nb (cent_type_i)
      bas_i = bas_i + 1
      call alloc ('n_bas', n_bas, bas_i)
      call alloc ('l_bas', l_bas, bas_i)
      call alloc ('m_bas', m_bas, bas_i)

    select case (trim(basis_fns_type_by_center_type (cent_type_i)%row(bas_c_i)))
     case ('1s');   n_bas (bas_i) = 1; l_bas (bas_i) = 0; m_bas (bas_i) = 0
     case ('2s');   n_bas (bas_i) = 2; l_bas (bas_i) = 0; m_bas (bas_i) = 0
     case ('3s');   n_bas (bas_i) = 3; l_bas (bas_i) = 0; m_bas (bas_i) = 0
     case ('4s');   n_bas (bas_i) = 4; l_bas (bas_i) = 0; m_bas (bas_i) = 0
     case ('5s');   n_bas (bas_i) = 5; l_bas (bas_i) = 0; m_bas (bas_i) = 0
     case ('2px');  n_bas (bas_i) = 2; l_bas (bas_i) = 1; m_bas (bas_i) = 1
     case ('2py');  n_bas (bas_i) = 2; l_bas (bas_i) = 1; m_bas (bas_i) = -1
     case ('2pz');  n_bas (bas_i) = 2; l_bas (bas_i) = 1; m_bas (bas_i) = 0
     case ('3px');  n_bas (bas_i) = 3; l_bas (bas_i) = 1; m_bas (bas_i) = 1
     case ('3py');  n_bas (bas_i) = 3; l_bas (bas_i) = 1; m_bas (bas_i) = -1
     case ('3pz');  n_bas (bas_i) = 3; l_bas (bas_i) = 1; m_bas (bas_i) = 0
     case ('4px');  n_bas (bas_i) = 4; l_bas (bas_i) = 1; m_bas (bas_i) = 1
     case ('4py');  n_bas (bas_i) = 4; l_bas (bas_i) = 1; m_bas (bas_i) = -1
     case ('4pz');  n_bas (bas_i) = 4; l_bas (bas_i) = 1; m_bas (bas_i) = 0
     case ('5px');  n_bas (bas_i) = 5; l_bas (bas_i) = 1; m_bas (bas_i) = 1
     case ('5py');  n_bas (bas_i) = 5; l_bas (bas_i) = 1; m_bas (bas_i) = -1
     case ('5pz');  n_bas (bas_i) = 5; l_bas (bas_i) = 1; m_bas (bas_i) = 0
     case ('3d0');  n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = 0
     case ('3d+1'); n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = 1
     case ('3d-1'); n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = -1
     case ('3d+2'); n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = 2
     case ('3d-2'); n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = -2
     case ('4d0');  n_bas (bas_i) = 4; l_bas (bas_i) = 2; m_bas (bas_i) = 0
     case ('4d+1'); n_bas (bas_i) = 4; l_bas (bas_i) = 2; m_bas (bas_i) = 1
     case ('4d-1'); n_bas (bas_i) = 4; l_bas (bas_i) = 2; m_bas (bas_i) = -1
     case ('4d+2'); n_bas (bas_i) = 4; l_bas (bas_i) = 2; m_bas (bas_i) = 2
     case ('4d-2'); n_bas (bas_i) = 4; l_bas (bas_i) = 2; m_bas (bas_i) = -2
     case ('5d0');  n_bas (bas_i) = 5; l_bas (bas_i) = 2; m_bas (bas_i) = 0
     case ('5d+1'); n_bas (bas_i) = 5; l_bas (bas_i) = 2; m_bas (bas_i) = 1
     case ('5d-1'); n_bas (bas_i) = 5; l_bas (bas_i) = 2; m_bas (bas_i) = -1
     case ('5d+2'); n_bas (bas_i) = 5; l_bas (bas_i) = 2; m_bas (bas_i) = 2
     case ('5d-2'); n_bas (bas_i) = 5; l_bas (bas_i) = 2; m_bas (bas_i) = -2
     case ('4f0');  n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = 0
     case ('4f+1'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = 1
     case ('4f-1'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = -1
     case ('4f+2'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = 2
     case ('4f-2'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = -2
     case ('4f+3'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = 3
     case ('4f-3'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = -3
     case ('5f0');  n_bas (bas_i) = 5; l_bas (bas_i) = 3; m_bas (bas_i) = 0
     case ('5f+1'); n_bas (bas_i) = 5; l_bas (bas_i) = 3; m_bas (bas_i) = 1
     case ('5f-1'); n_bas (bas_i) = 5; l_bas (bas_i) = 3; m_bas (bas_i) = -1
     case ('5f+2'); n_bas (bas_i) = 5; l_bas (bas_i) = 3; m_bas (bas_i) = 2
     case ('5f-2'); n_bas (bas_i) = 5; l_bas (bas_i) = 3; m_bas (bas_i) = -2
     case ('5f+3'); n_bas (bas_i) = 5; l_bas (bas_i) = 3; m_bas (bas_i) = 3
     case ('5f-3'); n_bas (bas_i) = 5; l_bas (bas_i) = 3; m_bas (bas_i) = -3
     case ('5g0');  n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 0
     case ('5g+1'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 1
     case ('5g-1'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = -1
     case ('5g+2'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 2
     case ('5g-2'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = -2
     case ('5g+3'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 3
     case ('5g-3'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = -3
     case ('5g+4'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 4
     case ('5g-4'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = -4

     case default
       call die (lhere, 'unknown basis function label >'+trim(basis_fns_type_by_center_type (cent_type_i)%row(bas_c_i))+'<')
    end select


    call alloc ('zex', zex, bas_i, nwf)
    zex (bas_i,1) = basis_fns_expo_by_center_type (cent_type_i)%row(bas_c_i)
    call alloc ('ictype_basis', ictype_basis, bas_i)
    ictype_basis (bas_i)= cent_type_i

    enddo
  enddo

! for now, ML_BAS needs to be at least 4 because of the code in basis_fns.f
  ML_BAS = max (4, maxval(l_bas))

  nbasis = bas_i

  call object_modified ('nbasis')
  call object_modified ('n_bas')
  call object_modified ('l_bas')
  call object_modified ('m_bas')
  call object_modified ('ictype_basis')
  call object_modified ('zex')


  end subroutine basis_functions_slater

!===========================================================================
  subroutine basis_functions_num
!---------------------------------------------------------------------------
! Description : read and set up numerical basis functions
!
! Created     : J. Toulouse, 16 Jun 2009
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'basis_functions_num'
  character(len=500) line
  character(len=max_string_len) basis_fns_type_temp
  integer basis_fns_rad_temp
  integer cent_type_i, bas_i, iostat, cent_i, bas_c_i

! begin
  call object_provide ('nctype')

  call object_alloc ('basis_fns_by_center_type_nb', basis_fns_by_center_type_nb, nctype)
  call object_alloc ('basis_fns_type_by_center_type', basis_fns_type_by_center_type, nctype)
  call object_alloc ('basis_fns_rad_by_center_type', basis_fns_rad_by_center_type, nctype)
  basis_fns_by_center_type_nb (:) = 0

  do
   read(unit_input,'(a)',iostat=iostat) line
   if(iostat < 0) then
     call die (lhere, 'error while reading basis functions')
   endif
!   write(6,*) 'line >',trim(line),'<'

!  convert to lower case
   call upplow (line)

!  skip if empty line
   if (trim(line) == '') cycle

!  exit when 'end' is read
   if (index(line,'end') /= 0) exit

!  new center type
   if (is_string_integer (line)) then
 
     cent_type_i = string_to_integer (line)

     if (cent_type_i > nctype) then
       call die (lhere, ' center type index = '+cent_type_i+' > number of center types ='+nctype)
     endif

     cycle
   endif

!  read basis function type and radial function
   basis_fns_by_center_type_nb (cent_type_i) = basis_fns_by_center_type_nb (cent_type_i) + 1
   
   read(line,*,iostat=iostat) basis_fns_type_temp , basis_fns_rad_temp
   if(iostat < 0) then
     call die (lhere, 'error while reading basis functions')
   endif
   
   
   call append (basis_fns_type_by_center_type (cent_type_i)%row, basis_fns_type_temp)
   call append (basis_fns_rad_by_center_type (cent_type_i)%row, basis_fns_rad_temp)
 
  enddo

  call alloc ('nbasis_ctype', nbasis_ctype, nctype)
  nbasis_ctype (1:nctype) = basis_fns_by_center_type_nb (1:nctype)
  mbasis_ctype = maxval (nbasis_ctype)
  call object_modified ('nbasis_ctype')
  call object_modified ('mbasis_ctype')
  call object_modified ('basis_fns_by_center_type_nb')
  call object_modified ('basis_fns_type_by_center_type')
  call object_modified ('basis_fns_rad_by_center_type')

! initialization
  bas_i = 0

  call object_provide ('ncent')
  call object_provide ('nwf')

  do cent_i = 1, ncent
    cent_type_i = iwctype (cent_i)
    
    do bas_c_i = 1, basis_fns_by_center_type_nb (cent_type_i)
      bas_i = bas_i + 1
      call alloc ('n_bas', n_bas, bas_i)
      call alloc ('l_bas', l_bas, bas_i)
      call alloc ('m_bas', m_bas, bas_i)

    select case (trim(basis_fns_type_by_center_type (cent_type_i)%row(bas_c_i)))
     case ('s');   n_bas (bas_i) = 1; l_bas (bas_i) = 0; m_bas (bas_i) = 0
     case ('px');  n_bas (bas_i) = 2; l_bas (bas_i) = 1; m_bas (bas_i) = 1
     case ('py');  n_bas (bas_i) = 2; l_bas (bas_i) = 1; m_bas (bas_i) = -1
     case ('pz');  n_bas (bas_i) = 2; l_bas (bas_i) = 1; m_bas (bas_i) = 0
     case ('d0');  n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = 0
     case ('d+1'); n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = 1
     case ('d-1'); n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = -1
     case ('d+2'); n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = 2
     case ('d-2'); n_bas (bas_i) = 3; l_bas (bas_i) = 2; m_bas (bas_i) = -2
     case ('f0');  n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = 0
     case ('f+1'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = 1
     case ('f-1'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = -1
     case ('f+2'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = 2
     case ('f-2'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = -2
     case ('f+3'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = 3
     case ('f-3'); n_bas (bas_i) = 4; l_bas (bas_i) = 3; m_bas (bas_i) = -3
     case ('g0');  n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 0
     case ('g+1'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 1
     case ('g-1'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = -1
     case ('g+2'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 2
     case ('g-2'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = -2
     case ('g+3'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 3
     case ('g-3'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = -3
     case ('g+4'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = 4
     case ('g-4'); n_bas (bas_i) = 5; l_bas (bas_i) = 4; m_bas (bas_i) = -4
     case ('h0');  n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = 0
     case ('h+1'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = 1
     case ('h-1'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = -1
     case ('h+2'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = 2
     case ('h-2'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = -2
     case ('h+3'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = 3
     case ('h-3'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = -3
     case ('h+4'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = 4
     case ('h-4'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = -4
     case ('h+5'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = 5
     case ('h-5'); n_bas (bas_i) = 6; l_bas (bas_i) = 5; m_bas (bas_i) = -5

     case default
       call die (lhere, 'unknown basis function label >'+trim(basis_fns_type_by_center_type (cent_type_i)%row(bas_c_i))+'<')
    end select 


    call alloc ('iwrwf', iwrwf, mbasis_ctype ,nctype)
    iwrwf (bas_c_i, cent_type_i) = basis_fns_rad_by_center_type (cent_type_i)%row(bas_c_i)
    call alloc ('ictype_basis', ictype_basis, bas_i)
    ictype_basis (bas_i)= cent_type_i
!   for now, a numerical basis function is indicated by zero exponent:
    call alloc ('zex', zex, bas_i, nwf)
    zex (bas_i,1) = 0.d0

    enddo
  enddo
   
! for now, ML_BAS needs to be at least 4 because of the code in basis_fns.f
  ML_BAS = max (4, maxval(l_bas))

  nbasis = bas_i

  call object_modified ('nbasis')
  call object_modified ('n_bas')
  call object_modified ('l_bas')
  call object_modified ('m_bas')
  call object_modified ('ictype_basis')
  call object_modified ('iwrwf')
  call object_modified ('zex')

  end subroutine basis_functions_num

! ==============================================================================
  subroutine norm_basis_bld
! ------------------------------------------------------------------------------
! Description   : Normalization constants of basis functions
!
! Taken from Cyrus Umrigar
! Set normalization of basis fns.
! In 3d:
! Norm of radial part: ((2*zeta)^{2n+1}/(2n)!)^{1/2} for Slaters (n>0).
!                      (2*(2*zeta)^(n+1/2)/Gamma(n+1/2))^{1/2} for gaussians
! where                Gamma(n+1/2) = integral {t^{n+1/2} e^-t dt} = (2n-1)!! sqrt(pi)/2^n for gaussians
!                      Gamma(1/2)=sqrt(pi), Gamma(a+1)=a*Gamma(a), Gamma(a)=(a-1)!
! obtained by Integrate[r^2 (r^{n-1} Exp[-zeta r])^2,{r,0,Infinity}] for Slaters
! obtained by Integrate[r^2 (r^{n-1} Exp[-zeta r^2])^2,{r,0,Infinity}] for Gaussians
! Norm of angular part: ((2*l+1)/(4*pi))^{1/2}.
! In 2d:
! Norm of radial part: ((2*zeta)^{2n}/(2n-1)!)^{1/2}.
! Norm of angular part: (min(m+1,2)/(2*pi))^{1/2}.
! If numr =0 or -1 we are using analytic basis functions and we use normalization
!                  for angular and radial parts.  The -1 is just to tell it to order
!                  basis functions by all the s's first, then all the p's etc.
!                  instead of 1s,2s,2p,3s,3p,3d,...
!         =1       we are using numerical basis functions and we use normalization
!                  for angular part only
! Whether one is using Slater or gaussian basis fns. was inputted by having
! n1s,n2s etc. be either > 0 or < 0.  Now use which_analytical_basis
!
! Created       : J. Toulouse, 16 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'norm_basis_bld'
  integer ib, ict, n, l, m, n1

! header
  if (header_exe) then

   call object_create ('norm_basis')

   call object_needed ('nbasis')
!   call object_needed ('ibasis')    !   ADG
   call object_needed ('ndim')
   call object_needed ('numr')
   call object_needed ('n_bas')
   call object_needed ('l_bas')
   call object_needed ('m_bas')
   call object_needed ('iwf')

   return

  endif

! begin

! allocation
  call object_alloc ('norm_basis', norm_basis, nbasis)

! provide zex for analytic basis functions
! must use object_provide_in_node and not object_provide 
! because the dependency of norm_basis on zex must be recorded
  if(numr <= 0) then
   call object_provide_in_node (lhere, 'zex')
  endif

  do ib = 1, nbasis
     ict=ictype_basis(ib)
     n = n_bas(ib)
     if (ndim == 3) then
       l = l_bas(ib)
       if(iwrwf2(ib).le.nrbas_analytical(ict)) then
         select case (trim(which_analytical_basis))   !fp
          case ('slater')   !fp
             norm_basis(ib)=sqrt((2*zex(ib,iwf))**(2*n+1)*(2*l+1)/(factorial(2*n)*4*pi))   !fp
          case ('gaussian')   !fp
             n1=abs(n)   !fp
             norm_basis(ib)=sqrt(2*(2*zex(ib,iwf))**(n1+0.5d0)*(2*l+1)/(gamma1(n1)*4*pi))   !fp
          case ('gauss-slater')   !fp
             write(*,*) "Testing: calculating gauss-slater normalization" !fp
             norm_basis(ib) = norm_gauss_slat_exp_1(n) * sqrt((2*l+1)/(4*pi)) * zex(ib,iwf)**(n+0.5d0)   !fp
          case default
             write(6,*) 'norm_basis_bld: Allowed basis types are slater gaussian gauss-slater!'
             stop 'norm_basis_bld: Allowed basis types are slater gaussian gauss-slater!'
         end select   !fp
        else
         norm_basis(ib)=sqrt((2*l+1)/(4*pi))
       endif
     elseif (ndim == 2) then
       m = m_bas(ib)
!      ibasis=4  ! testing !ADG
       if (numr <= 0 .and. ibasis < 4) then !ADG
         norm_basis(ib)=sqrt((2*zex(ib,iwf))**(2*n)*min(abs(m)+1,2)/(factorial(2*n-1)*2*pi))
! The following change is not necessary at this time and has not been tested.
!          elseif(numr<=0 .and. (ibasis>=4 .and. ibasis<=6)) then
           elseif(iwrwf2(ib)<=nrbas_analytical(ict) .and. (ibasis>=4 .and. ibasis<=6)) then
         norm_basis(ib)=sqrt(1/pi) !ADG
       else
! Warning: temporarily commented out diff norm for m=0
!        norm_basis(ib)=sqrt(min(abs(m)+1,2)/(2*pi))
         norm_basis(ib)=sqrt(1/pi)
       endif
     endif
  enddo ! ib

  end subroutine norm_basis_bld

! ==============================================================================
  subroutine phin_norm_bld
! ------------------------------------------------------------------------------
! Description   : Normalized basis functions
!
! Created       : J. Toulouse, 17 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i

! header
  if (header_exe) then

   call object_create ('phin_norm')

   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('phin')
   call object_needed ('norm_basis')

   return

  endif

! begin

! allocation
  call object_alloc ('phin_norm', phin_norm, nelec, nbasis)

  do bas_i = 1, nbasis
    phin_norm (1:nelec, bas_i) = norm_basis (bas_i) * phin (bas_i, 1:nelec)
  enddo ! bas_i

  end subroutine phin_norm_bld

! ==============================================================================
  subroutine phin_ortho_bld
! ------------------------------------------------------------------------------
! Description   : Orthonormalized basis functions
!
! Created       : J. Toulouse, 30 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_k

! header
  if (header_exe) then

   call object_create ('phin_ortho')

   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('phin')
   call object_needed ('basis_ovlp_m12')

   return

  endif

! begin

! allocation
  call object_alloc ('phin_ortho', phin_ortho, nelec, nbasis)

  do bas_i = 1, nbasis
      phin_ortho (1:nelec, bas_i) = 0.d0
    do bas_k = 1, nbasis
      phin_ortho (1:nelec, bas_i) = phin_ortho (1:nelec, bas_i) + basis_ovlp_m12 (bas_i, bas_k) * phin (bas_k, 1:nelec)
    enddo ! bas_k
  enddo ! bas_i

  end subroutine phin_ortho_bld

! ==============================================================================
  subroutine basis_fns_cent_bld
! ------------------------------------------------------------------------------
! Description   : centers of each basis functions
!
! Created       : J. Toulouse, 25 Jan 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer ib, ic, ict, ib2

! header
  if (header_exe) then

   call object_create ('basis_fns_cent')

   call object_needed ('ncent')
   call object_needed ('iwctype')
   call object_needed ('nbasis_ctype')
   call object_needed ('nbasis')

   return

  endif

! begin

! requirements
  if (numr > 0) then
   call die (here, 'numr='+numr+'>0. Implemented only for the case analytical functions numr <= 0.')
  endif

! allocation
  call object_alloc ('basis_fns_cent', basis_fns_cent, nbasis)

  ib=0
  do ic=1,ncent
     ict=iwctype(ic)
     do ib2=1,nbasis_ctype(ict)
        ib=ib+1
        basis_fns_cent (ib) = ic
     enddo
  enddo

  end subroutine basis_fns_cent_bld

! ==============================================================================
  subroutine basis_fns_type_bld
! ------------------------------------------------------------------------------
! Description   : type of basis functions
!
! Created       : J. Toulouse, 21 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i

! header
  if (header_exe) then

   call object_create ('basis_fns_type')

   call object_needed ('ncent')
   call object_needed ('nbasis')
   call object_needed ('n_bas')
   call object_needed ('l_bas')
   call object_needed ('m_bas')
   call object_needed ('basis_fns_cent')

   return

  endif

! begin

! allocation
  call object_alloc ('basis_fns_type', basis_fns_type, nbasis)

  do bas_i = 1, nbasis

!    principal quantum number
     basis_fns_type (bas_i) = trim(string (abs(n_bas (bas_i))))

!    angular and magnetic quantum numbers
     select case (l_bas (bas_i))
     case (0)
       basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'s'
     case (1)
       select case (m_bas (bas_i))
       case (1)
         basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'px'
       case (-1)
         basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'py'
       case (0)
         basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'pz'
       case default
         call die (here, 'wrong magnetic quantum number m_bas='+m_bas (bas_i))
       end select
     case (2)
       select case (m_bas (bas_i))
       case (0)
         basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'dzr'
       case (2)
         basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'dx2y2'
       case (-2)
         basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'dxy'
       case (1)
         basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'dxz'
       case (-1)
         basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'dyz'
       case default
         call die (here, 'wrong magnetic quantum number m_bas='+m_bas (bas_i))
       end select
     case (3)
       basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'f'
     case (4)
       basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'g'
     case (5)
       basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'h'
     case default
       call die (here, 'add symbol for angular momentum l_bas ='+l_bas(bas_i)+'.')
     end select

!    add center index
     if (ncent > 1) then
       basis_fns_type (bas_i) = trim(basis_fns_type (bas_i))//'('//trim(string(basis_fns_cent (bas_i)))//')'
     endif

  enddo

!  do bas_i = 1, nbasis
!    write(6,'(a,i3,a,a)') 'bas_i=',bas_i,' basis_fns_type=',trim(basis_fns_type (bas_i))
!  enddo ! dexp_i

  end subroutine basis_fns_type_bld

! ==============================================================================
  subroutine basis_fns_name_bld
! ------------------------------------------------------------------------------
! Description   : name of basis functions
!
! Created       : J. Toulouse, 21 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i

! header
  if (header_exe) then

   call object_create ('basis_fns_name')

   call object_needed ('nbasis')
   call object_needed ('basis_fns_type')

   return

  endif

! begin

! allocation
  call object_alloc ('basis_fns_name', basis_fns_name, nbasis)

  do bas_i = 1, nbasis
     basis_fns_name (bas_i) = trim(string2(bas_i))//'_'//trim(basis_fns_type (bas_i))
  enddo

  end subroutine basis_fns_name_bld

! ==============================================================================
  function slater_ovlp (n1, l1, m1, exp1, n2, l2, m2, exp2) result(result)
! ------------------------------------------------------------------------------
! Description   : overlap matrix of unormalized Slater basis functions
! Description   : only for one center
! Description   : Int r^(n1-1) exp(-z1 r) Yl1,m1 r^(n2-1) exp(-z2 r) Yl2,m2 d3r
!
! Created       : J. Toulouse, 17 Oct 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in)  :: n1, l1, m1, n2, l2, m2
  real(dp), intent(in) :: exp1, exp2

! output
  real(dp) :: result

! begin
  result = 0.d0

! orthogonality of spherical harmonics
  if (l1 /= l2 .or. m1 /= m2) return

  result = 4.d0 * pi * factorial (n1 + n2) / ( (exp1 + exp2)**(1 + n1 + n2) ) &   ! radial contribution
           * 1.d0 / (sqrt(2.d0 * l1 + 1.d0) * sqrt(2.d0 * l2 + 1.d0))             ! spherical harmonics contribution

  return
  end function slater_ovlp

! ==============================================================================
  subroutine basis_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : build overlap matrix of unormalized Slater basis functions
! Description   : only for one center
! Description   : Int r^n1 exp(-z1 r) Yl1,m1 r^n2 exp(-z2 r) Yl2,m2 d3r
!
! Created       : J. Toulouse, 17 Oct 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_j
  integer n_i, n_j, l_i, l_j, m_i, m_j
  real(dp) exp_i, exp_j

! header
  if (header_exe) then

   call object_create ('basis_ovlp')

   call object_needed ('nbasis')
   call object_needed ('n_bas')
   call object_needed ('l_bas')
   call object_needed ('m_bas')
   call object_needed ('zex')

   return

  endif

! begin

! requirements
  if (numr /= 0 .and. numr /= -1 .and. numr /= -2 .and. numr /= -3) then
   write(6,*) trim(here),': numr=',numr,' /= 0, - 1, -2, -3'
   write(6,*) trim(here),': implemented only for the case numr = 0 or -1 or -2 or -3'
   write(6,*) trim(here),': i.e. analytical slater functions with order: 1s, 2s, 2p, 3s, 3p, 3d, ... for 0'
   write(6,*) trim(here),': i.e. analytical slater functions with order: 1s, 2s, 3s, ...., 2p, 3p... for -1, -2, -3'
   call die (here)
  endif
  if (ncent /= 1 ) then
   call die (here, 'ncent='+ncent+' /= 1. Implemented only for one center')
  endif
! end requirements

! allocation
  call object_alloc ('basis_ovlp', basis_ovlp, nbasis, nbasis)

  do bas_i = 1, nbasis
    do bas_j = bas_i, nbasis

      n_i   = abs(n_bas (bas_i))
      n_j   = abs(n_bas (bas_j))
      l_i   = l_bas (bas_i)
      l_j   = l_bas (bas_j)
      m_i   = m_bas (bas_i)
      m_j   = m_bas (bas_j)
      exp_i = zex (bas_i, 1)
      exp_j = zex (bas_j, 1)

      basis_ovlp (bas_i, bas_j) = slater_ovlp (n_i, l_i, m_i, exp_i, n_j, l_j, m_j, exp_j)

      basis_ovlp (bas_j, bas_i) = basis_ovlp (bas_i, bas_j)

    enddo ! bas_j
  enddo ! bas_i

!  do bas_i = 1, nbasis
!    do bas_j = 1, nbasis
!     write(6,'(a,a,i,a,i,a,f)') trim(here),': bas_i=',bas_i,' bas_j=',bas_j,' basis_ovlp=',basis_ovlp (bas_i, bas_j)
!    enddo ! bas_j
!  enddo ! bas_i

  end subroutine basis_ovlp_bld

! ==============================================================================
  subroutine basis_ovlp_eig_bld
! ------------------------------------------------------------------------------
! Description   : Eigensystem of basis overlap matrix
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_j, lin_dep_nb
  real(dp) lin_dep_thres

! header
  if (header_exe) then

   call object_create ('basis_ovlp_eigvec')
   call object_create ('basis_ovlp_eigval')

   call object_needed ('nbasis')
   call object_needed ('basis_ovlp')

   return

  endif

! begin

! allocation
  call object_alloc ('basis_ovlp_eigvec', basis_ovlp_eigvec, nbasis, nbasis)
  call object_alloc ('basis_ovlp_eigval', basis_ovlp_eigval, nbasis)

  call eigensystem (basis_ovlp, basis_ovlp_eigvec, basis_ovlp_eigval, nbasis)

  write(6,'(a)') 'Eigenvalues of overlap matrix of unnormalized basis functions:'
  do bas_i = 1, nbasis
     write(6,'(a,i3,a,es15.8)') 'eigenvalue # ',bas_i,' : ', basis_ovlp_eigval (bas_i)
  enddo ! bas_i
  write(6,'(a)') 'Eigenvectors:'
  do bas_i = 1, nbasis
     write(6,'(a,i3,a,100f12.6)') 'eigenvector # ', bas_i,' :', (basis_ovlp_eigvec (bas_i,bas_j), bas_j = 1, nbasis)
  enddo ! bas_i

! check linear dependancies
  lin_dep_thres = 1.d-12
  lin_dep_nb = 0
  do bas_i = 1, nbasis
     if (basis_ovlp_eigval (bas_i) < lin_dep_thres) then
       lin_dep_nb = lin_dep_nb + 1
     endif
  enddo ! bas_i
  if (lin_dep_nb > 0) then
   write(6,'(a,i3,a,es15.8)') 'Warning: there are ',lin_dep_nb,' eigenvalues < ',lin_dep_thres
  endif

  end subroutine basis_ovlp_eig_bld

! ==============================================================================
  subroutine basis_ovlp_12_bld
! ------------------------------------------------------------------------------
! Description   : overlap matrix to the power 1/2 for symmetric orthonormalization
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_j, bas_k

! header
  if (header_exe) then

   call object_create ('basis_ovlp_12')

   call object_needed ('nbasis')
   call object_needed ('basis_ovlp_eigvec')
   call object_needed ('basis_ovlp_eigval')

   return

  endif

! begin

! allocation
  call object_alloc ('basis_ovlp_12', basis_ovlp_12, nbasis, nbasis)

  do bas_i = 1, nbasis
   do bas_j = 1, nbasis
     basis_ovlp_12 (bas_i, bas_j) = 0.d0
     do bas_k = 1, nbasis
      basis_ovlp_12 (bas_i, bas_j) = basis_ovlp_12 (bas_i, bas_j) + basis_ovlp_eigvec (bas_i, bas_k) * dsqrt(basis_ovlp_eigval (bas_k)) * basis_ovlp_eigvec (bas_j, bas_k)
     enddo ! bas_k
   enddo ! bas_j
  enddo ! bas_i

  end subroutine basis_ovlp_12_bld

! ==============================================================================
  subroutine basis_ovlp_m12_bld
! ------------------------------------------------------------------------------
! Description   : overlap matrix to the power -1/2 for symmetric orthonormalization
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_j, bas_k

! header
  if (header_exe) then

   call object_create ('basis_ovlp_m12')

   call object_needed ('nbasis')
   call object_needed ('basis_ovlp_eigvec')
   call object_needed ('basis_ovlp_eigval')

   return

  endif

! begin

! allocation
  call object_alloc ('basis_ovlp_m12', basis_ovlp_m12, nbasis, nbasis)

  do bas_i = 1, nbasis
   do bas_j = 1, nbasis
     basis_ovlp_m12 (bas_i, bas_j) = 0.d0
     do bas_k = 1, nbasis
      basis_ovlp_m12 (bas_i, bas_j) = basis_ovlp_m12 (bas_i, bas_j) + basis_ovlp_eigvec (bas_i, bas_k) * (1.d0/dsqrt(basis_ovlp_eigval (bas_k))) * basis_ovlp_eigvec (bas_j, bas_k)
     enddo ! bas_k
   enddo ! bas_j
  enddo ! bas_i

  end subroutine basis_ovlp_m12_bld

! ========================================================================
  function factorial (i) result(result)
!-------------------------------------------------------------------------
! Description   : factorial function
!
! Created       : J. Toulouse, 17 Oct 2005
!-------------------------------------------------------------------------
  implicit none

  integer i, j
  integer result

  if (i < 0) then
    call die ('factorial with negative argument='+i)
  endif

  result = 1

  do j = 2,i
     result = result * j
  enddo

  end function factorial

end module basis_mod

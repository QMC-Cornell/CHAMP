module projector

  use bspline_oo_module
  use, intrinsic :: iso_fortran_env, only: rk => real64

  implicit none

  ! Configuration
  integer, parameter :: TITLE_MAX_LENGTH = 16
  integer, parameter :: FILENAME_MAX_LENGTH = 64
  integer, parameter :: BUF_LENGTH = 64, BUF_LENGTH2 = 16
! real(rk), parameter :: E_T = -2.175_rk
  real(rk), parameter :: MU = 1.0_rk
! real(rk), parameter :: grid_z = 0.0_rk
  real(rk), parameter :: PI = 4.0_rk*atan(1._rk)

  ! Global variables
  character(TITLE_MAX_LENGTH) :: title
  integer, save :: num_type ! number of types of nuclei
  integer, save :: tot_nuclei
  integer, save, allocatable :: num_nuclei(:) ! number of nuclei each type
  integer, save, allocatable :: z_nuclei(:) ! Z of each type
  real(rk), save, allocatable :: pos_nuclei(:, :) ! position of nuclei
  integer :: num_elec, num_spin_up, num_spin_dn
  real(rk), allocatable :: pos_elec_init(:, :) ! initial position of electrons
  real(rk), allocatable :: pos_elec_final(:, :) ! final position of electrons
! integer :: num_grid
! real(rk) :: grid_xmin, grid_xmax, grid_ymin, grid_ymax
! real(rk) :: grid_dx, grid_dy
  real(rk), allocatable :: proj(:, :)
  real(rk), allocatable :: g_eZ_up(:, :), g_eZ_dn(:, :), g_eZ(:, :)
  integer, allocatable :: recur_used(:)
  integer, allocatable :: recur_track(:)
  real(rk) :: tau, E_T
  character(BUF_LENGTH2) :: tau_str
  type(bspline_2d) :: pot_ee
  type(bspline_2d), allocatable :: pot_eZ(:)
  integer :: antisym_method=3

! call get_cmd_input()
! call read_input()
! call read_potential()
! call evaluate_proj_grid()
! call output_proj()

  contains

! subroutine interface_projector(num_type, num_nuclei, pos_nuclei, z_nuclei, num_elec, num_spin_up, pos_elec_init, pos_elec_final, tau, E_T)
  subroutine interface_projector(nctype, ncent, iwctype, cent, znuc, nelec, nup, tau_in, etrial)
!   implicit none
    implicit real*8(a-h,o-z)
    integer, intent(in) :: nctype, ncent, iwctype(ncent), nelec, nup
    real(rk), intent(in) :: cent(3,ncent), znuc(ncent), tau_in, etrial
    integer i

    num_type=nctype; tot_nuclei=ncent; num_elec=nelec; num_spin_up=nup; num_spin_dn=nelec-nup; tau=tau_in; E_T=etrial

    write(6,'(''znuc='',9f8.4)') znuc(1:ncent)

    allocate(num_nuclei(num_type))
    allocate(z_nuclei(num_type))
    allocate(pos_nuclei(tot_nuclei, 3))
    allocate(pos_elec_init(num_elec, 3))
    allocate(pos_elec_final(num_elec, 3))
    !allocate(pot_eZ(num_type))
    !allocate(g_eZ(num_elec, num_elec))
    allocate(g_eZ_up(num_spin_up, num_spin_up))
    allocate(g_eZ_dn(num_spin_dn, num_spin_dn))
    allocate(recur_used(num_elec))
    allocate(recur_track(num_elec))

    num_type=1
    num_nuclei(1:num_type)=1
    do i=1,ncent
      if(i.ge.2) then
        if(znuc(i).eq.z_nuclei(iwctype(i-1))) then
          num_nuclei(num_type)=num_nuclei(num_type)+1
        else
          num_type=num_type+1
        endif
        write(6,'(''i,num_type='',9i5)') i,num_type
      endif
      !z_prev=znuc(i)
      z_nuclei(iwctype(i))=nint(znuc(i))
      write(6,'(''i,iwctype(i),z_nuclei(iwctype(i))'',9i5)') i,iwctype(i),z_nuclei(iwctype(i))
      pos_nuclei(i,1:3)=cent(1:3,i)
    enddo

    call read_potential()

  end subroutine interface_projector
!--------------------------------------------------------------------------------------------------
  function do_projector(xold, xnew)
!   implicit none
    implicit real*8(a-h,o-z)
    real(rk), intent(in) :: xold(3,num_elec), xnew(3,num_elec)
    integer i

    !pos_elec_init(1:3, 1:num_elec)=xold(1:3, 1:num_elec)
    !pos_elec_final(1:3, 1:num_elec)=xnew(1:3, 1:num_elec)
    do i=1,num_elec
      pos_elec_init(i, 1:3)=xold(1:3, i)
      pos_elec_final(i, 1:3)=xnew(1:3, i)
    enddo

    do_projector = evaluate_proj()

  end function do_projector
!--------------------------------------------------------------------------------------------------

! subroutine get_cmd_input()

!   implicit none

!   integer :: i
!   character(len = BUF_LENGTH) :: arg

!   if (iargc() == 0) then
!     write (*, '("Format: program_name antisym_method#")')
!     stop
!   endif

!   do i = 1, iargc()
!     call getarg(i, arg)
!     write (*, '("Command #", A, ": ", A)') trim(str(i)), arg

!     if (i == 1) then
!       read (arg, *) antisym_method
!       write (*, '("Anti Symmetrization Method: ", A)'), trim(str(antisym_method))
!     endif
!   enddo

! end subroutine get_cmd_input
!--------------------------------------------------------------------------------------------------

! subroutine read_input()

!   implicit none

!   integer :: i

!   write (*, '("==== Start Reading Input ====")')

!   ! read info of nuclei
!   read (*, '(A)') title
!   write (*, '("Title: ", A)') title
!   read (*, *) num_type
!   write (*, '("Number of Types: ", A)') trim(str(num_type))
!   allocate(num_nuclei(num_type))
!   allocate(z_nuclei(num_type))

!   read (*, *) (num_nuclei(i), i = 1, num_type)
!   read (*, *) (z_nuclei(i), i = 1, num_type)
!   tot_nuclei = 0
!   do i = 1, num_type
!     tot_nuclei = tot_nuclei + num_nuclei(i)
!     write (*, '("Type, Count, Z: ", 3I5)') i, num_nuclei(i), z_nuclei(i)
!   enddo

!   allocate(pos_nuclei(tot_nuclei, 3))
!   write (*, '("Nuclei Positions:")')
!   do i = 1, tot_nuclei
!     read (*, *) pos_nuclei(i, 1:3)
!     write (*, '("X, Y, Z: ", 3F10.6)') pos_nuclei(i, 1:3)
!   enddo

!   ! read info of electrons
!   read (*, *) num_elec, num_spin_up
!   write (*, '("Numbber of Electrons: ", A)') trim(str(num_elec))
!   num_spin_dn=num_elec-num_spin_up
!   write (*, '("Number of up, dn spin electrons: ", A)') trim(str(num_spin_up)), trim(str(num_spin_dn))
!   allocate(pos_elec_init(num_elec, 3))
!   allocate(pos_elec_final(num_elec, 3))
!   write (*, '("Initial Position of Electrons:")')
!   do i = 1, num_elec
!     read (*, *) pos_elec_init(i, 1:3)
!     write (*, '("X, Y, Z: ", 3F10.6)') pos_elec_init(i, 1:3)
!   enddo
!   write (*, '("Final Position of Electrons:")')
!   do i = 1, num_elec - 1
!     read (*, *) pos_elec_final(i, 1: 3)
!     write (*, '("X, Y, Z: ", 3F10.6)') pos_elec_final(i, 1:3)
!   enddo

!   ! read output grid
!   read (*, *) num_grid
!   read (*, *) grid_xmin, grid_xmax
!   read (*, *) grid_ymin, grid_ymax
!   write (*, '("Number of Grid Points: ", A)'), trim(str(num_grid))
!   write (*, '("x_min, x_max: ", 2F10.6)') grid_xmin, grid_xmax
!   write (*, '("y_min, y_max: ", 2F10.6)') grid_ymin, grid_ymax

!   ! read tau
!   read (*, '(A)') tau_str
!   read (tau_str, *) tau
!   write (*, '("Tau: ", F10.6)') tau

!   write (*, '("==== Finished Reading Input ====")')

! end subroutine read_input
!--------------------------------------------------------------------------------------------------

  subroutine read_potential()

    implicit none

    character(FILENAME_MAX_LENGTH) :: file_e, file_Z
    integer :: i

    write (*, '("==== Start Loading Potential ====")')
    file_e = get_filename(-1)
    call read_potential_sub(file_e, pot_ee)

    allocate(pot_eZ(num_type))
    do i = 1, num_type
      file_Z = get_filename(z_nuclei(i))
      call read_potential_sub(file_Z, pot_eZ(i))
    enddo

    write (*, '("==== Finished Loading Potential ====")')

  end subroutine read_potential
!--------------------------------------------------------------------------------------------------

  subroutine read_potential_sub(filename, bspline_obj)
    implicit none

    integer :: grid_in
    real(rk) :: r_0, r_n
    integer, parameter :: LINE_SKIP = 1
    real(rk), allocatable :: x_in(:), y_in(:), u_in(:, :)
    real(rk) :: tx, ty, tu ! temp variables
    character(FILENAME_MAX_LENGTH), intent(in) :: filename
    integer :: i, j
    type(bspline_2d) :: bspline_obj
    integer :: kx, ky
    integer :: iflag

    write (*, '("Loading Potential From: ", A)') trim(filename)

    ! read potential from filename
    open (1, file = trim(filename),status='old')
    rewind(1) ! move to the beginning of file
    read (1, *)
    read (1, *) grid_in, r_0, r_n

    allocate(x_in(grid_in))
    allocate(y_in(grid_in))
    allocate(u_in(grid_in, grid_in))

    do i = 1, LINE_SKIP
      read (1, *)
    enddo

    do i = 1, grid_in
      do j = 1, grid_in
        read (1, *) tx, ty, tu
        if (isnan(tu) .or. isnan(tx) .or. isnan(ty)) then
          write (*, '("Invalid Potential (tx, ty, tu): ", 3F15.6)') tx, ty, tu
          stop
        endif
        u_in(i, j) = tu
        if (i == 1) then
          y_in(j) = ty
        endif
        if (j == 1) then
          x_in(i) = tx
        endif
      enddo
    enddo

    close(1)

    ! spline interpolation
    kx = 3 ! use cubic spline
    ky = 3
    call bspline_obj%initialize(x_in, y_in, u_in, kx, ky, iflag)

    if(iflag == 1) then
      write (*, '("Interpolation Succeeded!")')
    else
      write (*, '("Interpolation Failed!")')
      call exit(0)
    endif
  end subroutine read_potential_sub
!--------------------------------------------------------------------------------------------------

  function str(num_in)
    integer, intent(in) :: num_in
    character(BUF_LENGTH) :: str
    write (str, *) num_in
    str = adjustl(str)
  end function str
!--------------------------------------------------------------------------------------------------

  function get_filename(Z_in)

    implicit none

    integer, intent(in) :: Z_in
    character(FILENAME_MAX_LENGTH) :: get_filename
    character(FILENAME_MAX_LENGTH) :: buf

    ! obtain folder name
    !write (buf1, '("Z=", A, "_tau=", A, "_grid=400")'), trim(str(Z_nuclei)), trim(tau_str)

    if(nint(log10(tau)).eq.0) tau_str='1'
    if(nint(log10(tau)).eq.1) tau_str='10'
    if(nint(log10(tau)).eq.2) tau_str='100'
    if(nint(log10(tau)).eq.-1) tau_str='.1'
    if(nint(log10(tau)).eq.-2) tau_str='.01'

    ! obtain file name
    if(Z_in > 0) then
      if (tau >= 100) then
        write (buf, '("u_", A, "_tau=", A, "_ngrid=400")') trim(str(Z_in)), trim(tau_str)
      else
        write (buf, '("u_", A, "_tau=", A, "_ngrid=200")') trim(str(Z_in)), trim(tau_str)
      endif
    elseif(Z_in == -1) then
      write (buf, '("u_e_tau=", A, "_ngrid=200")') trim(tau_str)
    endif

    ! combine folder and file name to filename with path
    write (get_filename, '(A)') trim(buf)

  end function get_filename
!--------------------------------------------------------------------------------------------------

! subroutine evaluate_proj_grid()
!   ! Calculate position corresponding to each grid point
!   ! Change the position of last electron to position calculated
!   ! Call evaluate_proj to evaluate proj
!   ! Save the result to proj

!   implicit none

!   integer :: i, j
!   real(rk) :: tx, ty
!   integer :: elec_id

!   write (*, '("Evaluate Projector...")')

!   elec_id = num_elec

!   if(.not. allocated(proj)) allocate(proj(num_grid, num_grid))
!  ! if(.not. allocated(g_eZ)) allocate(g_eZ(num_elec, num_elec))

!   if(.not. allocated(recur_used)) allocate(recur_used(num_elec))
!   if(.not. allocated(recur_track)) allocate(recur_track(num_elec))

!   grid_dx = (grid_xmax - grid_xmin) / num_grid
!   grid_dy = (grid_ymax - grid_ymin) / num_grid

!   do i = 1, num_grid
!     tx = grid_xmin + grid_dx * (i - 1)
!     do j = 1, num_grid
!       ty = grid_ymin + grid_dy * (j - 1)
!       pos_elec_final(num_elec, 1: 3) = (/tx, ty, grid_z/)
!       proj(i, j) = evaluate_proj()
!     enddo
!   enddo

! end subroutine evaluate_proj_grid
!--------------------------------------------------------------------------------------------------

  function evaluate_proj()
    ! Evaluate proj for the current configuration of particles

    implicit none

    real(rk) :: evaluate_proj
    real(rk) :: u_ee
    real(rk) :: g_eZ_det_up, g_eZ_det_dn, t

    integer :: i, j, shift
!   real(rk), save :: g_sav
    !!if(.not. allocated(g_eZ)) allocate(g_eZ(num_elec, num_elec))

! Set up separate g_eZ_up and g_eZ_dn matrices for all initial and final electron positions within a spin type.
! Note: In the case of method3, g_eZ is a misnomer because it includes parallel-spin g_ee.
    if(antisym_method.ge.0.and.antisym_method.le.3) then

      if(num_spin_up.ne.0) then
        !if(.not. allocated(g_eZ_up)) allocate(g_eZ_up(num_spin_up, num_spin_up))
        do i = 1, num_spin_up
          do j = 1, num_spin_up
            t = get_u_eZ(i, j)
            if(antisym_method.eq.3) then
              t = t + get_u_ee3(i, j, 1, num_spin_up)
            endif
            g_eZ_up(i, j) = exp(-t)
          enddo
        enddo
      endif

      if(num_spin_up.ne.num_elec) then
        !if(.not. allocated(g_eZ_dn)) allocate(g_eZ_dn(num_spin_dn, num_spin_dn))
        do i = 1, num_spin_dn
          do j = 1, num_spin_dn
            t = get_u_eZ(num_spin_up + i, num_spin_up + j)
            if(antisym_method.eq.3) then
              t = t + get_u_ee3(num_spin_up + i, num_spin_up + j, num_spin_up + 1, num_elec)
            endif
            g_eZ_dn(i, j) = exp(-t)
          enddo
        enddo
      endif

    endif

!   do i = 1, num_elec
!     do j = 1, num_elec
!       !g_eZ(i, j) = get_u_eZ(i, j)
!       t = get_u_eZ(i, j)
!       if(antisym_method.eq.3) then
!         t = t + get_u_ee3(i, j)
!       endif
!       g_eZ(i, j) = exp(-t)
!     enddo
!   enddo

    g_eZ_det_up = 1.0_rk
    g_eZ_det_dn = 1.0_rk
    ! matinv changes the matrix but since it is not called for method0, we do not need to save a copy.
    if(antisym_method.ge.1 .and. antisym_method.le.3) then
      if(num_spin_up.ne.0) then
        call matinv(g_eZ_up, num_spin_up, g_eZ_det_up)
        !write(6,*) 'g_eZ_det_up = ', g_eZ_det_up
      endif
      if(num_spin_up.ne.num_elec) then
        call matinv(g_eZ_dn, num_spin_dn, g_eZ_det_dn)
        !write(6,*) 'g_eZ_det_dn = ', g_eZ_det_dn
      endif
    endif

    ! Calculate based on which antisym method selected
    select case(antisym_method)
    case (0)
!     do i = 1, num_elec
!       recur_used(i) = 0
!     enddo
!     evaluate_proj = evaluate_proj_recur(1, recur_used, recur_track)
      u_ee = get_u_ee3_diff_spin()
      if (num_spin_up.ge.1) then
        recur_used(1:num_spin_up) = 0
        shift=0
        g_eZ_det_up = evaluate_proj_recur(1, recur_used, recur_track, shift, num_spin_up)
      else
        g_eZ_det_up = 1
      endif
      if (num_spin_dn.ge.1) then
        recur_used(1:num_spin_dn) = 0
        shift=num_spin_up
        g_eZ_det_dn = evaluate_proj_recur(1, recur_used, recur_track, shift, num_spin_dn)
      else
        g_eZ_det_dn = 1
      endif
!     write(6,'(''u_ee,g_eZ_det_up,g_eZ_det_dn='',9es12.4)') u_ee,g_eZ_det_up,g_eZ_det_dn
    case (1)
      u_ee = get_u_ee1()
    case (2)
      u_ee = get_u_ee2()
    case (3)
      u_ee = get_u_ee3_diff_spin()
    case default
      stop 'In evaluate_proj, only cases 0-3 are implemented'
    end select
    evaluate_proj = g_eZ_det_up * g_eZ_det_dn * exp(E_T*tau - u_ee)
    !write(6,'(''evaluate_proj, g_eZ_det_up, g_eZ_det_dn, exp(E_T*tau - u_ee), u_ee, E_T'',9es12.4)') evaluate_proj, g_eZ_det_up, g_eZ_det_dn, exp(E_T*tau - u_ee), u_ee, E_T

    ! if (u_ee_method == 3 .and. num_elec == 2) then
    !  evaluate_proj = (g_eZ(1, 1) * g_eZ(2, 2) * exp(-get_u_ee(3, 1)) - g_eZ(1, 2) * g_eZ(2, 1) * exp(-get_u_ee(3, 2))) * exp(E_T * tau)
    !  write (*, *) "get_u_ee3", -get_u_ee(3, 1), -get_u_ee(3, 2)
    !else
    !  u_ee = get_u_ee(u_ee_method, 0)
    !  evaluate_proj = g_eZ_det * exp(E_T * tau - u_ee)
    !endif

!   ifi(abs(g_eZ_det).gt.g_sav) then
!     write(6,'(''g_eZ_det, exp(E_T * tau - u_ee), 1/ (2*PI*tau)**(3._rk*num_elec/2._rk)'',9es12.4)') g_eZ_det, exp(E_T * tau - u_ee), 1/ (2*PI*tau)**(3._rk*num_elec/2._rk)
!     g_sav=abs(g_eZ_det)
!   endif

    evaluate_proj = evaluate_proj * (MU / (2 * PI * tau))**(1.5_rk * num_elec)

  end function evaluate_proj
!--------------------------------------------------------------------------------------------------

  recursive function evaluate_proj_recur(level, recur_used, recur_track, shift, n) result (ret)

    implicit none

    real(rk) :: ret, ret_prev
    integer :: level, shift, n
    integer :: i
    integer :: sn
    integer :: recur_used(n), recur_track(n)

    if (level > n) then
      ret = 1.0_rk
      do i = 1, n
        if(shift.eq.0) then
          ret = ret * g_eZ_up(recur_track(i), i)
        else
          ret = ret * g_eZ_dn(recur_track(i), i)
        endif
      enddo
      ret = ret * exp(-get_u_ee0(recur_track, shift, n))
      return
    endif

    sn = 1
    ret = 0.0_rk
    do i = 1, n
      if (recur_used(i) == 1) then
        cycle
      endif
      recur_used(i) = 1
      recur_track(level) = i
      ret_prev = evaluate_proj_recur(level + 1, recur_used, recur_track, shift, n)
      ret = ret + ret_prev * sn
      recur_used(i) = 0
      sn = -sn
    enddo

  end function evaluate_proj_recur
!--------------------------------------------------------------------------------------------------

  function get_u_eZ(elec_to, elec_from)

    implicit none

    integer, intent(in) :: elec_to, elec_from
    real(rk) :: get_u_eZ
    integer :: i, type_id, type_cnt
    real(rk) :: q, s, t
    real(rk) :: r_to_a(3), r_from_a(3)
    integer :: iflag

    get_u_eZ = 0.0_rk

    type_id = 1
    type_cnt = 0

    s = norm2(pos_elec_final(elec_to,1:3) - pos_elec_init(elec_from,1:3))

    do i = 1, tot_nuclei
      r_to_a = pos_elec_final(elec_to, 1:3) - pos_nuclei(i, 1:3)
      r_from_a = pos_elec_init(elec_from, 1:3) - pos_nuclei(i, 1:3)
      q = (norm2(r_to_a) + norm2(r_from_a)) / 2
      call pot_eZ(type_id)%evaluate(q + s / 2, q - s / 2, 0, 0, t, iflag)

      type_cnt = type_cnt + 1

      if (type_cnt == num_nuclei(type_id)) then
        type_id = type_id + 1
        type_cnt = 0
      endif

      get_u_eZ = get_u_eZ + t
    enddo

    get_u_eZ = get_u_eZ + MU * s**2 / (2 * tau)

    !get_u_eZ = exp(-get_u_eZ)

  end function get_u_eZ
!--------------------------------------------------------------------------------------------------

  function get_u_ee0(recur_track, shift, n)
    ! Evaluate u_ee according to:
    ! Umrigar 2015. "Observation on variational and projector Monte Carlo Method"
    ! With full antisymmetrization

    implicit none

    real(rk) :: get_u_ee0
    integer, intent(in) :: recur_track(:), shift, n
    integer :: i, track_i, j, track_j
    real(rk) :: q, s
    real(rk) :: t
    real(rk) :: r_ij(3), r_ij_prime(3)
    integer :: iflag

    get_u_ee0 = 0.0_rk

    do i = 1, n
      track_i = recur_track(i)
      do j = i + 1, n
        track_j = recur_track(j)
        r_ij = pos_elec_init(shift+j, 1:3) - pos_elec_init(shift+i, 1:3)
        r_ij_prime = pos_elec_final(shift+track_j, 1:3) - pos_elec_final(shift+track_i, 1:3)
        q = (norm2(r_ij) + norm2(r_ij_prime)) / 2
        s = norm2(r_ij - r_ij_prime)
        call pot_ee%evaluate(q + s / 2, q - s / 2, 0, 0, t, iflag)
        get_u_ee0 = get_u_ee0 + t
      enddo
    enddo

    ! Only deals with 2 electrons for now
    !if (num_elec .ne. 2) then
    !  return
    !endif

    !i = 1
    !j = 2
    !r_ij = pos_elec_init(j, 1:3) - pos_elec_init(i, 1:3)
    !r_ij_prime = pos_elec_final(j, 1:3) - pos_elec_final(i, 1:3)
    !q = (norm2(r_ij) + norm2(r_ij_prime)) / 2

    !if (arg .eq. 1) then
    !  s = norm2(r_ij - r_ij_prime)
    !else
    !  s = norm2(r_ij + r_ij_prime)
    !endif

    !call pot_ee%evaluate(q + s / 2, q - s / 2, 0, 0, t, iflag)
    !get_u_ee0 = t

  end function get_u_ee0
!--------------------------------------------------------------------------------------------------
  function get_u_ee1()
    ! Evaluate u_ee according to:
    ! Umrigar 2015. "Observation on variational and projector Monte Carlo Method", Eq. 17.

    implicit none

    real(rk) :: get_u_ee1
    integer :: i, j
    real(rk) :: q_prime, q, s
    real(rk) :: t_prime, t
    real(rk) :: r_ij(3), r_ij_prime(3)
    integer :: iflag

    get_u_ee1 = 0.0_rk
    s = 0.0_rk

    do i = 1, num_elec
      do j = i + 1, num_elec
        r_ij = pos_elec_init(j, 1:3) - pos_elec_init(i, 1:3)
        r_ij_prime = pos_elec_final(j, 1:3) - pos_elec_final(i, 1:3)
        q_prime = norm2(r_ij_prime)
        q = norm2(r_ij)
        call pot_ee%evaluate(q, q, 0, 0, t, iflag)
        call pot_ee%evaluate(q_prime, q_prime, 0, 0, t_prime, iflag)
        get_u_ee1 = get_u_ee1 + (t + t_prime) / 2
      enddo
    enddo

  end function get_u_ee1
!--------------------------------------------------------------------------------------------------

  function get_u_ee2()
    ! Evaluate u_ee according to:
    ! Umrigar 2015. "Observation on variational and projector Monte Carlo Method", Eq. 18.

    implicit none

    real(rk) :: get_u_ee2
    integer :: i, j
    real(rk) :: q, s
    real(rk) :: t
    real(rk) :: r_ij(3), r_ij_prime(3)
    integer :: iflag

    get_u_ee2 = 0.0_rk

    do i = 1, num_elec
      do j = i + 1, num_elec
        r_ij = pos_elec_init(j, 1:3) - pos_elec_init(i, 1:3)
        r_ij_prime = pos_elec_final(j, 1:3) - pos_elec_final(i, 1:3)
        q = (norm2(r_ij) + norm2(r_ij_prime)) / 2
        s = norm2(r_ij - r_ij_prime)
        call pot_ee%evaluate(q + s / 2, q - s / 2, 0, 0, t, iflag)
        get_u_ee2 = get_u_ee2 + t
      enddo
    enddo

  end function get_u_ee2
!--------------------------------------------------------------------------------------------------
  function get_u_ee3_diff_spin()
    ! Evaluate u_ee for opposite spin pairs according to method 1
    implicit none

    real(rk) :: get_u_ee3_diff_spin
    integer :: i, j
    real(rk) :: q, s
    real(rk) :: t
    real(rk) :: r_ij(3), r_ij_prime(3)
    integer :: iflag

    get_u_ee3_diff_spin = 0.0_rk
    s = 0.0_rk

    do i = 1, num_spin_up
      do j = num_spin_up + 1, num_elec
! This is the same as for method 1, but what we want is method 2
!       r_ij = pos_elec_init(j, 1:3) - pos_elec_init(i, 1:3)
!       r_ij_prime = pos_elec_final(j, 1:3) - pos_elec_final(i, 1:3)
!       q_prime = norm2(r_ij_prime)
!       q = norm2(r_ij)
!       call pot_ee%evaluate(q, q, 0, 0, t, iflag)
!       call pot_ee%evaluate(q_prime, q_prime, 0, 0, t_prime, iflag)
!       get_u_ee3_diff_spin_diff_spin = get_u_ee3_diff_spin_diff_spin + (t + t_prime) / 2
! Same as in method 2
        r_ij = pos_elec_init(j, 1:3) - pos_elec_init(i, 1:3)
        r_ij_prime = pos_elec_final(j, 1:3) - pos_elec_final(i, 1:3)
        q = (norm2(r_ij) + norm2(r_ij_prime)) / 2
        s = norm2(r_ij - r_ij_prime)
        call pot_ee%evaluate(q + s / 2, q - s / 2, 0, 0, t, iflag)
        get_u_ee3_diff_spin = get_u_ee3_diff_spin + t
      enddo
    enddo

  end function get_u_ee3_diff_spin
!
!--------------------------------------------------------------------------------------------------

  function get_u_ee3(elec_to, elec_from, start_index, end_index)
    ! Evaluate u_ee according to:
    ! Umrigar 2015. "Going beyond the fixed-node approximation in DMC"

    implicit none

    real(rk) :: get_u_ee3
    integer, intent(in) :: elec_to, elec_from, start_index, end_index
    integer :: k
    real(rk) :: q, s
    real(rk) :: t
    real(rk) :: r_ij(3), r_ij_prime(3)
    integer :: iflag

    get_u_ee3 = 0._rk
   !Loop not evaluated if start_index == end_index
    do k = start_index + 1, end_index
        r_ij = pos_elec_init(elec_from,1:3) - pos_elec_init(map(elec_from,k, start_index),1:3)
        r_ij_prime = pos_elec_final(elec_to, 1:3) - pos_elec_final(map(elec_to,k, start_index), 1:3)
        q = (norm2(r_ij) + norm2(r_ij_prime)) / 2
        s = norm2(r_ij - r_ij_prime)
        call pot_ee%evaluate(q + s / 2, q - s / 2, 0, 0, t, iflag)
        get_u_ee3 = get_u_ee3 + t
    enddo

    get_u_ee3 = .5_rk*get_u_ee3

  end function get_u_ee3
!-------------------------------------------------------------------------------------------------

  integer function map(i,k, start_index)
    ! Evaluate M according to:
    ! Umrigar 2015. "Going beyond the fixed-node approximation in DMC"

    implicit none

    integer, intent(in) :: i,k, start_index

    if(k.eq.i) then
      map=start_index
    else
      map=k
    endif

  end function map
!--------------------------------------------------------------------------------------------------

! subroutine output_proj()
!   ! Output calculated proj to {title}.proj

!   implicit none

!   character(FILENAME_MAX_LENGTH) :: filename
!   integer :: i, j
!   real(rk) :: tx, ty

!   write(filename, '(A, "_", A, "_M", A, ".proj")') trim(title), trim(tau_str), trim(str(antisym_method))

!   open(2, file = trim(filename))

!   write (2, '(3A20)') "x", "y", "proj"

!   do i = 1, num_grid
!     tx = grid_dx * (i - 1) + grid_xmin
!     do j = 1, num_grid
!       ty = grid_dy * (j - 1) + grid_ymin
!       write (2, '(3ES20.10)') tx, ty, proj(i, j)
!     enddo
!     write (2, *)
!   enddo

!   close(2)

!   write (*, '("Results Saved To: ", A)') trim(filename)

! end subroutine output_proj

end module projector

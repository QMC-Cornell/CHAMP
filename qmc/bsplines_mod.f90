! things to do:
!       move the multiple file reads and file writes into two single subroutines
!               read_file, write_file --need an interface that can handle
!               either real or complex arrays and varying numbers of arrays(1-5)
!       eliminate goto statements and numbered lines(by creating above
!               routines for file i/o)
!       eliminate unnecessary case l_need_phase = .true., nkvec=1
!       take advantage of antiperiodic boundary conditions in einspline
!       use FFT to get plane wave values on real space grid
module bsplines_mod

  use basic_tools_mod
  use mpi_mod
  use orbital_grid_mod, only: igrad_lap,ngrid_orbx,ngrid_orby,ngrid_orbz

! Declaration of global variables and default values
  logical                            :: l_need_phase = .false.
  integer                            :: norb_primitive
  integer                            :: nspline ! 1=orbitals, 2=orbitals' Laplacian, 3-5=orbitals' gradient
  integer(dp)                        :: bspline_coefficient_pointer(5)
  real(dp)                           :: rlatt_inv_transpose(3,3)
  real(dp)                           :: rlatt_inv_transpose_2(3,3)
  real(dp),parameter                 :: zero_tolerance = 1.d-5

  contains

!===========================================================================
  subroutine setup_bsplines_coefficients
!---------------------------------------------------------------------------
! Description : Setup variables and plane-wave coefficient array to pass to
!               einspline library for calculating B-spline coefficients
!
! Created     : W. Parker, 21 Oct 2008
!---------------------------------------------------------------------------
  implicit none
! need ndim, norb, nband
  include 'commons.h'

! local
  integer                             :: idum
  integer                             :: ix
  integer                             :: iy
  integer                             :: iz
  integer                             :: iorb
  integer                             :: ikvec
  integer                             :: igrid_bsplines
  integer                             :: ispline
  integer                             :: ngrid
  integer                             :: num_orb_exist
  integer                             :: ngrid_orbxf ! variables
  integer                             :: ngrid_orbyf ! read in
  integer                             :: ngrid_orbzf ! from file
  integer                             :: igrad_lapf  !
  integer, allocatable                :: ngrid_orb(:)
  integer, allocatable                :: ibc(:,:)
  real(dp)                            :: k_dot_r
  real(dp), allocatable               :: r(:)
  real(dp), allocatable               :: grid_orb(:)
  real(dp), allocatable               :: rnaught(:)
  real(dp), allocatable               :: rone(:)
  real(dp), allocatable               :: rbc(:,:)
  real(dp), allocatable               ::   orb_pw(:)
  real(dp), allocatable               ::  dorb_pw(:,:)
  real(dp), allocatable               :: ddorb_pw(:)
  real(dp), allocatable               ::   orb_pw_terms(:,:,:)
  real(dp), allocatable               ::  dorb_pw_terms(:,:,:,:)
  real(dp), allocatable               :: ddorb_pw_terms(:,:,:)
  real(dp), allocatable               :: plane_wave_orbitals_at_grid_points(:,:)
  real(dp), allocatable               :: plane_wave_orbitals_laplacian_at_grid_points(:,:)
  real(dp), allocatable               :: plane_wave_orbitals_gradient1_at_grid_points(:,:)
  real(dp), allocatable               :: plane_wave_orbitals_gradient2_at_grid_points(:,:)
  real(dp), allocatable               :: plane_wave_orbitals_gradient3_at_grid_points(:,:)
  complex(dp), allocatable            :: rbc_complex(:,:)
  complex(dp), allocatable            :: complex_plane_wave_orbitals_at_grid_points(:,:,:)
  complex(dp), allocatable            :: complex_plane_wave_orbitals_laplacian_at_grid_points(:,:,:)
  complex(dp), allocatable            :: complex_plane_wave_orbitals_gradient1_at_grid_points(:,:,:)
  complex(dp), allocatable            :: complex_plane_wave_orbitals_gradient2_at_grid_points(:,:,:)
  complex(dp), allocatable            :: complex_plane_wave_orbitals_gradient3_at_grid_points(:,:,:)
  character(len=max_string_len_rout), save :: lhere = 'setup_bsplines_coefficients'

  write(6,'(''ndim,norb'',20i5)') ndim,norb,ngrid_orbx,ngrid_orby,ngrid_orbz,nkvec
! allocation
  !number of grid points
  call alloc ('ngrid_orb', ngrid_orb, ndim)
  !position in real space (used here to evaluate plane waves at grid points
  call alloc ('r', r, ndim)
  !orbital value of plane waves at particular r
  call alloc (  'orb_pw',   orb_pw,       norb)
  !values of orbitals at edge of B-splines(use crystal coordinates rnaught=0 rone=1)
  call alloc ('rnaught', rnaught, ndim)
  call alloc ('rone', rone, ndim)
  !boundary conditions for B-splines
  call alloc ('ibc', ibc, ndim, 2)
  call alloc ('rbc', rbc, ndim, 2)

! begin
  write(6,*)'Setting up interpolating B-spline coefficients'
  idum=1

  if(ndim == 3) then
    !set number of grid points in each lattice vector direction based on input
    ngrid_orb(1)=ngrid_orbx
    ngrid_orb(2)=ngrid_orby
    ngrid_orb(3)=ngrid_orbz
    !total number of grid points
    ngrid=product(ngrid_orb)
    call alloc ('grid_orb', grid_orb, ndim)
  else
    call die (lhere, 'ndim /= 3 not supported in B-splines yet')
  endif

  select case (igrad_lap)
  ! Calculate gradient and Laplacian from orbital interpolation
  case (0)

    ! Cases (One k-point=Gamma; One k-point=G/2; Multiple k-points)
    if(nkvec == 1) then
      norb_primitive=nband(1)

      !real and imaginary terms of orbital value of plane waves at particular r
      call alloc (  'orb_pw_terms', orb_pw_terms, nkvec, norb_primitive, 2)

      nspline=norb
      !if Gamma,
      ! then no phase is necessary
      ! and we only need one set of B-spline coefficients for each orbital
      if (maxval(rkvec_shift) < zero_tolerance ) then
        if(ndim == 3) then
          allocate (plane_wave_orbitals_at_grid_points(0:ngrid-1,norb_primitive))
          ! Set up an evenly spaced grid in each direction from 0 to (ngrid-1)/ngrid

          !If oritals_num_bspline exists, read pw orbitals on grid from there and skip calculating them
          num_orb_exist=0
          open(4,file='orbitals_num_bspline',form='unformatted',status='old',err=10)
          num_orb_exist=1
          write(6,'(''Bspline interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
          read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,igrad_lapf
          if(ngrid_orbxf /= ngrid_orbx .or. ngrid_orbyf /= ngrid_orby .or.  ngrid_orbzf /= ngrid_orbz) then
             write(6,'(''orbital grids do not match, program, orbitals_num:'',3i4,x,3i4)') &
                  ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
             call die (lhere, 'orbital grids do not match; change ngrid_orbx, ngrid_orby and/or ngrid_orbz')
          endif
          if(igrad_lapf /= igrad_lap) then
             write(6, '(''choice of interpolations do not match, program, orbitals_num:'',i4,x,i4)') &
                  igrad_lap,igrad_lapf
             call die (lhere, 'choice of interpolations do not match; change igrad_lap')
          endif
          read(4) ((plane_wave_orbitals_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
          write(6,'(''Done reading orbitals_num_bspline'')')
          close(4)
          goto 20

  10      do ix=1,ngrid_orb(1)
            do iy=1,ngrid_orb(2)
              do iz=1,ngrid_orb(3)
                grid_orb(1)=(ix-one)/(ngrid_orb(1))
                grid_orb(2)=(iy-one)/(ngrid_orb(2))
                grid_orb(3)=(iz-one)/(ngrid_orb(3))

                ! Convert the current control point into Cartesian coordinates
                ! because the orbitals evaluation routine needs a point in Cartesian
                ! coordinates -- use primitive cell!!
                r(:)=matmul(rlatt,grid_orb(:))

                ! Evaluate grid point on primitive cell
                call orbitals_pw_primitive(r,orb_pw_terms,dorb_pw_terms,ddorb_pw_terms)

                ! Index the array in the way einspline expects
                igrid_bsplines=(((ix-1)*ngrid_orb(2)+(iy-1))*ngrid_orb(3)+(iz-1))

                do iorb=1,norb_primitive
                  if(ireal_imag(iorb) == 1) then !use real part of plane waves
                    orb_pw(iorb)= orb_pw_terms(1,iorb,1)
                  else if(ireal_imag(iorb) == 2) then !use imaginary part of plane waves
                    orb_pw(iorb)= orb_pw_terms(1,iorb,2)
                  endif
                enddo

                ! Store the plane-wave value in the array
                plane_wave_orbitals_at_grid_points(igrid_bsplines,:)=orb_pw

              enddo !iz
            enddo !iy
          enddo !ix

   20     call release ('orb_pw', orb_pw)
          call release ('orb_pw_terms', orb_pw_terms)

        if(index(mode,'mpi') == 0 .or. idtask == 0) then
          ! If orbitals_num_bspline does not exist and inum_orb is positive, write them to that file.
          if(num_orb_exist==0 .and. inum_orb.eq.8) then
            open(4,file='orbitals_num_bspline',form='unformatted',status='new',err=30)
            write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap
            write(4) ((plane_wave_orbitals_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
            write(6,'(''Done writing orbitals_num_bspline'')')
            close(4)
          endif
        endif ! mpi or idtask == 0

          !grid point coordinates run from 0 to 1
   30     rnaught(:) = 0.d0
          rone(:) = 1.d0
          !boundary conditions (periodic: ibc=0, antiperiodic: ibc=5)
          ibc = 0
          rbc = 0.d0
          call fcreate_multi_ubspline_3d_d(rnaught(1),rone(1),ngrid_orb(1), &
                                            rnaught(2),rone(2),ngrid_orb(2), &
                                            rnaught(3),rone(3),ngrid_orb(3), &
                                            ibc(1,1),rbc(1,1),ibc(1,2),rbc(1,2), &
                                            ibc(2,1),rbc(2,1),ibc(2,2),rbc(2,2), &
                                            ibc(3,1),rbc(3,1),ibc(3,2),rbc(3,2), &
                                            nspline,bspline_coefficient_pointer(1))
          do iorb=1,nspline
             ispline=iorb-1
             call fset_multi_ubspline_3d_d(bspline_coefficient_pointer(1), ispline, &
                                           plane_wave_orbitals_at_grid_points(:,iorb))
          enddo
          write (6,*) 'B-splines coefficients successfully set up'
        endif ! ndim == 3

        call release ('rbc', rbc)
        call release ('plane_wave_orbitals_at_grid_points', plane_wave_orbitals_at_grid_points)

        !if non-Gamma but single k-point,
        ! then phase is necessary
        ! but we only need one set of B-spline coefficients for each orbital
      else !maxval(rkvec_shift) > zero_tolerance

        l_need_phase = .true.
        write(6,*)'Non-gamma k-point in use, will use phase'
        write(6,*)'Will spline real and imaginary terms separately'
        call alloc ('rbc_complex', rbc_complex, ndim, 2)

        if(ndim == 3) then
        allocate (complex_plane_wave_orbitals_at_grid_points(0:ngrid-1,nkvec,norb_primitive))

        !If orbitals_num_bspline exists, read orbitals on grid from there and skip calculating them
        num_orb_exist=0
        open(4,file='orbitals_num_bspline',form='unformatted',status='old',err=50)
        num_orb_exist=1
        write(6,'(''Bspline interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
        read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,igrad_lapf
        if(ngrid_orbxf /= ngrid_orbx .or. ngrid_orbyf /= ngrid_orby .or.  ngrid_orbzf /= ngrid_orbz) then
          write(6,'(''orbital grids do not match, program, orbitals_num:'',3i4,x,3i4)') &
   &      ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
          call die (lhere, 'orbital grids do not match; change ngrid_orbx, ngrid_orby and/or ngrid_orbz')
        endif
        read(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
          write(6,'(''Done reading orbitals_num_bspline'')')
          close(4)
!       endif
        goto 60

        ! Set up an evenly spaced grid for the control points in each direction
        !        with values from 0 to 1
  50    do ix=1,ngrid_orb(1)
           do iy=1,ngrid_orb(2)
              do iz=1,ngrid_orb(3)
                  grid_orb(1)=(ix-one)/(ngrid_orb(1))
                  grid_orb(2)=(iy-one)/(ngrid_orb(2))
                  grid_orb(3)=(iz-one)/(ngrid_orb(3))

                 ! Convert the current control point into Cartesian coordinates
                 ! because the orbitals evaluation routine needs a point in Cartesian
                 ! coordinates -- use primitive cell!!
                 r(:)=matmul(rlatt,grid_orb(:))

                 ! Evaluate grid point on primitive cell
                 call orbitals_pw_primitive(r,orb_pw_terms,dorb_pw_terms,ddorb_pw_terms)

                 ! Index the array in the way einspline expects
                 igrid_bsplines=(((ix-1)*ngrid_orb(2)+(iy-1))*ngrid_orb(3)+(iz-1))

                 ! Store the plane-wave value in the array
                 complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,:,:)=&
                 &        cmplx(orb_pw_terms(:,:,1),orb_pw_terms(:,:,2))

              enddo !iz
           enddo !iy
        enddo !ix

  60    call release ('orb_pw', orb_pw)
        call release (  'orb_pw_terms',  orb_pw_terms)

        if(index(mode,'mpi') == 0 .or. idtask == 0) then
        ! If orbitals_num_bspline does not exist and inum_orb is positive, write them to that file.
        if(num_orb_exist==0 .and. inum_orb.eq.8) then
          open(4,file='orbitals_num_bspline',form='unformatted',status='new',err=70)
          write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap
          write(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
          write(6,'(''Done writing orbitals_num_bspline'')')
          close(4)
        endif
       endif ! mpi or idtask == 0

    !grid point coordinates run from 0 to 1
  70  rnaught(:) = 0.d0
    rone(:) = 1.d0
    !boundary conditions (periodic: ibc=0, antiperiodic: ibc=5)
    ibc = 0
    rbc_complex = (0.,0.)
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(1))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(1),ispline, &
                                     complex_plane_wave_orbitals_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    write (6,*) 'B-splines coefficients successfully set up'

    else
            call die (lhere, 'ndim /= 3 not supported in B-splines yet')
    endif

     call release ('rbc', rbc)
     call release ('complex_plane_wave_orbitals_at_grid_points', &
                    complex_plane_wave_orbitals_at_grid_points)

    endif ! maxval(rkvec_shift) < zero_tolerance
  else ! nkvec > 1

    call alloc ('rbc_complex', rbc_complex, ndim, 2)
    ! calculate number of orbitals in the primitive cell
    ! so phase properly multiplies the orbitals later
    !!! should do this better - check and make sure this is an integer !!!
    !! is this defined elsewhere? does it relate to nband? !!
    norb_primitive=maxval(nband(1:nkvec))
    !real and imaginary terms of orbital value of plane waves at particular r
    call alloc (  'orb_pw_terms', orb_pw_terms, nkvec, norb_primitive, 2)
    !need both real and imaginary parts of plane-wave sum splined separately
    !for each k-point
    nspline=sum(nband(1:nkvec))

    l_need_phase = .true.
    write(6,*)'More than one k-point exists, will use phase'
    write(6,*)'Will spline real and imaginary terms separately'

    if(ndim == 3) then

      allocate (complex_plane_wave_orbitals_at_grid_points(0:ngrid-1,nkvec,norb_primitive))

      !if orbitals_num_bspline exists, read orbitals on grid from there and skip calculating them
      num_orb_exist=0
      open(4,file='orbitals_num_bspline',form='unformatted',status='old',err=90)
      num_orb_exist=1
      write(6,'(''Bspline interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
      read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,igrad_lapf
      if(ngrid_orbxf /= ngrid_orbx .or. ngrid_orbyf /= ngrid_orby .or.  ngrid_orbzf /= ngrid_orbz) then
        write(6,'(''orbital grids do not match, program, orbitals_num:'',3i4,x,3i4)') &
   &    ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
        call die (lhere, 'orbital grids do not match')
      endif
      if(igrad_lapf /= igrad_lap) then
         write(6, '(''choice of interpolations do not match, program, orbitals_num:'',i4,x,i4)') &
              igrad_lap,igrad_lapf
         call die (lhere, 'choice of interpolations do not match; change igrad_lap')
      endif
      read(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
      write(6,'(''Done reading orbitals_num_bspline'')')
      close(4)
      goto 100

      ! Set up an evenly spaced grid for the control points in each direction
      !        with values from 0 to 1
90    do ix=1,ngrid_orb(1)
         do iy=1,ngrid_orb(2)
            do iz=1,ngrid_orb(3)
                grid_orb(1)=(ix-one)/(ngrid_orb(1))
                grid_orb(2)=(iy-one)/(ngrid_orb(2))
                grid_orb(3)=(iz-one)/(ngrid_orb(3))

               ! Convert the current control point into Cartesian coordinates
               ! because the orbitals evaluation routine needs a point in Cartesian
               ! coordinates -- use primitive cell!!
               r(:)=matmul(rlatt,grid_orb(:))

               ! Evaluate grid point on primitive cell
               call orbitals_pw_primitive(r,orb_pw_terms,dorb_pw_terms,ddorb_pw_terms)

               ! Index the array in the way einspline expects
               igrid_bsplines=(((ix-1)*ngrid_orb(2)+(iy-1))*ngrid_orb(3)+(iz-1))

               ! Store the plane-wave value in the array
               complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,:,:)=cmplx(orb_pw_terms(:,:,1),orb_pw_terms(:,:,2))

            enddo !iz
         enddo !iy
      enddo !ix

100 call release ('orb_pw', orb_pw)
    call release ('orb_pw_terms', orb_pw_terms)

    if(index(mode,'mpi') == 0 .or. idtask == 0) then
    !if orbitals_num_bspline does not exist and inum_orb is positive, write them to that file.
       if(num_orb_exist==0) then
         open(4,file='orbitals_num_bspline',form='unformatted',status='new',err=120)
         write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap
         write(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
         write(6,'(''Done writing orbitals_num_bspline'')')
         close(4)
       endif
    endif ! mpi or idtask == 0

    !grid point coordinates run from 0 to 1
120 rnaught(:) = 0.d0
    rone(:) = 1.d0
    !boundary conditions (periodic: ibc=0, antiperiodic: ibc=5)
    ibc = 0
    rbc_complex = (0.,0.)
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(1))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(1),ispline, &
                                     complex_plane_wave_orbitals_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    write (6,*) 'B-splines coefficients successfully set up'
    else
       call die (lhere, 'ndim != 3 not supported in B-splines yet')
    endif

  call release ('rbc_complex', rbc_complex)
  endif ! nkvec == 1

  ! Calculate gradient from orbital interpolation, interpolate Laplacian separately
  case (1)
    call alloc ('ddorb_pw', ddorb_pw,       norb)

    ! Cases (One k-point=Gamma; One k-point=G/2; Multiple k-points)
    if(nkvec == 1) then
      norb_primitive=nband(1)
      !real and imaginary terms of orbital value of plane waves at particular r
      call alloc (  'orb_pw_terms',   orb_pw_terms, nkvec, norb_primitive, 2)
      call alloc ('ddorb_pw_terms', ddorb_pw_terms, nkvec, norb_primitive, 2)
      nspline=norb
      !if Gamma,
      ! then no phase is necessary
      ! and we only need one set of B-spline coefficients for each orbital
      if (maxval(rkvec_shift) < zero_tolerance ) then
        if(ndim == 3) then
          allocate (plane_wave_orbitals_at_grid_points(0:ngrid-1,norb_primitive))
          allocate (plane_wave_orbitals_laplacian_at_grid_points(0:ngrid-1,norb_primitive))
          ! Set up an evenly spaced grid in each direction from 0 to (ngrid-1)/ngrid

          !If oritals_num_bspline exists, read pw orbitals on grid from there and skip calculating them
          num_orb_exist=0
          open(4,file='orbitals_num_bspline',form='unformatted',status='old',err=210)
          num_orb_exist=1
          write(6,'(''Bspline interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
          read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,igrad_lapf
          if(ngrid_orbxf /= ngrid_orbx .or. ngrid_orbyf /= ngrid_orby .or.  ngrid_orbzf /= ngrid_orbz) then
             write(6,'(''orbital grids do not match, program, orbitals_num:'',3i4,x,3i4)')&
                  ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
             call die (lhere, 'orbital grids do not match')
          endif
          if(igrad_lapf /= igrad_lap) then
             write(6, '(''choice of interpolations do not match, program, orbitals_num:'',i4,x,i4)')&
                  igrad_lap,igrad_lapf
             call die (lhere, 'choice of interpolations do not match; change igrad_lap')
          endif
          read(4) ((plane_wave_orbitals_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
          read(4) ((plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
          write(6,'(''Done reading orbitals_num_bspline'')')
          close(4)
          goto 220

 210      do ix=1,ngrid_orb(1)
            do iy=1,ngrid_orb(2)
              do iz=1,ngrid_orb(3)
                grid_orb(1)=(ix-one)/(ngrid_orb(1))
                grid_orb(2)=(iy-one)/(ngrid_orb(2))
                grid_orb(3)=(iz-one)/(ngrid_orb(3))

                ! Convert the current control point into Cartesian coordinates
                ! because the orbitals evaluation routine needs a point in Cartesian
                ! coordinates -- use primitive cell!!
                r(:)=matmul(rlatt,grid_orb(:))

                ! Evaluate grid point on primitive cell
                call orbitals_pw_primitive(r,orb_pw_terms,dorb_pw_terms,ddorb_pw_terms)

                ! Index the array in the way einspline expects
                igrid_bsplines=(((ix-1)*ngrid_orb(2)+(iy-1))*ngrid_orb(3)+(iz-1))

                do iorb=1,norb_primitive
                  if(ireal_imag(iorb) == 1) then !use real part of plane waves
                    orb_pw(iorb)= orb_pw_terms(1,iorb,1)
                    ddorb_pw(iorb)= ddorb_pw_terms(1,iorb,1)
                  else if(ireal_imag(iorb) == 2) then !use imaginary part of plane waves
                    orb_pw(iorb)= orb_pw_terms(1,iorb,2)
                    ddorb_pw(iorb)= ddorb_pw_terms(1,iorb,2)
                  endif
                enddo

                ! Store the plane-wave value and its Laplacian in arrays
                plane_wave_orbitals_at_grid_points(igrid_bsplines,:)=orb_pw
                plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,:)=ddorb_pw

              enddo !iz
            enddo !iy
          enddo !ix

  220     call release ('orb_pw', orb_pw)
          call release ('ddorb_pw', ddorb_pw)
          call release ('orb_pw_terms', orb_pw_terms)
          call release ('ddorb_pw_terms', ddorb_pw_terms)

          if(index(mode,'mpi') == 0 .or. idtask == 0) then

            !if orbitals_num_bspline does not exist and inum_orb is positive, write them to that file.
            if(num_orb_exist == 0) then
              open(4,file='orbitals_num_bspline',form='unformatted',status='new',err=230)
              write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap
              write(4) ((plane_wave_orbitals_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
              write(4) ((plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
              write(6,'(''Done writing orbitals_num_bspline'')')
              close(4)
            endif
         endif ! mpi or idtask == 0

          !grid point coordinates run from 0 to 1
  230     rnaught(:) = 0.d0
          rone(:) = 1.d0
          !boundary conditions (periodic: ibc=0, antiperiodic: ibc=5)
          ibc = 0
          rbc = 0.d0
          !create spline for orbitals
          call fcreate_multi_ubspline_3d_d(rnaught(1),rone(1),ngrid_orb(1), &
                                            rnaught(2),rone(2),ngrid_orb(2), &
                                            rnaught(3),rone(3),ngrid_orb(3), &
                                            ibc(1,1),rbc(1,1),ibc(1,2),rbc(1,2), &
                                            ibc(2,1),rbc(2,1),ibc(2,2),rbc(2,2), &
                                            ibc(3,1),rbc(3,1),ibc(3,2),rbc(3,2), &
                                            nspline,bspline_coefficient_pointer(1))
          do iorb=1,nspline
             ispline=iorb-1
             call fset_multi_ubspline_3d_d(bspline_coefficient_pointer(1), ispline, &
                                           plane_wave_orbitals_at_grid_points(:,iorb))
          enddo
          !create spline for orbitals' Laplacian
          call fcreate_multi_ubspline_3d_d(rnaught(1),rone(1),ngrid_orb(1), &
                                            rnaught(2),rone(2),ngrid_orb(2), &
                                            rnaught(3),rone(3),ngrid_orb(3), &
                                            ibc(1,1),rbc(1,1),ibc(1,2),rbc(1,2), &
                                            ibc(2,1),rbc(2,1),ibc(2,2),rbc(2,2), &
                                            ibc(3,1),rbc(3,1),ibc(3,2),rbc(3,2), &
                                            nspline,bspline_coefficient_pointer(2))
          do iorb=1,nspline
             ispline=iorb-1
             call fset_multi_ubspline_3d_d(bspline_coefficient_pointer(2), ispline, &
                                           plane_wave_orbitals_laplacian_at_grid_points(:,iorb))
          enddo
          write (6,*) 'B-splines coefficients successfully set up'
        endif ! ndim == 3

        call release ('rbc', rbc)
        call release ('plane_wave_orbitals_at_grid_points', plane_wave_orbitals_at_grid_points)
        call release ('plane_wave_orbitals_laplacian_at_grid_points', plane_wave_orbitals_laplacian_at_grid_points)
        !if non-Gamma but single k-point,
        ! then phase is necessary
        ! but we only need one set of B-spline coefficients for each orbital
      else !maxval(rkvec_shift) > zero_tolerance

        l_need_phase = .true.
        write(6,*)'Non-gamma k-point in use, will use phase'
        write(6,*)'Will spline real and imaginary terms separately'
        call alloc ('rbc_complex', rbc_complex, ndim, 2)

        if(ndim == 3) then
        allocate (complex_plane_wave_orbitals_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
        allocate (complex_plane_wave_orbitals_laplacian_at_grid_points(0:ngrid-1,nkvec,norb_primitive))

        !If orbitals_num_bspline exists, read orbitals on grid from there and skip calculating them
        num_orb_exist=0
        open(4,file='orbitals_num_bspline',form='unformatted',status='old',err=250)
        num_orb_exist=1
        write(6,'(''Bspline interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
        read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,igrad_lapf
        if(ngrid_orbxf /= ngrid_orbx .or. ngrid_orbyf /= ngrid_orby .or.  ngrid_orbzf /= ngrid_orbz) then
          write(6,'(''orbital grids do not match, program, orbitals_num:'',3i4,x,3i4)') &
               ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
          call die (lhere, 'orbital grids do not match')
        endif
        if(igrad_lapf /= igrad_lap) then
           write(6, '(''choice of interpolations do not match, program, orbitals_num:'',i4,x,i4)') &
                igrad_lap,igrad_lapf
           call die (lhere, 'choice of interpolations do not match; change igrad_lap')
        endif
        read(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
        read(4) (((complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
          write(6,'(''Done reading orbitals_num_bspline'')')
          close(4)
!       endif
        goto 260

        ! Set up an evenly spaced grid for the control points in each direction
        !        with values from 0 to 1
 250    do ix=1,ngrid_orb(1)
           do iy=1,ngrid_orb(2)
              do iz=1,ngrid_orb(3)
                  grid_orb(1)=(ix-one)/(ngrid_orb(1))
                  grid_orb(2)=(iy-one)/(ngrid_orb(2))
                  grid_orb(3)=(iz-one)/(ngrid_orb(3))

                 ! Convert the current control point into Cartesian coordinates
                 ! because the orbitals evaluation routine needs a point in Cartesian
                 ! coordinates -- use primitive cell!!
                 r(:)=matmul(rlatt,grid_orb(:))

                 ! Evaluate grid point on primitive cell
                 call orbitals_pw_primitive(r,orb_pw_terms,dorb_pw_terms,ddorb_pw_terms)

                 ! Index the array in the way einspline expects
                 igrid_bsplines=(((ix-1)*ngrid_orb(2)+(iy-1))*ngrid_orb(3)+(iz-1))

                 ! Store the plane-wave value in the array
                 complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,:,:)=&
                 &        cmplx(orb_pw_terms(:,:,1),orb_pw_terms(:,:,2))
                 complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,:,:)=&
                 &        cmplx(ddorb_pw_terms(:,:,1),ddorb_pw_terms(:,:,2))

              enddo !iz
           enddo !iy
        enddo !ix

 260    call release ('orb_pw', orb_pw)
        call release ('ddorb_pw', ddorb_pw)
        call release ('orb_pw_terms', orb_pw_terms)
        call release ('ddorb_pw_terms', ddorb_pw_terms)

        if(index(mode,'mpi') == 0 .or. idtask == 0) then

           !if orbitals_num_bspline does not exist and inum_orb is positive, write them to that file.
           if(num_orb_exist==0) then
             open(4,file='orbitals_num_bspline',form='unformatted',status='new',err=270)
             write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap
             write(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
             write(4) (((complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
             write(6,'(''Done writing orbitals_num_bspline'')')
             close(4)
           endif
        endif ! mpi or idtask == 0

    !grid point coordinates run from 0 to 1
 270  rnaught(:) = 0.d0
    rone(:) = 1.d0
    !boundary conditions (periodic: ibc=0, antiperiodic: ibc=5)
    ibc = 0
    rbc_complex = (0.,0.)
    !create splines for orbitals
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(1))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(1),ispline, &
                                     complex_plane_wave_orbitals_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    !create splines for orbitals' Laplacian
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(2))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(2),ispline, &
                                     complex_plane_wave_orbitals_laplacian_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    write (6,*) 'B-splines coefficients successfully set up'

    else
            call die (lhere, 'ndim != 3 not supported in B-splines yet')
    endif

    call release ('rbc', rbc)
    call release ('complex_plane_wave_orbitals_at_grid_points', &
                   complex_plane_wave_orbitals_at_grid_points)
    call release ('complex_plane_wave_orbitals_laplacian_at_grid_points', &
                   complex_plane_wave_orbitals_laplacian_at_grid_points)

    endif ! maxval(rkvec_shift) < zero_tolerance
  else ! nkvec > 1

    call alloc ('rbc_complex', rbc_complex, ndim, 2)
    ! calculate number of orbitals in the primitive cell
    ! so phase properly multiplies the orbitals later
    !!! should do this better - check and make sure this is an integer !!!
    !! is this defined elsewhere? does it relate to nband? !!
    norb_primitive=maxval(nband(1:nkvec))
    !real and imaginary terms of orbital value of plane waves at particular r
    call alloc (  'orb_pw_terms',   orb_pw_terms, nkvec, norb_primitive, 2)
    call alloc ('ddorb_pw_terms', ddorb_pw_terms, nkvec, norb_primitive, 2)
    !need both real and imaginary parts of plane-wave sum splined separately
    !for each k-point
    nspline=sum(nband(1:nkvec))

    l_need_phase = .true.
    write(6,*)'More than one k-point exists, will use phase'
    write(6,*)'Will spline real and imaginary terms separately'

    if(ndim == 3) then

      allocate (complex_plane_wave_orbitals_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
      allocate (complex_plane_wave_orbitals_laplacian_at_grid_points(0:ngrid-1,nkvec,norb_primitive))

!If orbitals_num_bspline exists, read orbitals on grid from there and skip calculating them
      num_orb_exist=0
        open(4,file='orbitals_num_bspline',form='unformatted',status='old',err=290)
        num_orb_exist=1
        write(6,'(''Bspline interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
        read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,igrad_lapf
        if(ngrid_orbxf /= ngrid_orbx .or. ngrid_orbyf /= ngrid_orby .or.  ngrid_orbzf /= ngrid_orbz) then
          write(6,'(''orbital grids do not match, program, orbitals_num:'',3i4,x,3i4)') &
               ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
          call die (lhere, 'orbital grids do not match')
        endif
        if(igrad_lapf /= igrad_lap) then
           write(6, '(''choice of interpolations do not match, program, orbitals_num:'',i4,x,i4)') &
                igrad_lap,igrad_lapf
           call die (lhere, 'choice of interpolations do not match; change igrad_lap')
        endif
        read(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
        read(4) (((complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
        write(6,'(''Done reading orbitals_num_bspline'')')
        close(4)
      goto 300

      ! Set up an evenly spaced grid for the control points in each direction
      !        with values from 0 to 1
290   do ix=1,ngrid_orb(1)
         do iy=1,ngrid_orb(2)
            do iz=1,ngrid_orb(3)
                grid_orb(1)=(ix-one)/(ngrid_orb(1))
                grid_orb(2)=(iy-one)/(ngrid_orb(2))
                grid_orb(3)=(iz-one)/(ngrid_orb(3))

               ! Convert the current control point into Cartesian coordinates
               ! because the orbitals evaluation routine needs a point in Cartesian
               ! coordinates -- use primitive cell!!
               r(:)=matmul(rlatt,grid_orb(:))

               ! Evaluate grid point on primitive cell
               call orbitals_pw_primitive(r,orb_pw_terms,dorb_pw_terms,ddorb_pw_terms)

               ! Index the array in the way einspline expects
               igrid_bsplines=(((ix-1)*ngrid_orb(2)+(iy-1))*ngrid_orb(3)+(iz-1))

               ! Store the plane-wave value in the array
               complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,:,:)=cmplx(orb_pw_terms(:,:,1),orb_pw_terms(:,:,2))
               complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,:,:)=cmplx(ddorb_pw_terms(:,:,1),ddorb_pw_terms(:,:,2))

            enddo !iz
         enddo !iy
      enddo !ix

  300 call release (  'orb_pw',   orb_pw)
      call release ('ddorb_pw', ddorb_pw)
      call release (  'orb_pw_terms',   orb_pw_terms)
      call release ('ddorb_pw_terms', ddorb_pw_terms)

      if(index(mode,'mpi') == 0 .or. idtask == 0) then
         !if orbitals_num_bspline does not exist and inum_orb is positive, write them to that file.
         if(num_orb_exist == 0) then
           open(4,file='orbitals_num_bspline',form='unformatted',status='new',err=320)
           write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap
           write(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
           write(4) (((complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
           write(6,'(''Done writing orbitals_num_bspline'')')
           close(4)
         endif
      endif ! mpi or idtask == 0

      !grid point coordinates run from 0 to 1
  320 rnaught(:) = 0.d0
      rone(:) = 1.d0
      !boundary conditions (periodic: ibc=0, antiperiodic: ibc=5)
      ibc = 0
      rbc_complex = (0.,0.)
      !create splines for orbitals
      call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                       rnaught(2),rone(2),ngrid_orb(2), &
                                       rnaught(3),rone(3),ngrid_orb(3), &
                                       ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                       ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                       ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                       nspline, &
                                       bspline_coefficient_pointer(1))
      do ikvec=1,nkvec
         do iorb=1,nband(ikvec)
            ispline=(ikvec-1)*nkvec+(iorb-1)
            call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(1),ispline, &
                                       complex_plane_wave_orbitals_at_grid_points(:,ikvec,iorb))
         enddo
      enddo
      !create splines for orbitals' Laplacian
      call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                       rnaught(2),rone(2),ngrid_orb(2), &
                                       rnaught(3),rone(3),ngrid_orb(3), &
                                       ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                       ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                       ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                       nspline, &
                                       bspline_coefficient_pointer(2))
      do ikvec=1,nkvec
         do iorb=1,nband(ikvec)
            ispline=(ikvec-1)*nkvec+(iorb-1)
            call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(2),ispline, &
                                       complex_plane_wave_orbitals_laplacian_at_grid_points(:,ikvec,iorb))
         enddo
      enddo
      write (6,*) 'B-splines coefficients successfully set up'
      else
         call die (lhere, 'ndim != 3 not supported in B-splines yet')
      endif

    call release ('rbc_complex', rbc_complex)
    call release ('complex_plane_wave_orbitals_at_grid_points', &
                   complex_plane_wave_orbitals_at_grid_points)
    call release ('complex_plane_wave_orbitals_laplacian_at_grid_points', &
                   complex_plane_wave_orbitals_laplacian_at_grid_points)
    endif ! nkvec == 1

  ! Calculate orbitals, gradient and Laplacian separately
  case (2)
    call alloc ( 'dorb_pw', dorb_pw, ndim, norb)
    call alloc ('ddorb_pw',ddorb_pw,       norb)

    ! Cases (One k-point=Gamma; One k-point=G/2; Multiple k-points)
    if(nkvec == 1) then
      norb_primitive=nband(1)
      !real and imaginary terms of orbital value of plane waves at particular r
      call alloc (  'orb_pw_terms',   orb_pw_terms, nkvec,       norb_primitive, 2)
      call alloc ( 'dorb_pw_terms',  dorb_pw_terms, nkvec, ndim, norb_primitive, 2)
      call alloc ('ddorb_pw_terms', ddorb_pw_terms, nkvec,       norb_primitive, 2)
      nspline=norb
      !if Gamma,
      ! then no phase is necessary
      ! and we only need one set of B-spline coefficients for each orbital
      if (maxval(rkvec_shift) < zero_tolerance ) then
        if(ndim == 3) then
          allocate (plane_wave_orbitals_at_grid_points(0:ngrid-1,norb_primitive))
          allocate (plane_wave_orbitals_gradient1_at_grid_points(0:ngrid-1,norb_primitive))
          allocate (plane_wave_orbitals_gradient2_at_grid_points(0:ngrid-1,norb_primitive))
          allocate (plane_wave_orbitals_gradient3_at_grid_points(0:ngrid-1,norb_primitive))
          allocate (plane_wave_orbitals_laplacian_at_grid_points(0:ngrid-1,norb_primitive))
          !set up an evenly spaced grid in each direction from 0 to (ngrid-1)/ngrid

          !if oritals_num_bspline exists, read pw orbitals on grid from there and skip calculating them
          num_orb_exist=0
          open(4,file='orbitals_num_bspline',form='unformatted',status='old',err=410)
          num_orb_exist=1
          write(6,'(''Bspline interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
          read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,igrad_lapf
          if(ngrid_orbxf /= ngrid_orbx .or. ngrid_orbyf /= ngrid_orby .or.  ngrid_orbzf /= ngrid_orbz) then
             write(6,'(''orbital grids do not match, program, orbitals_num:'',3i4,x,3i4)') &
                  ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
             call die (lhere, 'orbital grids do not match; change ngrid_orbx, ngrid_orby and/or ngrid_orbz')
          endif
          if(igrad_lapf /= igrad_lap) then
             write(6, '(''choice of interpolations do not match, program, orbitals_num:'',i4,x,i4)') &
                  igrad_lap,igrad_lapf
             call die (lhere, 'choice of interpolations do not match; change igrad_lap')
          endif
          read(4) ((plane_wave_orbitals_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
          read(4) ((plane_wave_orbitals_gradient1_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
          read(4) ((plane_wave_orbitals_gradient2_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
          read(4) ((plane_wave_orbitals_gradient3_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
          read(4) ((plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
          write(6,'(''Done reading orbitals_num_bspline'')')
          close(4)
          goto 420

 410      do ix=1,ngrid_orb(1)
            do iy=1,ngrid_orb(2)
              do iz=1,ngrid_orb(3)
                grid_orb(1)=(ix-one)/(ngrid_orb(1))
                grid_orb(2)=(iy-one)/(ngrid_orb(2))
                grid_orb(3)=(iz-one)/(ngrid_orb(3))

                ! Convert the current control point into Cartesian coordinates
                ! because the orbitals evaluation routine needs a point in Cartesian
                ! coordinates -- use primitive cell!!
                r(:)=matmul(rlatt,grid_orb(:))

                ! Evaluate grid point on primitive cell
                call orbitals_pw_primitive(r,orb_pw_terms,dorb_pw_terms,ddorb_pw_terms)

                ! Index the array in the way einspline expects
                igrid_bsplines=(((ix-1)*ngrid_orb(2)+(iy-1))*ngrid_orb(3)+(iz-1))

                do iorb=1,norb_primitive
                  if(ireal_imag(iorb) == 1) then !use real part of plane waves
                    orb_pw(iorb)    =   orb_pw_terms(1,  iorb,1)
                    dorb_pw(:,iorb) =  dorb_pw_terms(1,:,iorb,1)
                    ddorb_pw(iorb)  = ddorb_pw_terms(1,  iorb,1)
                  else if(ireal_imag(iorb) == 2) then !use imaginary part of plane waves
                    orb_pw(iorb)    =   orb_pw_terms(1,  iorb,2)
                    dorb_pw(:,iorb) =  dorb_pw_terms(1,:,iorb,2)
                    ddorb_pw(iorb)  = ddorb_pw_terms(1,  iorb,2)
                  endif
                enddo

                ! Store the plane-wave value in the array
                plane_wave_orbitals_at_grid_points(igrid_bsplines,:)=orb_pw
                plane_wave_orbitals_gradient1_at_grid_points(igrid_bsplines,:)=dorb_pw(1,:)
                plane_wave_orbitals_gradient2_at_grid_points(igrid_bsplines,:)=dorb_pw(2,:)
                plane_wave_orbitals_gradient3_at_grid_points(igrid_bsplines,:)=dorb_pw(3,:)
                plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,:)=ddorb_pw

              enddo !iz
            enddo !iy
          enddo !ix

  420     call release (  'orb_pw',   orb_pw)
          call release ( 'dorb_pw',  dorb_pw)
          call release ('ddorb_pw', ddorb_pw)
          call release (  'orb_pw_terms',   orb_pw_terms)
          call release ( 'dorb_pw_terms',  dorb_pw_terms)
          call release ('ddorb_pw_terms', ddorb_pw_terms)

          if(index(mode,'mpi') == 0 .or. idtask == 0) then
            !if orbitals_num_bspline does not exist and inum_orb is positive, write them to that file.
            if(num_orb_exist == 0) then
              open(4,file='orbitals_num_bspline',form='unformatted',status='new',err=430)
              write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap
              write(4) ((plane_wave_orbitals_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
              write(4) ((plane_wave_orbitals_gradient1_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
              write(4) ((plane_wave_orbitals_gradient2_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
              write(4) ((plane_wave_orbitals_gradient3_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
              write(4) ((plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,iorb),igrid_bsplines=0,ngrid-1),iorb=1,norb)
              write(6,'(''Done writing orbitals_num_bspline'')')
              close(4)
            endif
          endif ! mpi or idtask == 0

          !grid point coordinates run from 0 to 1
  430     rnaught(:) = 0.d0
          rone(:) = 1.d0
          !boundary conditions (periodic: ibc=0, antiperiodic: ibc=5)
          ibc = 0
          rbc = 0.d0
          !create spline for orbitals
          call fcreate_multi_ubspline_3d_d(rnaught(1),rone(1),ngrid_orb(1), &
                                            rnaught(2),rone(2),ngrid_orb(2), &
                                            rnaught(3),rone(3),ngrid_orb(3), &
                                            ibc(1,1),rbc(1,1),ibc(1,2),rbc(1,2), &
                                            ibc(2,1),rbc(2,1),ibc(2,2),rbc(2,2), &
                                            ibc(3,1),rbc(3,1),ibc(3,2),rbc(3,2), &
                                            nspline,bspline_coefficient_pointer(1))
          do iorb=1,nspline
             ispline=iorb-1
             call fset_multi_ubspline_3d_d(bspline_coefficient_pointer(1), ispline, &
                                           plane_wave_orbitals_at_grid_points(:,iorb))
          enddo
          !create spline for orbitals' Laplacian
          call fcreate_multi_ubspline_3d_d(rnaught(1),rone(1),ngrid_orb(1), &
                                            rnaught(2),rone(2),ngrid_orb(2), &
                                            rnaught(3),rone(3),ngrid_orb(3), &
                                            ibc(1,1),rbc(1,1),ibc(1,2),rbc(1,2), &
                                            ibc(2,1),rbc(2,1),ibc(2,2),rbc(2,2), &
                                            ibc(3,1),rbc(3,1),ibc(3,2),rbc(3,2), &
                                            nspline,bspline_coefficient_pointer(2))
          do iorb=1,nspline
             ispline=iorb-1
             call fset_multi_ubspline_3d_d(bspline_coefficient_pointer(2), ispline, &
                                           plane_wave_orbitals_laplacian_at_grid_points(:,iorb))
          enddo
          !create spline for first component of orbitals' gradient
          call fcreate_multi_ubspline_3d_d(rnaught(1),rone(1),ngrid_orb(1), &
                                            rnaught(2),rone(2),ngrid_orb(2), &
                                            rnaught(3),rone(3),ngrid_orb(3), &
                                            ibc(1,1),rbc(1,1),ibc(1,2),rbc(1,2), &
                                            ibc(2,1),rbc(2,1),ibc(2,2),rbc(2,2), &
                                            ibc(3,1),rbc(3,1),ibc(3,2),rbc(3,2), &
                                            nspline,bspline_coefficient_pointer(3))
          do iorb=1,nspline
             ispline=iorb-1
             call fset_multi_ubspline_3d_d(bspline_coefficient_pointer(3), ispline, &
                                           plane_wave_orbitals_gradient1_at_grid_points(:,iorb))
          enddo
          !create spline for second component of orbitals' gradient
          call fcreate_multi_ubspline_3d_d(rnaught(1),rone(1),ngrid_orb(1), &
                                            rnaught(2),rone(2),ngrid_orb(2), &
                                            rnaught(3),rone(3),ngrid_orb(3), &
                                            ibc(1,1),rbc(1,1),ibc(1,2),rbc(1,2), &
                                            ibc(2,1),rbc(2,1),ibc(2,2),rbc(2,2), &
                                            ibc(3,1),rbc(3,1),ibc(3,2),rbc(3,2), &
                                            nspline,bspline_coefficient_pointer(4))
          do iorb=1,nspline
             ispline=iorb-1
             call fset_multi_ubspline_3d_d(bspline_coefficient_pointer(4), ispline, &
                                           plane_wave_orbitals_gradient2_at_grid_points(:,iorb))
          enddo
          !create spline for third component of orbitals' gradient
          call fcreate_multi_ubspline_3d_d(rnaught(1),rone(1),ngrid_orb(1), &
                                            rnaught(2),rone(2),ngrid_orb(2), &
                                            rnaught(3),rone(3),ngrid_orb(3), &
                                            ibc(1,1),rbc(1,1),ibc(1,2),rbc(1,2), &
                                            ibc(2,1),rbc(2,1),ibc(2,2),rbc(2,2), &
                                            ibc(3,1),rbc(3,1),ibc(3,2),rbc(3,2), &
                                            nspline,bspline_coefficient_pointer(5))
          do iorb=1,nspline
             ispline=iorb-1
             call fset_multi_ubspline_3d_d(bspline_coefficient_pointer(5), ispline, &
                                           plane_wave_orbitals_gradient3_at_grid_points(:,iorb))
          enddo
          write (6,*) 'B-splines coefficients successfully set up'

        call release ('plane_wave_orbitals_gradient1_at_grid_points', plane_wave_orbitals_gradient1_at_grid_points)
        call release ('plane_wave_orbitals_gradient2_at_grid_points', plane_wave_orbitals_gradient2_at_grid_points)
        call release ('plane_wave_orbitals_gradient3_at_grid_points', plane_wave_orbitals_gradient3_at_grid_points)

        endif ! ndim == 3

        call release ('rbc', rbc)
        call release ('plane_wave_orbitals_at_grid_points', plane_wave_orbitals_at_grid_points)
        call release ('plane_wave_orbitals_laplacian_at_grid_points', plane_wave_orbitals_laplacian_at_grid_points)

        !if non-Gamma but single k-point,
        ! then phase is necessary
        ! but we only need one set of B-spline coefficients for each orbital
      else !maxval(rkvec_shift) > zero_tolerance

        l_need_phase = .true.
        write(6,*)'Non-gamma k-point in use, will use phase'
        write(6,*)'Will spline real and imaginary terms separately'
        call alloc ('rbc_complex', rbc_complex, ndim, 2)

        if(ndim == 3) then
        allocate (complex_plane_wave_orbitals_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
        allocate (complex_plane_wave_orbitals_laplacian_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
        allocate (complex_plane_wave_orbitals_gradient1_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
        allocate (complex_plane_wave_orbitals_gradient2_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
        allocate (complex_plane_wave_orbitals_gradient3_at_grid_points(0:ngrid-1,nkvec,norb_primitive))

        !if orbitals_num_bspline exists, read orbitals on grid from there and skip calculating them
        num_orb_exist=0
        open(4,file='orbitals_num_bspline',form='unformatted',status='old',err=450)
        num_orb_exist=1
        write(6,'(''Bspline interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
        read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,igrad_lapf
        if(ngrid_orbxf /= ngrid_orbx .or. ngrid_orbyf /= ngrid_orby .or.  ngrid_orbzf /= ngrid_orbz) then
          write(6,'(''orbital grids do not match, program, orbitals_num:'',3i4,x,3i4)') &
               ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
          call die (lhere, 'orbital grids do not match; change ngrid_orbx, ngrid_orby and/or ngrid_orbz')
        endif
        if(igrad_lapf /= igrad_lap) then
           write(6, '(''choice of interpolations do not match, program, orbitals_num:'',i4,x,i4)') &
                igrad_lap,igrad_lapf
           call die (lhere, 'choice of interpolations do not match; change igrad_lap')
        endif
        read(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
        read(4) (((complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
        read(4) (((complex_plane_wave_orbitals_gradient1_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
        read(4) (((complex_plane_wave_orbitals_gradient2_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
        read(4) (((complex_plane_wave_orbitals_gradient3_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
        write(6,'(''Done reading orbitals_num_bspline'')')
        close(4)
!       endif
        goto 460

        ! Set up an evenly spaced grid for the control points in each direction
        !        with values from 0 to 1
 450    do ix=1,ngrid_orb(1)
           do iy=1,ngrid_orb(2)
              do iz=1,ngrid_orb(3)
                  grid_orb(1)=(ix-one)/(ngrid_orb(1))
                  grid_orb(2)=(iy-one)/(ngrid_orb(2))
                  grid_orb(3)=(iz-one)/(ngrid_orb(3))

                 ! Convert the current control point into Cartesian coordinates
                 ! because the orbitals evaluation routine needs a point in Cartesian
                 ! coordinates -- use primitive cell!!
                 r(:)=matmul(rlatt,grid_orb(:))

                 ! Evaluate grid point on primitive cell
                 call orbitals_pw_primitive(r,orb_pw_terms,dorb_pw_terms,ddorb_pw_terms)

                 ! Index the array in the way einspline expects
                 igrid_bsplines=(((ix-1)*ngrid_orb(2)+(iy-1))*ngrid_orb(3)+(iz-1))

                 ! Store the plane-wave value in the array
                 complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,:,:)=&
                 &        cmplx(  orb_pw_terms(:,  :,1),  orb_pw_terms(:,  :,2))
                 complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,:,:)=&
                 &        cmplx(ddorb_pw_terms(:,  :,1),ddorb_pw_terms(:,  :,2))
                 complex_plane_wave_orbitals_gradient1_at_grid_points(igrid_bsplines,:,:)=&
                 &        cmplx( dorb_pw_terms(:,1,:,1), dorb_pw_terms(:,1,:,2))
                 complex_plane_wave_orbitals_gradient2_at_grid_points(igrid_bsplines,:,:)=&
                 &        cmplx( dorb_pw_terms(:,2,:,1), dorb_pw_terms(:,2,:,2))
                 complex_plane_wave_orbitals_gradient3_at_grid_points(igrid_bsplines,:,:)=&
                 &        cmplx( dorb_pw_terms(:,3,:,1), dorb_pw_terms(:,3,:,2))

              enddo !iz
           enddo !iy
        enddo !ix

 460    call release ('orb_pw_terms', orb_pw_terms)
        call release ('ddorb_pw_terms', ddorb_pw_terms)
        call release ('dorb_pw_terms', dorb_pw_terms)

        if(index(mode,'mpi') == 0 .or. idtask == 0) then
          !if orbitals_num_bspline does not exist and inum_orb is positive, write them to that file.
          if(num_orb_exist==0) then
            open(4,file='orbitals_num_bspline',form='unformatted',status='new',err=470)
            write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap
            write(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
            write(4) (((complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
            write(4) (((complex_plane_wave_orbitals_gradient1_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
            write(4) (((complex_plane_wave_orbitals_gradient2_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
            write(4) (((complex_plane_wave_orbitals_gradient3_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb)
            write(6,'(''Done writing orbitals_num_bspline'')')
            close(4)
          endif
       endif !mpi or idtask == 0

    !grid point coordinates run from 0 to 1
 470  rnaught(:) = 0.d0
    rone(:) = 1.d0
    !boundary conditions (periodic: ibc=0, antiperiodic: ibc=5)
    ibc = 0
    rbc_complex = (0.,0.)
    !create splines for orbitals
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(1))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(1),ispline, &
                                     complex_plane_wave_orbitals_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    !create splines for orbitals' Laplacian
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(2))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(2),ispline, &
                                     complex_plane_wave_orbitals_laplacian_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    !create splines for first component of orbitals' gradient
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(3))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(3),ispline, &
                                     complex_plane_wave_orbitals_gradient1_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    !create splines for second component of orbitals' gradient
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(4))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(4),ispline, &
                                     complex_plane_wave_orbitals_gradient2_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    !create splines for third component of orbitals' gradient
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(5))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(5),ispline, &
                                     complex_plane_wave_orbitals_gradient3_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    write (6,*) 'B-splines coefficients successfully set up'

    else
            call die (lhere, 'ndim != 3 not supported in B-splines yet')
    endif

     call release ('rbc_complex', rbc_complex)
     call release ('complex_plane_wave_orbitals_at_grid_points', &
                    complex_plane_wave_orbitals_at_grid_points)
     call release ('complex_plane_wave_orbitals_laplacian_at_grid_points', &
                    complex_plane_wave_orbitals_laplacian_at_grid_points)
     call release ('complex_plane_wave_orbitals_gradient1_at_grid_points', &
                    complex_plane_wave_orbitals_gradient1_at_grid_points)
     call release ('complex_plane_wave_orbitals_gradient2_at_grid_points', &
                    complex_plane_wave_orbitals_gradient2_at_grid_points)
     call release ('complex_plane_wave_orbitals_gradient3_at_grid_points', &
                    complex_plane_wave_orbitals_gradient3_at_grid_points)

    endif ! maxval(rkvec_shift) < zero_tolerance
  else ! nkvec > 1

    call alloc ('rbc_complex', rbc_complex, ndim, 2)
    ! calculate number of orbitals in the primitive cell
    ! so phase properly multiplies the orbitals later
    norb_primitive=maxval(nband(1:nkvec))
    !need both real and imaginary parts of plane-wave sum splined separately
    !for each k-point
    nspline=sum(nband(1:nkvec))

    !real and imaginary terms of orbital value of plane waves at particular r
    call alloc (  'orb_pw_terms',   orb_pw_terms, nkvec,       norb_primitive, 2)
    call alloc ( 'dorb_pw_terms',  dorb_pw_terms, nkvec, ndim, norb_primitive, 2)
    call alloc ('ddorb_pw_terms', ddorb_pw_terms, nkvec,       norb_primitive, 2)

    l_need_phase = .true.
    write(6,*)'More than one k-point exists, will use phase'
    write(6,*)'Will spline real and imaginary terms separately'

    if(ndim == 3) then

      allocate (complex_plane_wave_orbitals_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
      allocate (complex_plane_wave_orbitals_laplacian_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
      allocate (complex_plane_wave_orbitals_gradient1_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
      allocate (complex_plane_wave_orbitals_gradient2_at_grid_points(0:ngrid-1,nkvec,norb_primitive))
      allocate (complex_plane_wave_orbitals_gradient3_at_grid_points(0:ngrid-1,nkvec,norb_primitive))

      !if orbitals_num_bspline exists, read orbitals on grid from there and skip calculating them
      num_orb_exist=0
      open(4,file='orbitals_num_bspline',form='unformatted',status='old',err=490)
      num_orb_exist=1
      write(6,'(''Bspline interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
      read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,igrad_lapf
      if((ngrid_orbxf /= ngrid_orbx) .or. (ngrid_orbyf /= ngrid_orby) .or. (ngrid_orbzf /= ngrid_orbz)) then
        write(6,'(''orbital grids do not match, program, orbitals_num:'',3i4,x,3i4)') &
   &    ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
        call die (lhere, 'orbital grids do not match; change ngrid_orbx, ngrid_orby and/or ngrid_orbz')
      endif
      if(igrad_lapf /= igrad_lap) then
         write(6, '(''choice of interpolations do not match, program, orbitals_num:'',i4,x,i4)') &
              igrad_lap,igrad_lapf
         call die (lhere, 'choice of interpolations do not match; change igrad_lap')
      endif
      read(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
      read(4) (((complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
      read(4) (((complex_plane_wave_orbitals_gradient1_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
      read(4) (((complex_plane_wave_orbitals_gradient2_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
      read(4) (((complex_plane_wave_orbitals_gradient3_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
      write(6,'(''Done reading orbitals_num_bspline'')')
      close(4)
      goto 500

      ! Set up an evenly spaced grid for the control points in each direction
      !        with values from 0 to 1
490   do ix=1,ngrid_orb(1)
         do iy=1,ngrid_orb(2)
            do iz=1,ngrid_orb(3)
                grid_orb(1)=(ix-one)/(ngrid_orb(1))
                grid_orb(2)=(iy-one)/(ngrid_orb(2))
                grid_orb(3)=(iz-one)/(ngrid_orb(3))

               ! Convert the current control point into Cartesian coordinates
               ! because the orbitals evaluation routine needs a point in Cartesian
               ! coordinates -- use primitive cell!!
               r(:)=matmul(rlatt,grid_orb(:))

               ! Evaluate grid point on primitive cell
               call orbitals_pw_primitive(r,orb_pw_terms,dorb_pw_terms,ddorb_pw_terms)

               ! Index the array in the way einspline expects
               igrid_bsplines=(((ix-1)*ngrid_orb(2)+(iy-1))*ngrid_orb(3)+(iz-1))

               ! Store the plane-wave value in the array
               complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,:,:)=&
                             cmplx(  orb_pw_terms(:,  :,1),  orb_pw_terms(:,  :,2))
               complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,:,:)=&
                             cmplx(ddorb_pw_terms(:,  :,1),ddorb_pw_terms(:,  :,2))
               complex_plane_wave_orbitals_gradient1_at_grid_points(igrid_bsplines,:,:)=&
                             cmplx( dorb_pw_terms(:,1,:,1), dorb_pw_terms(:,1,:,2))
               complex_plane_wave_orbitals_gradient2_at_grid_points(igrid_bsplines,:,:)=&
                             cmplx( dorb_pw_terms(:,2,:,1), dorb_pw_terms(:,2,:,2))
               complex_plane_wave_orbitals_gradient3_at_grid_points(igrid_bsplines,:,:)=&
                             cmplx( dorb_pw_terms(:,3,:,1), dorb_pw_terms(:,3,:,2))

            enddo !iz
         enddo !iy
      enddo !ix

500 call release (  'orb_pw_terms',   orb_pw_terms)
    call release ('ddorb_pw_terms', ddorb_pw_terms)
    call release ( 'dorb_pw_terms',  dorb_pw_terms)

    if(index(mode,'mpi') == 0 .or. idtask == 0) then
      !if orbitals_num_bspline does not exist and inum_orb is positive, write them to that file.
      if(num_orb_exist == 0) then
        open(4,file='orbitals_num_bspline',form='unformatted',status='new',err=520)
        write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap
        write(4) (((complex_plane_wave_orbitals_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
        write(4) (((complex_plane_wave_orbitals_laplacian_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
        write(4) (((complex_plane_wave_orbitals_gradient1_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
        write(4) (((complex_plane_wave_orbitals_gradient2_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
        write(4) (((complex_plane_wave_orbitals_gradient3_at_grid_points(igrid_bsplines,ikvec,iorb),igrid_bsplines=0,ngrid-1),ikvec=1,nkvec),iorb=1,norb_primitive)
        write(6,'(''Done writing orbitals_num_bspline'')')
        close(4)
      endif
    endif !mpi or idtask == 0

    !grid point coordinates run from 0 to 1
520 rnaught(:) = 0.d0
    rone(:) = 1.d0
    !boundary conditions (periodic: ibc=0, antiperiodic: ibc=5)
    ibc = 0
    rbc_complex = (0.,0.)
    !create spline for orbitals
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(1))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(1),ispline, &
                                     complex_plane_wave_orbitals_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    !create spline for orbitals' Laplacian
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(2))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(2),ispline, &
                                     complex_plane_wave_orbitals_laplacian_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    !create spline for first component of orbitals' gradient
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(3))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(3),ispline, &
                                     complex_plane_wave_orbitals_gradient1_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    !create spline for second component of orbitals' gradient
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(4))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(4),ispline, &
                                     complex_plane_wave_orbitals_gradient2_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    !create spline for third component of orbitals' gradient
    call fcreate_multi_ubspline_3d_z(rnaught(1),rone(1),ngrid_orb(1), &
                                     rnaught(2),rone(2),ngrid_orb(2), &
                                     rnaught(3),rone(3),ngrid_orb(3), &
                                     ibc(1,1),rbc_complex(1,1),ibc(1,2),rbc_complex(1,2), &
                                     ibc(2,1),rbc_complex(2,1),ibc(2,2),rbc_complex(2,2), &
                                     ibc(3,1),rbc_complex(3,1),ibc(3,2),rbc_complex(3,2), &
                                     nspline, &
                                     bspline_coefficient_pointer(5))
    do ikvec=1,nkvec
       do iorb=1,nband(ikvec)
          ispline=(ikvec-1)*nkvec+(iorb-1)
          call fset_multi_ubspline_3d_z(bspline_coefficient_pointer(5),ispline, &
                                     complex_plane_wave_orbitals_gradient3_at_grid_points(:,ikvec,iorb))
       enddo
    enddo
    write (6,*) 'B-splines coefficients successfully set up'
    else
       call die (lhere, 'ndim /= 3 not supported in B-splines yet')
    endif

  call release ('rbc_complex', rbc_complex)
  call release ('complex_plane_wave_orbitals_at_grid_points',&
                 complex_plane_wave_orbitals_at_grid_points)
  call release ('complex_plane_wave_orbitals_laplacian_at_grid_points',&
                 complex_plane_wave_orbitals_laplacian_at_grid_points)
  call release ('complex_plane_wave_orbitals_gradient1_at_grid_points',&
                 complex_plane_wave_orbitals_gradient1_at_grid_points)
  call release ('complex_plane_wave_orbitals_gradient2_at_grid_points',&
                 complex_plane_wave_orbitals_gradient2_at_grid_points)
  call release ('complex_plane_wave_orbitals_gradient3_at_grid_points',&
                 complex_plane_wave_orbitals_gradient3_at_grid_points)
  endif ! nkvec == 1
  case default
       call die (lhere, 'igrad_lap must be 0, 1 or 2')
  end select

  call release ('r', r)
  call release ('rnaught', rnaught)
  call release ('rone', rone)
  call release ('ngrid_orb', ngrid_orb)
  call release ('grid_orb', grid_orb)
  call release ('ibc', ibc)


  ! get primitive cell lattice vector matrix transpose
  rlatt_inv_transpose=transpose(rlatt_inv)

  ! get primitive cell lattice vector matrix transpose times 2
  rlatt_inv_transpose_2=matmul(rlatt_inv,rlatt_inv_transpose)

  end subroutine setup_bsplines_coefficients

! ==============================================================================
  subroutine evaluate_bsplines_function_only(r_cart,orb_bsplines)
! ------------------------------------------------------------------------------
! Description   : Evaluate B-splines for one electron in all orbitals
! Description   : without any derivatives
!
! Created       : W. Parker, 21 Oct 2008
! ------------------------------------------------------------------------------
  implicit none
  !need rlatt_inv,nkvec,rkvec,rkvec_shift,nband
  include 'commons.h'

! input
  real(dp),intent(in)                 :: r_cart(ndim)

!! output
  real(dp),intent(out)                :: orb_bsplines(norb)

! local
  integer                             :: i
  integer                             :: iorb
  integer                             :: ikvec
  integer                             :: ispline
  real(dp)                            :: r(ndim)
  real(dp)                            :: orb_tmp(nspline)
  real(dp)                            :: k_dot_r
  real(dp)                            :: real_phase
  real(dp)                            :: imaginary_phase
  real(dp)                            :: cos_k(nkvec)
  real(dp)                            :: sin_k(nkvec)
  real(dp)                            :: dcos_k(ndim,nkvec)
  real(dp)                            :: dsin_k(ndim,nkvec)
  real(dp)                            :: ddcos_k(nkvec)
  real(dp)                            :: ddsin_k(nkvec)
  complex(dp)                         :: complex_phase
  complex(dp)                         :: orb_tmp_complex(nspline)

! get primitive cell coordinates
  r=matmul(rlatt_inv,r_cart)
  !make sure values are between zero and one
  r=mod(r+abs(int(r))+one,one)

  if(l_need_phase) then
    call feval_multi_ubspline_3d_z(bspline_coefficient_pointer(1), &
                                   r(1), r(2), r(3), orb_tmp_complex)
    iorb=0
    call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,r_cart,ng1d_sim, &
                     cos_k,sin_k,dcos_k,dsin_k,ddcos_k,ddsin_k,rkvec_shift)
    do ikvec=1,nkvec
       complex_phase=cmplx(cos_k(ikvec),sin_k(ikvec))
       do i=1,nband(ikvec)
          iorb=iorb+1
          ispline=(ikvec-1)*nkvec+i
          if(ireal_imag(iorb) == 1) then !use real part
            orb_bsplines(iorb)   =dble(orb_tmp_complex(ispline)*complex_phase)
          else if(ireal_imag(iorb) == 2) then !use imaginary part
            orb_bsplines(iorb)   =aimag(orb_tmp_complex(ispline)*complex_phase)
          else !use both real and imaginary parts
            orb_bsplines(iorb)   =dble(orb_tmp_complex(ispline)*complex_phase)
            iorb=iorb+1
            orb_bsplines(iorb)   =aimag(orb_tmp_complex(ispline)*complex_phase)
          endif
       enddo !orbitals
    enddo !k-points
  else !no phase (Gamma)
    call feval_multi_ubspline_3d_d(bspline_coefficient_pointer(1), &
                                   r(1), r(2), r(3), orb_tmp)
    orb_bsplines=orb_tmp
  endif !need phase

  end subroutine evaluate_bsplines_function_only

! ==============================================================================
  subroutine evaluate_bsplines_with_derivatives(r_cart,orb_bsplines,dorb_bsplines,ddorb_bsplines)
! ------------------------------------------------------------------------------
! Description : Evaluate B-splines for one electron in all orbitals
! Description : returning gradient and Laplacian
!
! Created     : W. Parker, 21 Oct 2008
! ------------------------------------------------------------------------------

  implicit none
  !need ndim,rlatt_inv,rkvec,nkvec,glatt_sim,rknorm,rkvec,kvec,ng1d_sim,rkvec_shift,ireal_imag
  include 'commons.h'

! input
  real(dp), intent(in)                :: r_cart(ndim)

! output
  real(dp), intent(out)                :: orb_bsplines(norb)
  real(dp), intent(out)                :: dorb_bsplines(ndim,norb)
  real(dp), intent(out)                :: ddorb_bsplines(norb)

! local
  integer                             :: i
  integer                             :: iorb
  integer                             :: ikvec
  integer                             :: ispline
  real(dp)                            :: r(ndim)
  real(dp)                            :: k_dot_r
  real(dp)                            :: orb_tmp(nspline)
  real(dp)                            :: dorb_tmp(ndim,nspline)
  real(dp)                            :: hessian_bsplines(ndim*ndim,nspline)
  real(dp)                            :: hessian_tmp(ndim,ndim)
  real(dp)                            :: ddorb_tmp(nspline)
  real(dp)                            :: cos_k(nkvec)
  real(dp)                            :: sin_k(nkvec)
  real(dp)                            :: dcos_k(ndim,nkvec)
  real(dp)                            :: dsin_k(ndim,nkvec)
  real(dp)                            :: ddcos_k(nkvec)
  real(dp)                            :: ddsin_k(nkvec)
  complex(dp)                         :: orb_tmp_complex(nspline)
  complex(dp)                         :: dorb_tmp_complex(ndim,nspline)
  complex(dp)                         :: ddorb_tmp_complex(nspline)
  complex(dp)                         :: hessian_bsplines_complex(ndim*ndim,nspline)
  complex(dp)                         :: complex_phase
  complex(dp)                         :: dcomplex_phase(ndim)
  complex(dp)                         :: ddcomplex_phase
  complex(dp)                         :: hessian_tmp_complex(ndim,ndim)

! get primitive cell coordinates
  r=matmul(rlatt_inv,r_cart)
!make sure values are between zero and one
  r=mod(r+abs(int(r))+one,one)

  iorb=0
  select case (igrad_lap)
  case (0) !orbitals interpolated, gradient and Laplacian calculated
    if(l_need_phase) then
      call feval_multi_ubspline_3d_z_vgh(bspline_coefficient_pointer(1),&
                                         r(1), r(2), r(3),              &
                                         orb_tmp_complex,               &
                                         dorb_tmp_complex,              &
                                         hessian_bsplines_complex)
      dorb_tmp_complex=matmul(rlatt_inv_transpose,dorb_tmp_complex)
      !hessian has entries 1= d2f/dx2,  2= d2f/dxdy, 3=d2f/dxdz,
      !                    4= d2f/dydx, 5= d2f/dy2,  6=d2f/dydz,
      !                    7= d2f/dzdx, 8= d2f/dzdy, 9=d2f/dzdz.
      do i=1,nspline
         hessian_tmp_complex=reshape(hessian_bsplines_complex(:,i),(/ndim,ndim/))
         ddorb_tmp_complex(i)=sum(rlatt_inv_transpose_2*hessian_tmp_complex)
      enddo
      call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,r_cart,ng1d_sim, &
                        cos_k,sin_k,dcos_k,dsin_k,ddcos_k,ddsin_k,rkvec_shift)
      do ikvec=1,nkvec
        complex_phase=cmplx(cos_k(ikvec),sin_k(ikvec))
        dcomplex_phase(:)=cmplx(dcos_k(:,ikvec),dsin_k(:,ikvec))
        ddcomplex_phase=cmplx(ddcos_k(ikvec),ddsin_k(ikvec))
        do i=1,nband(ikvec)
           iorb=iorb+1
           ispline=(ikvec-1)*nkvec+i
           if(ireal_imag(iorb) == 1) then !use real part
             orb_bsplines(iorb)   =dble(orb_tmp_complex(ispline)*complex_phase)
             dorb_bsplines(:,iorb)=dble(orb_tmp_complex(ispline)*dcomplex_phase(:)+&
                                        dorb_tmp_complex(:,ispline)*complex_phase)
             ddorb_bsplines(iorb) =dble(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                 2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                      ddorb_tmp_complex(ispline)   *complex_phase)
           else if(ireal_imag(iorb) == 2) then!use imaginary part
              orb_bsplines(iorb)   =aimag(orb_tmp_complex(ispline)   *complex_phase)
              dorb_bsplines(:,iorb)=aimag(orb_tmp_complex(ispline)  *dcomplex_phase(:)+&
                                          dorb_tmp_complex(:,ispline)*complex_phase)
              ddorb_bsplines(iorb) =aimag(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                   2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                        ddorb_tmp_complex(ispline)   *complex_phase)
           else !use real and imaginary parts
              orb_bsplines(iorb)   =dble(orb_tmp_complex(ispline)   *complex_phase)
              dorb_bsplines(:,iorb)=dble(orb_tmp_complex(ispline)  *dcomplex_phase(:)+&
                                        dorb_tmp_complex(:,ispline) *complex_phase)
              ddorb_bsplines(iorb) =dble(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                  2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                       ddorb_tmp_complex(ispline)   *complex_phase)
              iorb=iorb+1
              orb_bsplines(iorb)   =aimag(orb_tmp_complex(ispline)   *complex_phase)
              dorb_bsplines(:,iorb)=aimag(orb_tmp_complex(ispline)  *dcomplex_phase(:)+&
                                         dorb_tmp_complex(:,ispline) *complex_phase)
              ddorb_bsplines(iorb) =aimag(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                   2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                        ddorb_tmp_complex(ispline)   *complex_phase)
           endif !real or imaginary part
        enddo !orbitals k-point
      enddo !k-vectors
    else !no phase (Gamma)
       call feval_multi_ubspline_3d_d_vgh(bspline_coefficient_pointer(1),&
                                          r(1),r(2),r(3),                &
                                          orb_tmp,dorb_tmp,              &
                                          hessian_bsplines)
       orb_bsplines=orb_tmp
       dorb_bsplines=matmul(rlatt_inv_transpose,dorb_tmp)
       do i=1,nspline
          hessian_tmp=reshape(hessian_bsplines(:,i),(/ndim,ndim/))
          ddorb_bsplines(i)=sum(rlatt_inv_transpose_2*hessian_tmp)
       enddo
    endif !need phase
  case (1) !orbitals and Laplacian interpolated, gradient calculated
    if(l_need_phase) then
      call feval_multi_ubspline_3d_z_vg(bspline_coefficient_pointer(1),&
                                        r(1), r(2), r(3),              &
                                        orb_tmp_complex,               &
                                        dorb_tmp_complex)
      dorb_tmp_complex=matmul(rlatt_inv_transpose,dorb_tmp_complex)
      call feval_multi_ubspline_3d_z(bspline_coefficient_pointer(2),&
                                       r(1), r(2), r(3),              &
                                       ddorb_tmp_complex)
      call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,r_cart,ng1d_sim, &
                        cos_k,sin_k,dcos_k,dsin_k,ddcos_k,ddsin_k,rkvec_shift)
      do ikvec=1,nkvec
        complex_phase=cmplx(cos_k(ikvec),sin_k(ikvec))
        dcomplex_phase(:)=cmplx(dcos_k(:,ikvec),dsin_k(:,ikvec))
        ddcomplex_phase=cmplx(ddcos_k(ikvec),ddsin_k(ikvec))
        do i=1,nband(ikvec)
           iorb=iorb+1
           ispline=(ikvec-1)*nkvec+i
           if(ireal_imag(iorb) == 1) then !use real part
             orb_bsplines(iorb)   =dble(orb_tmp_complex(ispline)*complex_phase)
             dorb_bsplines(:,iorb)=dble(orb_tmp_complex(ispline)*dcomplex_phase(:)+&
                                        dorb_tmp_complex(:,ispline)*complex_phase)
             ddorb_bsplines(iorb) =dble(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                 2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                      ddorb_tmp_complex(ispline)   *complex_phase)
           else if(ireal_imag(iorb) == 2) then!use imaginary part
              orb_bsplines(iorb)   =aimag(orb_tmp_complex(ispline)   *complex_phase)
              dorb_bsplines(:,iorb)=aimag(orb_tmp_complex(ispline)  *dcomplex_phase(:)+&
                                          dorb_tmp_complex(:,ispline)*complex_phase)
              ddorb_bsplines(iorb) =aimag(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                   2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                        ddorb_tmp_complex(ispline)   *complex_phase)
           else !use real and imaginary parts
              orb_bsplines(iorb)   =dble(orb_tmp_complex(ispline)   *complex_phase)
              dorb_bsplines(:,iorb)=dble(orb_tmp_complex(ispline)  *dcomplex_phase(:)+&
                                        dorb_tmp_complex(:,ispline) *complex_phase)
              ddorb_bsplines(iorb) =dble(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                  2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                       ddorb_tmp_complex(ispline)   *complex_phase)
              iorb=iorb+1
              orb_bsplines(iorb)   =aimag(orb_tmp_complex(ispline)   *complex_phase)
              dorb_bsplines(:,iorb)=aimag(orb_tmp_complex(ispline)  *dcomplex_phase(:)+&
                                         dorb_tmp_complex(:,ispline) *complex_phase)
              ddorb_bsplines(iorb) =aimag(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                   2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                        ddorb_tmp_complex(ispline)   *complex_phase)
           endif !real or imaginary part
        enddo !orbitals k-point
      enddo !k-vectors
    else !no phase (Gamma)
       call feval_multi_ubspline_3d_d_vg(bspline_coefficient_pointer(1), &
                                         r(1), r(2), r(3),               &
                                         orb_tmp, dorb_tmp)
       call feval_multi_ubspline_3d_d(bspline_coefficient_pointer(2), &
                                        r(1), r(2), r(3),               &
                                        ddorb_tmp)
       orb_bsplines=orb_tmp
       dorb_bsplines=matmul(rlatt_inv_transpose,dorb_tmp)
       ddorb_bsplines=ddorb_tmp
    endif !need phase
  case (2) !orbitals, gradient and Laplacian interpolated
    if(l_need_phase) then
      call feval_multi_ubspline_3d_z(bspline_coefficient_pointer(1),&
                                     r(1), r(2), r(3),              &
                                     orb_tmp_complex)
      call feval_multi_ubspline_3d_z(bspline_coefficient_pointer(2),&
                                     r(1), r(2), r(3),              &
                                     ddorb_tmp_complex)
      call feval_multi_ubspline_3d_z(bspline_coefficient_pointer(3),&
                                     r(1), r(2), r(3),              &
                                     dorb_tmp_complex(1,:))
      call feval_multi_ubspline_3d_z(bspline_coefficient_pointer(4),&
                                     r(1), r(2), r(3),              &
                                     dorb_tmp_complex(2,:))
      call feval_multi_ubspline_3d_z(bspline_coefficient_pointer(5),&
                                     r(1), r(2), r(3),              &
                                     dorb_tmp_complex(3,:))
      call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,r_cart,ng1d_sim, &
                        cos_k,sin_k,dcos_k,dsin_k,ddcos_k,ddsin_k,rkvec_shift)
      do ikvec=1,nkvec
        complex_phase=cmplx(cos_k(ikvec),sin_k(ikvec))
        dcomplex_phase(:)=cmplx(dcos_k(:,ikvec),dsin_k(:,ikvec))
        ddcomplex_phase=cmplx(ddcos_k(ikvec),ddsin_k(ikvec))
        do i=1,nband(ikvec)
           iorb=iorb+1
           ispline=(ikvec-1)*nkvec+i
           if(ireal_imag(iorb) == 1) then !use real part
             orb_bsplines(iorb)   =dble(orb_tmp_complex(ispline)*complex_phase)
             dorb_bsplines(:,iorb)=dble(orb_tmp_complex(ispline)*dcomplex_phase(:)+&
                                        dorb_tmp_complex(:,ispline)*complex_phase)
             ddorb_bsplines(iorb) =dble(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                 2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                      ddorb_tmp_complex(ispline)   *complex_phase)
           else if(ireal_imag(iorb) == 2) then!use imaginary part
              orb_bsplines(iorb)   =aimag(orb_tmp_complex(ispline)   *complex_phase)
              dorb_bsplines(:,iorb)=aimag(orb_tmp_complex(ispline)  *dcomplex_phase(:)+&
                                          dorb_tmp_complex(:,ispline)*complex_phase)
              ddorb_bsplines(iorb) =aimag(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                   2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                        ddorb_tmp_complex(ispline)   *complex_phase)
           else !use real and imaginary parts
              orb_bsplines(iorb)   =dble(orb_tmp_complex(ispline)   *complex_phase)
              dorb_bsplines(:,iorb)=dble(orb_tmp_complex(ispline)  *dcomplex_phase(:)+&
                                        dorb_tmp_complex(:,ispline) *complex_phase)
              ddorb_bsplines(iorb) =dble(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                  2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                       ddorb_tmp_complex(ispline)   *complex_phase)
              iorb=iorb+1
              orb_bsplines(iorb)   =aimag(orb_tmp_complex(ispline)   *complex_phase)
              dorb_bsplines(:,iorb)=aimag(orb_tmp_complex(ispline)  *dcomplex_phase(:)+&
                                         dorb_tmp_complex(:,ispline) *complex_phase)
              ddorb_bsplines(iorb) =aimag(orb_tmp_complex(ispline) *ddcomplex_phase    +&
                                   2*sum(dorb_tmp_complex(:,ispline)*dcomplex_phase(:))+&
                                        ddorb_tmp_complex(ispline)   *complex_phase)
           endif !real or imaginary part
        enddo !orbitals k-point
      enddo !k-vectors
    else !no phase (Gamma)
       call feval_multi_ubspline_3d_d(bspline_coefficient_pointer(1), &
                                      r(1), r(2), r(3),               &
                                      orb_tmp)
       call feval_multi_ubspline_3d_d(bspline_coefficient_pointer(2), &
                                      r(1), r(2), r(3),               &
                                      ddorb_tmp)
       call feval_multi_ubspline_3d_d(bspline_coefficient_pointer(3), &
                                      r(1), r(2), r(3),               &
                                      dorb_tmp(1,:))
       call feval_multi_ubspline_3d_d(bspline_coefficient_pointer(4), &
                                      r(1), r(2), r(3),               &
                                      dorb_tmp(2,:))
       call feval_multi_ubspline_3d_d(bspline_coefficient_pointer(5), &
                                      r(1), r(2), r(3),               &
                                      dorb_tmp(3,:))
       orb_bsplines=orb_tmp
       dorb_bsplines=dorb_tmp
       ddorb_bsplines=ddorb_tmp
    endif !need phase
  end select

  end subroutine evaluate_bsplines_with_derivatives

! ==============================================================================
  subroutine orbitals_pw_primitive(x,orb_pw,dorb_pw,ddorb_pw)
! ------------------------------------------------------------------------------
! Description : Evaluate the plane wave orbitals on the primitive cell alone
!
! Created     : W. Parker, 21 Oct 2008
! ------------------------------------------------------------------------------
  implicit none
  !need ndim,ngvec,glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,ng1d,
  !     rkvec_shift,k_inv,ireal_imag,nband
  include 'commons.h'

! input
  real(dp),intent(in)     :: x(ndim)

! output
  real(dp),intent(out)    ::   orb_pw(nkvec,norb_primitive,2)
  real(dp),intent(out)    ::  dorb_pw(:,:,:,:)
  real(dp),intent(out)    :: ddorb_pw(:,:,:)

! local
  integer                             :: iorb
  integer                             :: jorb
  integer                             :: ikvec
  integer                             :: iband
  integer                             :: ig
  real(dp)                            ::   cos_rp
  real(dp)                            ::   sin_rm
  real(dp)                            ::   cos_ip
  real(dp)                            ::   sin_im
  real(dp)                            ::  dcos_rp(ndim)
  real(dp)                            ::  dsin_rm(ndim)
  real(dp)                            ::  dcos_ip(ndim)
  real(dp)                            ::  dsin_im(ndim)
  real(dp)                            :: ddcos_rp
  real(dp)                            :: ddsin_rm
  real(dp)                            :: ddcos_ip
  real(dp)                            :: ddsin_im
  real(dp),allocatable                ::   cos_g(:)
  real(dp),allocatable                ::   sin_g(:)
  real(dp),allocatable                ::  dcos_g(:,:)
  real(dp),allocatable                ::  dsin_g(:,:)
  real(dp),allocatable                :: ddcos_g(:)
  real(dp),allocatable                :: ddsin_g(:)

! allocation
  call alloc (  'cos_g',   cos_g, ngvec)
  call alloc (  'sin_g',   sin_g, ngvec)
  call alloc ( 'dcos_g',  dcos_g, ndim, ngvec)
  call alloc ( 'dsin_g',  dsin_g, ndim, ngvec)
  call alloc ('ddcos_g', ddcos_g, ngvec)
  call alloc ('ddsin_g', ddsin_g, ngvec)

  !get exp(i G.r)
  call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,&
                    gvec,igvec,ngvec_orb,x,ng1d,  &
                    cos_g,sin_g,dcos_g,dsin_g,    &
                    ddcos_g,ddsin_g,rkvec_shift)

  iorb=0
  jorb=0

  select case (igrad_lap)
  case (0) !orbitals only

    do ikvec=1,nkvec
       do iband=1,nband(ikvec)
          jorb=jorb+1
          cos_rp=0
          sin_rm=0
          cos_ip=0
          sin_im=0
          !plane wave sum for single orbital
          do ig=2,ngvec_orb
             cos_rp=cos_rp+cos_g(ig)*c_rp(ig,jorb)
             sin_rm=sin_rm+sin_g(ig)*c_rm(ig,jorb)
             cos_ip=cos_ip+cos_g(ig)*c_ip(ig,jorb)
             sin_im=sin_im+sin_g(ig)*c_im(ig,jorb)
          enddo

          !add correct sine and cosine sums with real and imaginary parts
          !of first coefficient (for G={0,0,0})
          orb_pw(ikvec,iband,1)=c_rp(1,jorb)+cos_rp-sin_im
          orb_pw(ikvec,iband,2)=c_ip(1,jorb)+cos_ip+sin_rm

       enddo !bands
    enddo !k-vectors

  case (1) !orbitals and Laplacian

    do ikvec=1,nkvec
       do iband=1,nband(ikvec)
          jorb=jorb+1

          !zero the sum variables for each band
          cos_rp=0
          sin_rm=0
          cos_ip=0
          sin_im=0
          ddcos_rp=0
          ddsin_rm=0
          ddcos_ip=0
          ddsin_im=0

          !plane wave sum for single orbital
          do ig=2,ngvec_orb
             cos_rp=cos_rp+cos_g(ig)*c_rp(ig,jorb)
             sin_rm=sin_rm+sin_g(ig)*c_rm(ig,jorb)
             cos_ip=cos_ip+cos_g(ig)*c_ip(ig,jorb)
             sin_im=sin_im+sin_g(ig)*c_im(ig,jorb)
             ddcos_rp=ddcos_rp+ddcos_g(ig)*c_rp(ig,jorb)
             ddsin_rm=ddsin_rm+ddsin_g(ig)*c_rm(ig,jorb)
             ddcos_ip=ddcos_ip+ddcos_g(ig)*c_ip(ig,jorb)
             ddsin_im=ddsin_im+ddsin_g(ig)*c_im(ig,jorb)
          enddo

          !add correct sine and cosine sums with real and imaginary parts
          !of first coefficient (for G={0,0,0})
          orb_pw(ikvec,iband,1)  =c_rp(1,jorb)+cos_rp-sin_im
          orb_pw(ikvec,iband,2)  =c_ip(1,jorb)+cos_ip+sin_rm
          ddorb_pw(ikvec,iband,1)=ddcos_rp-ddsin_im
          ddorb_pw(ikvec,iband,2)=ddcos_ip+ddsin_rm

       enddo !bands
    enddo !k-vectors

  case (2) !orbitals, gradient and Laplacian

    do ikvec=1,nkvec
       do iband=1,nband(ikvec)
          jorb=jorb+1

          !zero the sum variables for each band
          cos_rp=0
          sin_rm=0
          cos_ip=0
          sin_im=0
          dcos_rp=0
          dsin_rm=0
          dcos_ip=0
          dsin_im=0
          ddcos_rp=0
          ddsin_rm=0
          ddcos_ip=0
          ddsin_im=0

          !plane wave sum for single orbital
          do ig=2,ngvec_orb
             cos_rp=cos_rp+cos_g(ig)*c_rp(ig,jorb)
             sin_rm=sin_rm+sin_g(ig)*c_rm(ig,jorb)
             cos_ip=cos_ip+cos_g(ig)*c_ip(ig,jorb)
             sin_im=sin_im+sin_g(ig)*c_im(ig,jorb)

             dcos_rp(:)=dcos_rp(:)+dcos_g(:,ig)*c_rp(ig,jorb)
             dsin_rm(:)=dsin_rm(:)+dsin_g(:,ig)*c_rm(ig,jorb)
             dcos_ip(:)=dcos_ip(:)+dcos_g(:,ig)*c_ip(ig,jorb)
             dsin_im(:)=dsin_im(:)+dsin_g(:,ig)*c_im(ig,jorb)

             ddcos_rp=ddcos_rp+ddcos_g(ig)*c_rp(ig,jorb)
             ddsin_rm=ddsin_rm+ddsin_g(ig)*c_rm(ig,jorb)
             ddcos_ip=ddcos_ip+ddcos_g(ig)*c_ip(ig,jorb)
             ddsin_im=ddsin_im+ddsin_g(ig)*c_im(ig,jorb)
          enddo

          !add correct sine and cosine sums with real and imaginary parts
          !of first coefficient (for G={0,0,0})
          orb_pw(ikvec,iband,1)   =c_rp(1,jorb)+cos_rp   -  sin_im
          orb_pw(ikvec,iband,2)   =c_ip(1,jorb)+cos_ip   +  sin_rm
          dorb_pw(ikvec,:,iband,1)=            dcos_rp(:)- dsin_im(:)
          dorb_pw(ikvec,:,iband,2)=            dcos_ip(:)+ dsin_rm(:)
          ddorb_pw(ikvec,iband,1) =           ddcos_rp   -ddsin_im
          ddorb_pw(ikvec,iband,2) =           ddcos_ip   +ddsin_rm

       enddo !bands
    enddo !k-vectors
  end select

! deallocation
  call release (  'cos_g',   cos_g)
  call release (  'sin_g',   sin_g)
  call release ( 'dcos_g',  dcos_g)
  call release ( 'dsin_g',  dsin_g)
  call release ('ddcos_g', ddcos_g)
  call release ('ddsin_g', ddsin_g)


  end subroutine orbitals_pw_primitive

end module bsplines_mod

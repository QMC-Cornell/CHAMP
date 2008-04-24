Module periodic_jastrow_mod

  use all_tools_mod
  use crystal_symmetries_mod
  use eloc_mod
  use pjasen_mod
  use pjasee_mod

  implicit none

  integer                                :: i_pjas_en_read = 0
  integer                                :: i_pjas_ee_read = 0
  real(dp), allocatable                  :: gn_pjas  (:)

# if defined (PATHSCALE)
  real(dp)                               :: pjas_ee_read(max_double_array_len)
  real(dp)                               :: pjas_en_read(max_double_array_len)
# else
  real(dp), allocatable                  :: pjas_ee_read(:)
  real(dp), allocatable                  :: pjas_en_read(:)
# endif


contains


  !===========================================================================
  subroutine  periodic_jastrow_menu
!---------------------------------------------------------------------------
! Description : menu for periodic jastrow
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    include 'common_pjase.h'
    include 'commons_jaspar6.h'
    integer                              :: icutjasc
    real(dp)                             :: cutjasc
    common /jas_c_cut/ cutjasc, icutjasc


    icutjasc = 0
    cutjasc= 0.0_dp

    ! loop over menu lines
    do

       call get_next_word (word)

       select case(trim(word))

       case ('help')
          write(6,'(2a)') 'help for periodic_jastrow menu:'
          write(6,'(2a)') 'do_pjasen = true to include e-n  periodic jastrow term'
          write(6,'(2a)') 'do_pjasee = true to include e-e  periodic jastrow term'
          write(6,'(2a)') 'ecut_en:  cuttoff to generate the stars in primitive cell'
          write(6,'(2a)') 'ecut_ee:  cuttoff to generate the stars in simulation cell'
          write(6,'(2a)') 'nstar_en:  actual number of stars included in e-n periodic jastrow term'
          write(6,'(2a)') 'nstar_ee:  actual number of stars included in e-e periodic jastrow term'
          write(6,'(2a)') 'pjas_en_read (list): initial starting point fo e-n periodic jastrow term'
          write(6,'(2a)') 'pjas_ee_read (list): initial starting point fo e-e periodic jastrow term'
          write(6,'(2a)') 'sym_file (list): Names of two existing filenames with symmetries of the primitive and simulation cells'
          write(6,'(2a)') 'end'

       case ('icutjasc')
          call get_next_value (icutjasc)

       case ('cutjasc')
          call get_next_value (cutjasc)

!!! commented for now
!!$       case ('sym_file')
!!$          call get_next_value_list ("symmetry_file", symmetry_file, iread_sym_file )
!!$          read_sym_file = 1
!!$          if (iread_sym_file .ne. 2) stop "two symmetry input files needed"

       case ('nstar_en')
          call get_next_value (nstar_en)

       case ('nstar_ee')
          call get_next_value (nstar_ee)

       case ('inversion')
          call get_next_value (inversion)
          if (.not. inversion) n_inv = 2

       case ('shift_cos')
          call get_next_value (shift_cos)
          shift_cos = 4* atan (1.0) * shift_cos
          writE(6,*) "shift_cos = ", shift_cos

       case ('do_pjasen')
          call get_next_value (do_pjasen)

       case ('do_pjasee')
          call get_next_value (do_pjasee)

       case ('ecut_en')
          call get_next_value (ecut_en)

       case ('ecut_ee')
          call get_next_value (ecut_ee)

       case ('pjas_en_read')
          call get_next_value_list ('pjas_en_read', pjas_en_read, i_pjas_en_read)

       case ('pjas_ee_read')
          call get_next_value_list ('pjas_ee_read', pjas_ee_read, i_pjas_ee_read)

       case ('end')
          exit

       case default
          write(6,'(3a)') ': unknown keyword = ',trim(word)
          call die ('periodic jastrow menu')
       end select

    enddo ! end loop over menu lines

    if (cutjasc .eq. 0) cutjasc= cutjas_en



    ido_pjasen = 0
    if (do_pjasen) then
       write(6,*) "e-n periodic jastrow is included "
       ido_pjasen = 1
    else
       nstar_en = 0
    endif


    ido_pjasee = 0
    if (do_pjasee) then
       write(6,*) "e-e periodic jastrow is included "
       ido_pjasee = 1
    else
       nstar_ee = 0
    endif

    ido_pjas= 0  !! common block
    if (do_pjasen .or. do_pjasee) then
       do_pjas= .true.
       ido_pjas= 1
    endif

!!!!!!!!!! although not used they have to be defined

    call object_modified_by_index (cos_star_ee_index)

    call object_modified_by_index (cos_star_en_index)

    call object_modified_by_index (star_en_index)

    call object_modified_by_index (grad_star_en_index)

    call object_modified_by_index (star_ee_index)

    call object_modified_by_index (grad_star_ee_index)


    if (ido_pjas == 0 ) return


    call planewaves_interface
    write(6,*) "done planewave_interface "

    write(6,*) "do_pjasen  = ", do_pjasen

    if (nstar_en > 0) then
       if (nstar_en > nstar) nstar_en = nstar
       nstar= nstar_en
       write(6,*) "Number of e-n stars included = ", nstar
    endif


    if (nstar_ee > 0) then
       if (nstar_ee > nstar) nstar_ee = nstar_sim
       write(6,*) "Number of e-e stars included = ", nstar_ee
    endif
    call object_provide ('param_pjas_nb')


  end subroutine periodic_jastrow_menu



  subroutine pjas_init_bld
!---------------------------------------------------------------------------
! Description : initial subroutine to define arrays and dimensions
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    integer                              :: ipar, i, j
    real(dp)                             :: sum

    if (header_exe) then
       call object_create ('param_pjas_nb', param_pjas_nb_index)
       call object_create ('param_pjasen_nb', param_pjasen_nb_index)
       call object_create ('param_pjasee_nb', param_pjasee_nb_index)
       call object_needed ('ndim')
       call object_needed ('nelec')
       return
    endif

    if (l_opt_pjasen .and. .not. do_pjasen) then
       stop "do_pjasen and pjasen inconsistent keywords"
    endif

    if (l_opt_pjasee .and. .not. do_pjasee) then
       stop "do_pjasee and pjasee inconsistent keywords"
    endif

    param_pjasen_nb = nstar_en * n_inv

    param_pjasee_nb = nstar_ee

    param_pjas_nb = param_pjasee_nb + param_pjasen_nb

    call object_alloc ("pjas_parms", pjas_parms, param_pjas_nb, mwf )
    pjas_parms = 0

    if (i_pjas_en_read >0) then
       if (i_pjas_en_read > param_pjasen_nb) i_pjas_en_read = param_pjasen_nb
       do i=1, i_pjas_en_read
          pjas_parms(i,1) = pjas_en_read(i)
       enddo
       write(6,'(a,100F10.6)') "Reading pjas_parms_en", pjas_en_read(1:i_pjas_en_read)
    endif



    if (i_pjas_ee_read >0) then
       if (i_pjas_ee_read > param_pjasee_nb) i_pjas_ee_read = param_pjasee_nb
       do i=1, i_pjas_ee_read
          pjas_parms(i+ param_pjasen_nb,1) = pjas_ee_read(i)
       enddo
       write(6,'(a,100F10.6)') "Reading pjas_parms_ee", pjas_ee_read(1:i_pjas_ee_read)
    endif


!!! !! this is introduced to fix compiling mpi version
    ndim_pj= ndim


!!! allocations of global arrays

    call object_alloc ("cos_star_ee ", cos_star_ee , nelec * (nelec-1)/2,  nstar_ee)
    call object_alloc ("cos_star_en ", cos_star_en , nelec,  nstar_en)
    call object_modified_by_index (cos_star_ee_index)
    call object_modified_by_index (cos_star_en_index)

    if (param_pjasee_nb >0) then
       call object_alloc ('dpsi_pjasee', dpsi_pjasee, param_pjasee_nb)
       if (nopt_iter > 0) then
          call object_alloc ('lap_dpsi_pjasee', lap_dpsi_pjasee, nelec, param_pjasee_nb)
          call object_alloc ('grad_dpsi_pjasee', grad_dpsi_pjasee, ndim_pj, nelec , param_pjasee_nb)
       endif
    endif

    if (param_pjasen_nb >0) then
       call object_alloc ('dpsi_pjasen', dpsi_pjasen, param_pjasen_nb)
       if (nopt_iter > 0) then
          call object_alloc ('lap_dpsi_pjasen', lap_dpsi_pjasen, nelec, param_pjasen_nb)
          call object_alloc ('grad_dpsi_pjasen', grad_dpsi_pjasen, ndim_pj, nelec , param_pjasen_nb)
       endif
    endif

    call object_modified_by_index (dpsi_pjasee_index)
    call object_modified_by_index (grad_dpsi_pjasee_index)
    call object_modified_by_index (lap_dpsi_pjasee_index)

    call object_modified_by_index (dpsi_pjasen_index)
    call object_modified_by_index (grad_dpsi_pjasen_index)
    call object_modified_by_index (lap_dpsi_pjasen_index)


    if (index(mode,'fit') /= 0) then
       write(6,*) "Long range periodic jastrow is not implemented with fit"
       stop "Long range periodic jastrow is not implemented with fit"
    endif

    !! allocate global variables used in pjas_deriv_jas and pjas_jas
    call object_alloc ( "pjasfso", pjasfso, nelec, nelec)
    call object_alloc ( "pjasfijo", pjasfijo, ndim, nelec, nelec)
    call object_alloc ( "pjasd2ijo", pjasd2ijo,  nelec, nelec)
    call object_alloc ( "pjasfjo", pjasfjo, ndim, nelec)

    !! allocate global variables used in pjas_deriv_jas
    if (nopt_iter > 0) then
       call object_alloc ( "pengo", pengo, nelec, nelec, param_pjasen_nb)
       call object_alloc ( "peego", peego, nelec, nelec, param_pjasee_nb)
    endif
  end subroutine pjas_init_bld


  subroutine deloc_pot_nloc_pjas_bld
!---------------------------------------------------------------------------
! Description : derivatives of the potential part of the local energy with respect to the
!               periodic Jastrow parameters
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'

    ! local
    integer ist

    ! header
    if (header_exe) then
       call object_create ('deloc_pot_nloc_pjas',  deloc_pot_nloc_pjas_index)
       call object_needed ('dvpsp_pjas')
       call object_needed ('dpsi_pjas')
       call object_needed ('psid_pjas')
       call object_needed ('eloc_pot_nloc')
       return
    endif

    ! allocations
    call object_alloc ('deloc_pot_nloc_pjas',deloc_pot_nloc_pjas, param_pjas_nb)

    do ist = 1 , param_pjas_nb
       deloc_pot_nloc_pjas (ist)= dvpsp_pjas (ist)/psid_pjas - eloc_pot_nloc * dpsi_pjas (ist)
    enddo

  end subroutine deloc_pot_nloc_pjas_bld


  subroutine dpsi_pjas_bld
!---------------------------------------------------------------------------
! Description : Logarithmic derivatives of the wavefunction with respect to the peridoic Jastrow parameters
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    integer                              :: ipar, i, j
    real(dp)                             :: sum

    if (header_exe) then
       call object_create ('dpsi_pjas', dpsi_pjas_index)
       call object_needed ('dpsi_pjasen')
       call object_needed ('dpsi_pjasee')
       call object_needed ('ndim')
       call object_needed ('nelec')
       return
    endif

    call object_alloc ('dpsi_pjas', dpsi_pjas, param_pjas_nb)

    if (l_opt_pjasen) then
       call object_provide ('dpsi_pjasen')
       dpsi_pjas (1: param_pjasen_nb)= dpsi_pjasen
    endif

    if (l_opt_pjasee) then
       call object_provide ('dpsi_pjasee')
       dpsi_pjas (param_pjasen_nb+ 1: param_pjas_nb)= dpsi_pjasee
    endif

  end subroutine dpsi_pjas_bld



  subroutine deloc_pjas_bld
!---------------------------------------------------------------------------
! Description : derivatives of the local energy with respect to the periodic Jastrow parameters.
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    integer                              :: ieta
    include 'commons.h'
    ! header
    if (header_exe) then
       call object_create ('deloc_pjas',deloc_pjas_index)
       call object_needed ('deloc_pot_nloc_pjas')
       call object_needed ('nelec')
       call object_needed ('ndim')
       call object_needed ('param_pjasen_nb')
       call object_needed ('param_pjasee_nb')
       return
    endif

    call object_alloc ('deloc_pjas', deloc_pjas, param_pjas_nb)

    if (param_pjasen_nb >0) then
       call object_provide_by_index (deloc_pjasen_index)
       do ieta=1, param_pjasen_nb
          deloc_pjas (ieta)=  deloc_pjasen (ieta)
       enddo
    endif

    if (param_pjasee_nb >0) then
       call object_provide_by_index (deloc_pjasee_index)
       do ieta = 1, param_pjasee_nb
          deloc_pjas (ieta + param_pjasen_nb)=  deloc_pjasee (ieta)
       enddo
    endif

    if (nloc >0 ) then
       call object_provide_by_index (deloc_pot_nloc_pjas_index)
       deloc_pjas =  deloc_pjas + deloc_pot_nloc_pjas
    endif


  end subroutine deloc_pjas_bld


  subroutine copy_pjas
!---------------------------------------------------------------------------
! Description : copy pjas_parms  used in opt. step
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    integer                              :: i, ifr
    do ifr =2 , nforce
       do i=1, param_pjas_nb !param_pjasen_nb
          pjas_parms (i,ifr) =  pjas_parms (i,1)
       end do
    enddo
  end subroutine copy_pjas



  subroutine save_pjas
!---------------------------------------------------------------------------
! Description : save pjas parms  used in opt. step
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    integer i
    call object_alloc ('pjas_parms_sav', pjas_parms_sav, param_pjas_nb)
    do i=1, param_pjas_nb
       pjas_parms_sav(i)=  pjas_parms (i,1)
    enddo
  end subroutine save_pjas


  subroutine restore_pjas
!---------------------------------------------------------------------------
! Description : restore pjas parms used in opt. step
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    integer i

    do i=1, param_pjas_nb
       pjas_parms(i,1)=  pjas_parms_sav (i)
    enddo
  end subroutine restore_pjas


  subroutine deriv_nonloc_pjas ( iel, rvec, value) !, gn_pjas)
!---------------------------------------------------------------------------
! Description : drivers for  deriv_nonloc_pjas_en and     deriv_nonloc_pjas_ee
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    include "common_jasn.h"

    integer                              :: iel, ist
    real(dp)                             :: rvec (3,melec)
    real(dp)                             :: xvec (3), value1, temp , value2, value

    call object_alloc ("gn_pjas", gn_pjas, param_pjas_nb)

    if (do_pjasen) then
       xvec = rvec(:,iel)
       call deriv_nonloc_pjas_en ( iel, xvec, value1)
       value = value + value1
       gn_pjas (1:param_pjasen_nb) =  gn_pjasen (1:param_pjasen_nb)
    endif

    if (do_pjasee) then
       call deriv_nonloc_pjas_ee ( iel, rvec, value2)
       value = value + value2
       gn_pjas (param_pjasen_nb+1:param_pjas_nb) =  gn_pjasee (1:param_pjasee_nb)
    endif

!!! update  fsumn
    fsumn = fsumn + value

  end subroutine deriv_nonloc_pjas


  subroutine nonloc_pjas ( iel, rvec, value)
!---------------------------------------------------------------------------
! Description : driver for nonloc_pjas_en  and nonloc_pjas_ee
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include "commons.h"
    include "common_jasn.h"
    real(dp)                             :: rvec (3,melec)
    integer                              :: iel, ist
    real(dp)                             :: xvec (3), value1, temp, value2, value

    if (do_pjasen) then
       xvec = rvec(:,iel)
       call nonloc_pjas_en ( iel, xvec, value1)
       value = value + value1
     endif

     if (do_pjasee) then
        call nonloc_pjas_ee ( iel, rvec, value2)
        value = value + value2
     endif

     !! update global fsumn
     !! subtract was done with fsuno
     fsumn = fsumn + value

  end subroutine nonloc_pjas


  subroutine  pjas_jas (xvec, rvec, v, d2, div_vj, fsum)
!---------------------------------------------------------------------------
! Description : drivers for pjasen_jas and pjasee_jas
!!! input are the contributions from the short range jastrow
!!! add contributions from long range jastrow
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    real (dp)                            :: rvec (3, melec*(melec-1)/2)
    real (dp)                            :: xvec (3, melec)
    real (dp)                            :: v (3,melec)
    real (dp)                            :: fsum, d2
    real (dp)                            :: div_vj (melec)
    integer                              :: i
    real(dp)                             :: pjasd2o_l, pjasfsumo_l

!!! common block variables
    pjasfijo = 0
    pjasfso = 0
    pjasd2ijo = 0
    pjasfjo = 0

!!! local variables
    pjasd2o_l = 0
    pjasfsumo_l = 0

    if (do_pjasen) then

       call pjasen_jas (xvec, v, pjasd2o_l, div_vj, pjasfsumo_l)

    endif


    if (do_pjasee) then

       call pjasee_jas (rvec, v, pjasd2o_l, div_vj, pjasfsumo_l)

    endif

    pjasd2o = pjasd2o_l
    pjasfsumo = pjasfsumo_l
    fsum = fsum + pjasfsumo_l
    d2  = d2 + pjasd2o_l

  end subroutine pjas_jas



  subroutine  pjas_deriv_jas (xvec,rvec, v, d2, div_vj, fsum)
!---------------------------------------------------------------------------
! Description : driver for pjas_deriv_jas_en and pjas_deriv_jas_ee
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    real (dp)                            :: rvec (3, melec*(melec-1)/2)
    real (dp)                            :: xvec (3, melec)
    real (dp)                            :: v (3,melec), div_vj (melec)
    real (dp)                            :: fsum, d2
    integer                              :: i
    real(dp)                             :: pjasd2o_l, pjasfsumo_l


    pengo = 0
    peego = 0
    pjasfijo = 0
    pjasfso = 0
    pjasd2ijo = 0
    pjasfjo = 0

!!! local variables
    pjasd2o_l = 0
    pjasfsumo_l = 0

    if (do_pjasen) then

       call pjasen_deriv_jas (xvec, v, pjasd2o_l,  div_vj, pjasfsumo_l )

    endif

    if (do_pjasee) then

       call pjasee_deriv_jas (rvec, v, pjasd2o_l,  div_vj, pjasfsumo_l )

    endif

!!! add to input
    pjasd2o = pjasd2o_l
    pjasfsumo = pjasfsumo_l

    fsum = fsum + pjasfsumo_l
    d2  = d2 + pjasd2o_l

  end subroutine pjas_deriv_jas



  subroutine pjas_jas_e (iel, xvec, rvec,  v, value)
!---------------------------------------------------------------------------
! Description : same as  pjas_deriv_jas when one electron move
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    include 'commons.h'
    include 'common_jasn.h'
    integer                              :: iel
    real (dp)                            :: rvec (3, melec*(melec-1)/2)
    real (dp)                            :: xvec (3, melec)
    real (dp)                            :: v (3,melec)
    real (dp)                            :: fsum, value
    integer                              :: i


    fsum = 0

    if (do_pjasen) then

       call pjasen_jas_e (iel, xvec, fsum)

    endif


    if (do_pjasee) then

       call pjasee_jas_e (iel, rvec, fsum)

    endif

    do i=1,nelec
       v (:, i) = fjn (1:3,i)
    enddo
    fsumn = fsumn + fsum
    value = value + fsum
  end subroutine pjas_jas_e



end Module periodic_Jastrow_mod

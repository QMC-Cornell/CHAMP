Module crystal_symmetries_mod

  use all_tools_mod
  use pjas_setup_mod
  implicit none


  !! symmetry file to be read
  !! format as in output of pw code
  integer                                :: iread_sym_file
  integer                                :: read_sym_file= 0
  character(len=max_string_len), allocatable  :: symmetry_file (:)

  real(dp)                               :: ecut_ee, ecut_en


  logical                                :: do_pjasee= .false.
  logical                                :: do_pjasen= .false.
  logical                                :: do_pjas= .false.
  integer                                :: maxkv, maxkv_sim


contains




  subroutine planewaves_interface
!---------------------------------------------------------------------------
! Description : interface for reading symmetry operations and defining the stars
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit  none
    include 'commons.h'

    logical                              :: iprt
    integer                              :: i , j, k
    integer                              :: prim_sim
    character                            :: isname (48) * 45
    integer, allocatable                 :: isymrel (:,:,:)

    write(6,*)
    write(6,*) "========================================================"
    write(6,*) "================ Periodic Jastrow setup ================"
    write(6,*) "========================================================"

!!$    write(6,*) "champ gvecs "
!!$    write(6,*) "Gvecs are ", ngvec
!!$    do i=1, ngvec
!!$       write(6,'(3I6,6F10.6)') igvec(:,i), gnorm(i), gvec(:,i)
!!$    enddo
!!$
!!$    write(6,*) "nkvec are ", nkvec
!!$    write(6,*) "rkvec are ", rkvec(:,1:nkvec)
!!$    write(6,*) "rkvec shift", rkvec_shift(:)
!!$    write(6,*) "========================= champ gvecs ============  "
!!$


    if (do_pjasen) then
       write(6,*)
       write(6,*) "Setup for primitive cell symmetries...."
       write(6,*)
       prim_sim=1 !! primitive cell

       call generate_read_symmetry (prim_sim,nsym_crys)
       call create_pw_lattice(prim_sim,nstar,nbasis_pw,nsym_crys,ecut_en, &
            &  rlatt,glatt,vcell,symrel,tnons)

       maxkv = int (maxval (real (kvec_pjasen)))

    endif


    if (do_pjasee) then
       write(6,*)
       write(6,*) "Setup for simulations  cell symmetries...."
       write(6,*)
       prim_sim=2 !! simulation cell

       call generate_read_symmetry (prim_sim,nsym_crys_sim)
       call create_pw_lattice(prim_sim,nstar_sim,nbasis_pw_sim,nsym_crys_sim,ecut_ee, &
            &  rlatt_sim,glatt_sim,vcell_sim,symrel_sim,tnons_sim)

       maxkv_sim = int (maxval (real (kvec_pjasee)))

    endif

    write(6,*)
    write(6,*) "========================================================"
    write(6,*) "============== End Periodic Jastrow setup =============="
    write(6,*) "========================================================"
    write(6,*)

  end subroutine planewaves_interface



  subroutine check_rot_symmetries (nsym,sym)
!---------------------------------------------------------------------------
! Description : A minimal check that the symmetry opertions are viable
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    integer,intent(in)                   :: nsym
    !arrays
    integer,intent(in)                   :: sym(3,3,nsym)
    integer                              :: determinant(nsym)

    call symdet(determinant,nsym,sym)

    call chkgrp(nsym, sym)

  end subroutine check_rot_symmetries


  subroutine symdet(determinant,nsym,sym)
!---------------------------------------------------------------------------
! Description : checks that the determinant of the symmetry operation is +-1                                                            !
!               adpated from abinit
!
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    !Arguments ------------------------------------
    !scalars
    integer,intent(in)                   :: nsym
    !arrays
    integer,intent(in)                   :: sym(3,3,nsym)
    integer,intent(out)                  :: determinant(nsym)

    !Local variables-------------------------------
    !scalars
    integer                              :: det,isym
    character(len=500)                   :: message

    do isym=1,nsym
       det=sym(1,1,isym)*sym(2,2,isym)*sym(3,3,isym)+&
            &     sym(2,1,isym)*sym(3,2,isym)*sym(1,3,isym)+&
            &     sym(1,2,isym)*sym(2,3,isym)*sym(3,1,isym) - &
            &    (sym(3,1,isym)*sym(2,2,isym)*sym(1,3,isym)+&
            &     sym(2,1,isym)*sym(1,2,isym)*sym(3,3,isym)+&
            &     sym(3,2,isym)*sym(2,3,isym)*sym(1,1,isym))
       if (abs(det)/=1) then
          write(* , * )  ' symdet: ERROR -',&
               & '  Abs(determinant) for symmetry number',isym,&
               &                                ' is',det
          stop "error in symdet"
       end if
       determinant(isym)=det
    end do
  end subroutine symdet


  subroutine chkgrp(nsym, symrel)
!---------------------------------------------------------------------------
! Description : checks that symmetry operations belong to a group by checking closure property (any pairs give an element of the group)
!                   adpated from abinit  WAS
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none

    !Arguments ------------------------------------
    !scalars
    integer,intent(in)                   :: nsym
    !arrays
    !  integer,intent(in)                :: symafm(nsym),symrel(3,3,nsym)
    integer,intent(in)                   :: symrel(3,3,nsym)

    !Local variables-------------------------------
    !scalars
    integer                              :: ii,isym,jj,jsym,kk,ksym,symafmchk,testeq
    character(len=500)                   :: message
    !arrays
    integer                              :: chk(3,3)

    do isym=1,nsym
       do jsym=1,nsym

          !  Compute the product of the two symmetries
          do ii=1,3
             do jj=1,3
                chk(ii,jj)=0
                do kk=1,3
                   chk(ii,jj)=chk(ii,jj)+&
                        &                 symrel(ii,kk,jsym)*symrel(kk,jj,isym)
                end do
             end do
          end do
          !   symafmchk=symafm(jsym)*symafm(isym)
          symafmchk =1

          !  Check that product array is one of original symmetries
          do ksym=1,nsym
             testeq=1
             do ii=1,3
                do jj=1,3
                   if(chk(ii,jj)/=symrel(ii,jj,ksym))testeq=0
                end do
             end do
             !    if(symafmchk/=symafm(ksym))testeq=0
             !   The test is positive
             if (testeq==1) exit
          end do

          !  The test is positive
          if(testeq==1)exit

          write(* , * )  '  Error: product of symmetries',isym,jsym,' is not in group.'

          ! End loop on jsym. Note that an "exit" instruction is present inside the loop
       end do

       !End loop on isym
    end do

  end subroutine chkgrp
  !!***



  subroutine generate_read_symmetry (prim_sim,nsym_crys)
!---------------------------------------------------------------------------
! Description : read the symmetry operations now (generation not implemented)
!
! Created     : W. A. Al-Saidi, June 2007
!---------------------------------------------------------------------------
    implicit none
    integer                              :: prim_sim, nsym_crys
    integer                              :: isym, ii
    integer                              :: tmp(9)
    logical                              :: exist_file
    character*9                          :: file1,file2,file3

!! always read the files

    !! use fixed names for these files for now because of a problem on Jaguar
    read_sym_file = 1
    file1 = "pri_sym"  !!for primitive lattice
    file2 = "sim_sym"  !!for simulation lattice
	

    if (read_sym_file==1) then

       if (prim_sim == 1) then
          file3 = file1		
       else
          file3 = file2
       endif

!!$       write(6,*) "symmetry operations file = ", trim (symmetry_file (prim_sim))
!!$       inquire (file=symmetry_file (prim_sim),exist=exist_file)
!!$       if (.not. exist_file) then
!!$          write(6,*) "symmetry file does not exist: file name :", symmetry_file (prim_sim)
!!$          call die ("generate_read_symmetry")
!!$       endif
!!$       open(18, file=symmetry_file(prim_sim))

       write(6,*) "symmetry operations file = ", file3
       inquire (file=file3,exist=exist_file)
       if (.not. exist_file) then
          write(6,*) "symmetry file does not exist: file name :", file3
          call die ("generate_read_symmetry")
       endif
       open(18, file=file3)



       read(18,*) nsym_crys

       if (prim_sim==1) then
          call object_alloc ("symrel", symrel, 3, 3, nsym_crys)
          call object_alloc ("tnons", tnons, 3,  nsym_crys)
       else
          call object_alloc ("symrel", symrel_sim, 3, 3, nsym_crys)
          call object_alloc ("tnons", tnons_sim, 3,  nsym_crys)
       endif


       write(6,*) "Number of symmetry operataions= ", nsym_crys

       do isym = 1, nsym_crys

          !          write(6,*) "reading ", isym

          if (prim_sim==1) then
             read(18,*) ii, tmp(1:9), tnons (:,isym)
             symrel(1,1,isym) = tmp(1)
             symrel(1,2,isym) = tmp(2)
             symrel(1,3,isym) = tmp(3)

             symrel(2,1,isym) = tmp(4)
             symrel(2,2,isym) = tmp(5)
             symrel(2,3,isym) = tmp(6)

             symrel(3,1,isym) = tmp(7)
             symrel(3,2,isym) = tmp(8)
             symrel(3,3,isym) = tmp(9)
          else
             read(18,*) ii, tmp(1:9), tnons_sim (:,isym)
             symrel_sim(1,1,isym) = tmp(1)
             symrel_sim(1,2,isym) = tmp(2)
             symrel_sim(1,3,isym) = tmp(3)

             symrel_sim(2,1,isym) = tmp(4)
             symrel_sim(2,2,isym) = tmp(5)
             symrel_sim(2,3,isym) = tmp(6)

             symrel_sim(3,1,isym) = tmp(7)
             symrel_sim(3,2,isym) = tmp(8)
             symrel_sim(3,3,isym) = tmp(9)
          endif

       enddo

       close (18)
    endif


!!$    write(6,'(a)')' symrel :'
!!$    do ii=1,nsym_crys
!!$       write(6,'(9i4)') int(symrel(:,:,ii))
!!$    end do
!!$
!!$    write(6,'(a)')' tnons :'
!!$    do ii=1,nsym_crys
!!$       write(6,'(9F10.6)') tnons(:,ii)
!!$    end do

    if (prim_sim==1) then
       call check_rot_symmetries (nsym_crys,int(symrel(:,:,:)))
    else
       call check_rot_symmetries (nsym_crys,int(symrel_sim(:,:,:)))
    endif

  end subroutine generate_read_symmetry


end Module crystal_symmetries_mod

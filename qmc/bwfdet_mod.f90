MODULE bwfdet_mod
!-----------------------------------------------------------------!
! BWFDET                                                          !
! ------                                                          !
! Setup and evaluation of orbitals, gradients and Laplacians for  !
! periodic system expanded in blip basis set.                     !
!                                                                 !
! DA 3-2001                                                       !
! modified for CHAMP by William Parker Nov. 2005 - 2006           !
! Modified: J. Toulouse, 22 Mar 2007: in order to avoid conflicts !
! with other commons, the following variable have been renamed:   !
! nkvec -> nkvec_bwfdet                                           !
! nband -> nband_bwfdet                                           !
! ndet -> ndet_bwfdet                                             !
! kvec -> kvec_bwfdet                                             !
!-----------------------------------------------------------------!

 USE constants_mod
 USE mpi_mod
 USE orbital_grid_mod, only: igrad_lap
 IMPLICIT NONE
 INTEGER,DIMENSION(3) :: nrbwf ! Blip grid
 INTEGER :: nwvec,nkvec_bwfdet,maxband,nemax,ndet_bwfdet
 INTEGER,PARAMETER :: maxgidx=50,nspin=2
 INTEGER,PARAMETER :: r2s_length=80
 INTEGER,DIMENSION(:),ALLOCATABLE :: nele
 INTEGER,DIMENSION(:,:),ALLOCATABLE :: gmap,nband_bwfdet,boccband
 INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: iprom_repl_idx,iadd_idx
 INTEGER,DIMENSION(:,:,:,:),ALLOCATABLE :: isub_idx
 REAL(dp) norm,painv(3,3)
 REAL(dp),PARAMETER :: one_over_twopi=1.d0/(3.14159265358979324d0*2.d0)
!numbers less than tolerance considered zero
 REAL(dp),PARAMETER :: tolerance=1.d-5
 REAL(dp),DIMENSION(:),ALLOCATABLE :: kdotr
 REAL(dp),DIMENSION(:,:),ALLOCATABLE :: kvec_bwfdet,kvec2
 REAL(dp),DIMENSION(:,:,:),ALLOCATABLE :: eigenvalue
! Blip coefficients
 REAL(dp),ALLOCATABLE :: avc(:,:,:,:), avc2(:,:,:,:)
 REAL(dp),ALLOCATABLE :: avclap(:,:,:,:), avclap2(:,:,:,:)
 REAL(dp),ALLOCATABLE :: avcgrad1(:,:,:,:), avcgrad12(:,:,:,:)
 REAL(dp),ALLOCATABLE :: avcgrad2(:,:,:,:), avcgrad22(:,:,:,:)
 REAL(dp),ALLOCATABLE :: avcgrad3(:,:,:,:), avcgrad32(:,:,:,:)
 COMPLEX(dp),PARAMETER :: zi=(0.d0,1.d0)
 COMPLEX(dp),ALLOCATABLE :: cavc(:,:,:,:,:), cavc2(:,:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: cavclap(:,:,:,:,:), cavclap2(:,:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: cavcgrad1(:,:,:,:,:), cavcgrad12(:,:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: cavcgrad2(:,:,:,:,:), cavcgrad22(:,:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: cavcgrad3(:,:,:,:,:), cavcgrad32(:,:,:,:,:)
 COMPLEX(dp),DIMENSION(:),ALLOCATABLE :: zdum,lzdum,ztemp
 COMPLEX(dp),DIMENSION(:,:),ALLOCATABLE :: gzdum
 LOGICAL spin_polarized,pwreal,open_unit(99)
 LOGICAL,DIMENSION(:),ALLOCATABLE :: lkcalc,lkpair,lkedge
 LOGICAL,DIMENSION(:,:,:),ALLOCATABLE :: use_real_part
 PRIVATE
 PUBLIC readbwf,bwfdet_main,bwfdet_setup,bwfdet_wrapper,blip3dgamma_w,blip3dgamma_gl


CONTAINS


 SUBROUTINE readbwf(eionion)
!-------------------------------------------------------------------------!
! READBWF                                                                 !
! =======                                                                 !
! Read blip basis wave function plus geometry from the bwfn.data          !
!                                                                         !
! DA/MDT 3.2001                                                           !
!                                                                         !
! Modifications                                                           !
! =============                                                           !
! NDD  1.2005  Read bwfn.data on each node rather than on master & bcast. !
! William Parker Nov. 2005 - 2006    Modified for CHAMP                   !
! Richard Hennig  4.2008 Read bwfn.data on master & bcast.                !
!-------------------------------------------------------------------------!
 USE mpi_mod
 use all_modules_mod
 IMPLICIT NONE
 INTEGER i,j,k,io,band,ig,ialloc,idum,ispin,num_spins,num_electrons,n1,n2,n3,ierr,nbasisbwf
 INTEGER,DIMENSION(:),ALLOCATABLE :: atno
 REAL(dp) teionion,pa1(3),pa2(3),pa3(3),plane_wave_cutoff,total_energy,&
   &kinetic_energy,local_potential_energy,non_local_potential_energy,&
   &electron_electron_energy
 REAL(dp),INTENT(out) :: eionion
 REAL(dp),DIMENSION(:,:),ALLOCATABLE :: gvecwf,basisbwf
 CHARACTER(80)sline,ltitle,gtitle,dtitle,code,method,functional,pseudo_type,tmpr

#ifdef MPI
 INTEGER, DIMENSION(3) :: blen, indices, types
 INTEGER MPI_bwfsizetype, mpi_err

 TYPE bwfsize_type
    sequence
    logical spin_polarized, pwreal
    integer nbasisbwf, nwvec, nrbwf(3), nkvec_bwfdet
 END TYPE bwfsize_type

 TYPE(bwfsize_type) :: bwf
! Define the MPI type
!  spin_polarized
 blen(1)    = 2;
 indices(1) = 0;
 types(1)   = MPI_LOGICAL;
!  nbasisbwf, nwvec, nrbwf(3), nkvec_bwfdet
 blen(2)    = 6;
 indices(2) = 2*sizeof(MPI_LOGICAL);
 types(2)   = MPI_INTEGER;
!  upper bound
 blen(3)    = 1;
 indices(3) = 2*sizeof(MPI_LOGICAL) +  6*sizeof(MPI_INTEGER);
 types(3)   = MPI_UB;

 CALL MPI_Type_struct(3, blen, indices, types, MPI_bwfsizetype, ierr);
 if(ierr/=0)call errstop('READBWF','Error in MPI_TYPE_STRUCT')
 CALL MPI_TYPE_COMMIT(MPI_bwfsizetype, ierr)
 if(ierr/=0)call errstop('READBWF','Error in MPI_TYPE_COMMIT')
#endif

 if(idtask == 0) then
    ialloc=0

    call open_units(io,ierr)
    if(ierr/=0)call errstop('READBWF','Unable to find free i/o unit')

    open(io,file='bwfn.data',status='old',err=10)

! Rapid scan of file to get relevant array dimensions. Allocate arrays.
    call skip(io,15)
    read(io,*,end=20,err=30)spin_polarized               ; call skip(io,18)
    write(6,'(''spin_polarized='',l2)') spin_polarized   ; call systemflush(6)
    read(io,*,end=20,err=30)nbasisbwf                    ; call skip(io,nbasisbwf+9)
    write(6,'(''nbasisbwf='',i3)') nbasisbwf             ; call systemflush(6)
    read(io,*,end=20,err=30)nwvec                        ; call skip(io,nwvec+2)
    write(6,'(''nwvec='',i6)') nwvec                     ; call systemflush(6)
    read(io,*,end=20,err=30)nrbwf(1),nrbwf(2),nrbwf(3)   ; call skip(io,4)
    write(6,'(''nrbwf='',3i5)') (nrbwf(i),i=1,3)         ; call systemflush(6)
    read(io,*,end=20,err=30)nkvec_bwfdet
    write(6,'(''nkvec_bwfdet='',i5)') nkvec_bwfdet       ; call systemflush(6)
    if(nkvec_bwfdet /= nkvec) then
       write(6,*)'Warning: number of k-vectors in bwfn.data is different than what CHAMP expects'
       write(6,*)'Number of k-points in bwfn.data:',nkvec_bwfdet
       write(6,*)'Number of k-points CHAMP expects:',nkvec
       if(nkvec_bwfdet < nkvec) call errstop('READBWF','Not enough k-points in bwfn.data')
    endif

    allocate(atno(nbasisbwf),basisbwf(3,nbasisbwf),gvecwf(3,nwvec),kvec_bwfdet(3,nkvec_bwfdet),&
         &nband_bwfdet(nkvec_bwfdet,2),boccband(nkvec_bwfdet,2),stat=ialloc)
    if(ialloc/=0)call errstop('READBWF','Allocation problem.')
    num_spins=1 ; if(spin_polarized)num_spins=2

    do k=1,nkvec_bwfdet
       call skip(io,1)
       read(io,*,end=20,err=30)idum,nband_bwfdet(k,1),nband_bwfdet(k,2)
       write(6,'(''idum, nband_bwfdet='',9i5)') k, idum, nband_bwfdet(k,1), nband_bwfdet(k,2) ; call systemflush(6)
       do ispin=1,num_spins
          do band=1,nband_bwfdet(k,ispin)
             if(igrad_lap.eq.0) then
               call skip(io,3)
             else
               call skip(io,4)
             endif
             if(k==1.and.band==1)then
                read(io,fmt='(a)',end=20,err=30)sline
                write(6,'(''ispin, band, sline'',2i5,a)') ispin, band, sline ; call systemflush(6)
                if(scan(sline,",")/=0)then
                   pwreal=.false. ! No inversion symmetry --> complex PW coefficients
                else
                   pwreal=.true.  ! Inversion symmetry --> real PW coefficients
                endif
                backspace(io)
             endif
             if(igrad_lap.eq.0) then
               call skip(io,nrbwf(1)*nrbwf(2)*nrbwf(3))
             else
               call skip(io,2*nrbwf(1)*nrbwf(2)*nrbwf(3)+1)
             endif
          enddo ! bands
       enddo ! spins
    enddo ! k
 endif

#ifdef MPI
! Broadcast array dimensions so other processors start allocating
 if(idtask == 0) then
    bwf = bwfsize_type(spin_polarized,  pwreal, nbasisbwf, nwvec, nrbwf, nkvec_bwfdet)
 endif
 CALL MPI_BCAST(bwf, 1, MPI_bwfsizetype, 0, MPI_COMM_WORLD,ierr)
 if(ierr/=0)call errstop('READBWF','Error in MPI_BCAST')

 if (idtask /= 0) then
    spin_polarized = bwf%spin_polarized
    pwreal         = bwf%pwreal
    nbasisbwf      = bwf%nbasisbwf
    nwvec          = bwf%nwvec
    nrbwf          = bwf%nrbwf
    nkvec_bwfdet   = bwf%nkvec_bwfdet
    allocate(kvec_bwfdet(3,nkvec_bwfdet), nband_bwfdet(nkvec_bwfdet,2),boccband(nkvec_bwfdet,2),stat=ialloc)
    if(ialloc/=0)call errstop('READBWF','Allocation problem.')
    num_spins=1 ; if(spin_polarized)num_spins=2
 endif

! Broadcast nband_bwfdet
 CALL MPI_BCAST(nband_bwfdet, nkvec_bwfdet*2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
 if(ierr/=0)call errstop('READBWF','Error in MPI_BCAST')
#endif

 maxband=maxval(nband_bwfdet(:,:))
 allocate(eigenvalue(maxband,nkvec_bwfdet,2),stat=ialloc)
 if(ialloc/=0)call errstop('READBWF','Eigenvalue allocation problem.')

 select case (igrad_lap)
 case(0)
   if(pwreal)then
      allocate(avc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVC allocation problem 1.')
      avc=0.d0
      if(spin_polarized)then
         allocate(avc2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVC allocation problem 2.')
         avc2=0.d0
      endif
   else
      allocate(cavc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVC allocation problem 3.')
      cavc=0.d0
      if(spin_polarized)then
         allocate(cavc2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVC allocation problem 4.')
         cavc2=0.d0
      endif
   endif
 case(1)
   if(pwreal)then
      allocate(avc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVC allocation problem 1.')
      avc=0.d0
      allocate(avclap(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCLAP allocation problem 1.')
      avclap=0.d0
      if(spin_polarized)then
         allocate(avc2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVC allocation problem 2.')
         avc2=0.d0
         allocate(avclap2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCLAP allocation problem 2.')
         avclap2=0.d0
      endif
   else
      allocate(cavc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVC allocation problem 3.')
      cavc=0.d0
      allocate(cavclap(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCLAP allocation problem 3.')
      cavclap=0.d0
      if(spin_polarized)then
         allocate(cavc2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVC allocation problem 4.')
         cavc2=0.d0
         allocate(cavclap2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCLAP allocation problem 4.')
         cavclap2=0.d0
      endif
   endif
 case(2)
   if(pwreal)then
      allocate(avc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVC allocation problem 1.')
      avc=0.d0
      allocate(avclap(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCLAP allocation problem 1.')
      avclap=0.d0
      allocate(avcgrad1(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCGRAD1 allocation problem 1.')
      avcgrad1=0.d0
      allocate(avcgrad2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCGRAD2 allocation problem 1.')
      avcgrad2=0.d0
      allocate(avcgrad3(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCGRAD3 allocation problem 1.')
      avcgrad3=0.d0
      if(spin_polarized)then
         allocate(avc2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVC allocation problem 2.')
         avc2=0.d0
         allocate(avclap2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCLAP allocation problem 2.')
         avclap2=0.d0
         allocate(avcgrad12(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCGRAD1 allocation problem 2.')
         avcgrad12=0.d0
         allocate(avcgrad22(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCGRAD2 allocation problem 2.')
         avcgrad22=0.d0
         allocate(avcgrad32(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCGRAD3 allocation problem 2.')
         avcgrad32=0.d0
      endif
   else
      allocate(cavc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVC allocation problem 3.')
      cavc=0.d0
      allocate(cavclap(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCLAP allocation problem 3.')
      cavclap=0.d0
      allocate(cavcgrad1(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCGRAD1 allocation problem 3.')
      cavcgrad1=0.d0
      allocate(cavcgrad2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCGRAD2 allocation problem 3.')
      cavcgrad2=0.d0
      allocate(cavcgrad3(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
      if(ialloc/=0)call errstop('READBWF','AVCGRAD3 allocation problem 3.')
      cavcgrad3=0.d0
      if(spin_polarized)then
         allocate(cavc2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVC allocation problem 4.')
         cavc2=0.d0
         allocate(cavclap2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCLAP allocation problem 4.')
         cavclap2=0.d0
         allocate(cavcgrad12(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCGRAD1 allocation problem 4.')
         cavcgrad12=0.d0
         allocate(cavcgrad22(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCGRAD2 allocation problem 4.')
         cavcgrad22=0.d0
         allocate(cavcgrad32(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),stat=ialloc)
         if(ialloc/=0)call errstop('READBWF','AVCGRAD3 allocation problem 4.')
         cavcgrad32=0.d0
      endif
   endif
 end select

 if(idtask == 0) then
    rewind(io)

! Detailed read of bwfn.data:

! Title and basic info about plane-wave DFT calc
    read(io,'(a)',err=30,end=30)ltitle                           ; call skip(io,4)
    read(io,'(a)',err=30,end=30)code                            ; call skip(io,1)
    read(io,'(a)',err=30,end=30)method                          ; call skip(io,1)
    read(io,'(a)',err=30,end=30)functional                      ; call skip(io,1)
    read(io,'(a)',err=30,end=30)pseudo_type                     ; call skip(io,1)
    read(io,*,err=30,end=30)plane_wave_cutoff                   ; call skip(io,1)
    read(io,*,err=30,end=30)spin_polarized                      ; call skip(io,1)
    read(io,*,err=30,end=30)total_energy                        ; call skip(io,1)
    read(io,*,err=30,end=30)kinetic_energy                      ; call skip(io,1)
    read(io,*,err=30,end=30)local_potential_energy              ; call skip(io,1)
    read(io,*,err=30,end=30)non_local_potential_energy          ; call skip(io,1)
    read(io,*,err=30,end=30)electron_electron_energy            ; call skip(io,1)
    read(io,*,err=30,end=30)teionion                            ; call skip(io,1)
    read(io,*,err=30,end=30)num_electrons                       ; call skip(io,4)
    write(6,'(''num_electrons='',i5)') num_electrons            ; call systemflush(6)

    gtitle=ltitle ; dtitle=ltitle

! Geometry
    read(io,*,err=30)nbasisbwf                                     ; call skip(io,1)
    if(nbasisbwf /= ncent) then
       write(6,*)
       write(6,'(''Number of atoms from input:            '',      i10)')ncent
       write(6,'(''Number of atoms from bwfn.data:        '',      i10)')nbasisbwf
       call errstop('READBWF','Number of atoms differs in input and bwfn.data')
    endif
    do i=1,nbasisbwf
       read(io,*,err=30)atno(i),(basisbwf(j,i),j=1,3)
       do j=1,3
          if(cent(j,i)-basisbwf(j,i) > tolerance) then
             write(6,*)
             write(6,'(''Coordinates of atom '',i10,'' in  input:    '')')i
             write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')cent(:,i)
             write(6,'(''Coordinates of atom '',i10,'' in  bwfn.data:'')')i
             write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')basisbwf(:,i)
             call errstop('READBWF','Atom positions differ in input and bwfn.data')
          endif
       enddo
    enddo
    call skip(io,1)
    read(io,*,err=30)pa1
    do j=1,3
       if(rlatt(j,1)-pa1(j) > tolerance) then
          write(6,*)
          write(6,'(''First primitive lattice vector in bwfn.data different from input'')')
          write(6,'(''Primitive lattice vector in  input:    '')')
          write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')rlatt(:,1)
          write(6,'(''Primitive lattice vector in  bwfn.data:'')')
          write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')pa1
          call errstop('READBWF','Lattice vectors differ in input and bwfn.data')
       endif
    enddo
    read(io,*,err=30)pa2
    do j=1,3
       if(rlatt(j,2)-pa2(j) > tolerance) then
          write(6,*)
          write(6,'(''Second primitive lattice vector in bwfn.data different from input'')')
          write(6,'(''Primitive lattice vector in  input:    '')')
          write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')rlatt(:,2)
          write(6,'(''Primitive lattice vector in  bwfn.data:'')')
          write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')pa2
          call errstop('READBWF','Lattice vectors differ in input and bwfn.data')
       endif
    enddo
    read(io,*,err=30)pa3                                        ; call skip(io,4)
    do j=1,3
       if(rlatt(j,3)-pa3(j) > tolerance) then
          write(6,*)
          write(6,'(''Third primitive lattice vector in bwfn.data different from input'')')
          write(6,'(''Primitive lattice vector in  input:    '')')
          write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')rlatt(:,3)
          write(6,'(''Primitive lattice vector in  bwfn.data:'')')
          write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')pa3
          call errstop('READBWF','Lattice vectors differ in input and bwfn.data')
       endif
    enddo

! G vectors
    read(io,*,err=30)nwvec                                      ; call skip(io,1)
    do ig=1,nwvec
       read(io,*,err=30)gvecwf(1,ig),gvecwf(2,ig),gvecwf(3,ig)
    enddo; call skip(io,1)
    read(io,*,err=30)nrbwf
 endif

 if(idtask == 0) then
    call skip(io,4)
! k points, numbers of bands, eigenvalues, orbital coefficients
    read(io,*,err=30)nkvec_bwfdet
    do k=1,nkvec_bwfdet
       call skip(io,1)
       read(io,*,err=30)idum,nband_bwfdet(k,1),nband_bwfdet(k,2),kvec_bwfdet(1,k),kvec_bwfdet(2,k),kvec_bwfdet(3,k)
       write(6,'(/,''idum, nband, kvec='',3i6,3f10.6)') idum,nband_bwfdet(k,1),nband_bwfdet(k,2),kvec_bwfdet(1,k),kvec_bwfdet(2,k),kvec_bwfdet(3,k) ; call systemflush(6)
       if(kvec_bwfdet(1,k)-rkvec(1,k) > tolerance .or. kvec_bwfdet(2,k)-rkvec(2,k) > tolerance .or. kvec_bwfdet(3,k)-rkvec(3,k) > tolerance) then
          write(6,*)
          write(6,'(''K-point '',i4,'' in bwfn.data different from what CHAMP expects'')')k
          write(6,'(''K-point CHAMP expects:   '')')
          write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')rkvec(:,k)
          write(6,'(''K-point in bwfn.data:    '')')
          write(6,'(f12.8,'' '',f12.8,'' '',f12.8)')kvec_bwfdet(:,k)
       endif
! Why do cavc and cavc2 have k-vector indices on them when avc and avc2 do not?
       do ispin=1,num_spins ! 2 if spin_polarized, 1 if not
          do band=1,nband_bwfdet(k,ispin)
             call skip(io,1)
             read(io,*)i,j,eigenvalue(band,k,ispin)                    ; call skip(io,1)
             write(6,'(''kvec, band, eigenvalue'',2i5,f10.6)') i,j,eigenvalue(band,k,ispin) ; call systemflush(6)
             write(6,'(''igrad_lap,nrbwf(1),nrbwf(2),nrbwf(3)='',9i5)') igrad_lap, nrbwf(1),nrbwf(2),nrbwf(3) ; call systemflush(6)
             select case(igrad_lap)
             case(0)
               do n1=0,nrbwf(1)-1
                  do n2=0,nrbwf(2)-1
                     do n3=0,nrbwf(3)-1
                        if(pwreal)then
                           if(ispin==1)then
                              read(io,*,end=20,err=30)avc(n1,n2,n3,band)
                           else
                              read(io,*,end=20,err=30)avc2(n1,n2,n3,band)
                           endif
                           !write(6,*) avc(n1,n2,n3,band) ; call systemflush(6)
                        else
                           if(ispin==1)then
                              read(io,*,end=20,err=30)cavc(n1,n2,n3,band,k)
                           else
                              read(io,*,end=20,err=30)cavc2(n1,n2,n3,band,k)
                           endif
                           !write(6,*) cavc(n1,n2,n3,band,k) ; call systemflush(6)
                        endif
                     enddo
                  enddo
               enddo
             case(1)
               call skip(io,1) !line= 'Orbital'
               write(6,'(''reading orbitals before lap'',i10)') nrbwf(2)*nrbwf(1)*nrbwf(3)
               write(6,*) pwreal ; call systemflush(6)
               do n1=0,nrbwf(1)-1
                  do n2=0,nrbwf(2)-1
                     do n3=0,nrbwf(3)-1
                        if(pwreal)then
                           if(ispin==1)then
                              read(io,*,end=20,err=30)avc(n1,n2,n3,band)
                           else
                              read(io,*,end=20,err=30)avc2(n1,n2,n3,band)
                           endif
                           !write(6,*) avc(n1,n2,n3,band) ; call systemflush(6)
                        else
                           if(ispin==1)then
                              read(io,*,end=20,err=30)cavc(n1,n2,n3,band,k)
                           else
                              read(io,*,end=20,err=30)cavc2(n1,n2,n3,band,k)
                           endif
                           !write(6,*) cavc(n1,n2,n3,band,k) ; call systemflush(6)
                        endif
                     enddo
                  enddo
               enddo
               call skip(io,1) !line= 'Laplacian'
               write(6,'(''reading  lap'')') ; call systemflush(6)
               do n1=0,nrbwf(1)-1
                  do n2=0,nrbwf(2)-1
                     do n3=0,nrbwf(3)-1
                        if(pwreal)then
                           if(ispin==1)then
                              read(io,*,end=20,err=30)avclap(n1,n2,n3,band)
                           else
                              read(io,*,end=20,err=30)avclap2(n1,n2,n3,band)
                           endif
                           !write(6,*) avclap(n1,n2,n3,band) ; call systemflush(6)
                        else
                           if(ispin==1)then
                              !read(io,*,end=20,err=30)cavclap(n1,n2,n3,band,k)
                              read(io,*)cavclap(n1,n2,n3,band,k)
                           else
                              read(io,*,end=20,err=30)cavclap2(n1,n2,n3,band,k)
                           endif
                           !write(6,*) cavclap(n1,n2,n3,band,k) ; call systemflush(6)
                        endif
                     enddo
                  enddo
               enddo
             case(2)
               call skip(io,1) !line= 'Orbital'
               do n1=0,nrbwf(1)-1
                  do n2=0,nrbwf(2)-1
                     do n3=0,nrbwf(3)-1
                        if(pwreal)then
                           if(ispin==1)then
                              read(io,*,end=20,err=30)avc(n1,n2,n3,band)
                           else
                              read(io,*,end=20,err=30)avc2(n1,n2,n3,band)
                           endif
                        else
                           if(ispin==1)then
                              read(io,*,end=20,err=30)cavc(n1,n2,n3,band,k)
                           else
                              read(io,*,end=20,err=30)cavc2(n1,n2,n3,band,k)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               call skip(io,1) !line= 'Laplacian'
               do n1=0,nrbwf(1)-1
                  do n2=0,nrbwf(2)-1
                     do n3=0,nrbwf(3)-1
                        if(pwreal)then
                           if(ispin==1)then
                              read(io,*,end=20,err=30)avclap(n1,n2,n3,band)
                           else
                              read(io,*,end=20,err=30)avclap2(n1,n2,n3,band)
                           endif
                        else
                           if(ispin==1)then
                              read(io,*,end=20,err=30)cavclap(n1,n2,n3,band,k)
                           else
                              read(io,*,end=20,err=30)cavclap2(n1,n2,n3,band,k)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               call skip(io,1) !line= 'Gradient - a1 direction'
               do n1=0,nrbwf(1)-1
                  do n2=0,nrbwf(2)-1
                     do n3=0,nrbwf(3)-1
                        if(pwreal)then
                           if(ispin==1)then
                              read(io,*,end=20,err=30)avcgrad1(n1,n2,n3,band)
                           else
                              read(io,*,end=20,err=30)avcgrad12(n1,n2,n3,band)
                           endif
                        else
                           if(ispin==1)then
                              read(io,*,end=20,err=30)cavcgrad1(n1,n2,n3,band,k)
                           else
                              read(io,*,end=20,err=30)cavcgrad12(n1,n2,n3,band,k)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               call skip(io,1) !line= 'Gradient - a2 direction'
               do n1=0,nrbwf(1)-1
                  do n2=0,nrbwf(2)-1
                     do n3=0,nrbwf(3)-1
                        if(pwreal)then
                           if(ispin==1)then
                              read(io,*,end=20,err=30)avcgrad2(n1,n2,n3,band)
                           else
                              read(io,*,end=20,err=30)avcgrad22(n1,n2,n3,band)
                           endif
                        else
                           if(ispin==1)then
                              read(io,*,end=20,err=30)cavcgrad2(n1,n2,n3,band,k)
                           else
                              read(io,*,end=20,err=30)cavcgrad22(n1,n2,n3,band,k)
                           endif
                        endif
                     enddo
                  enddo
               enddo
               call skip(io,1) !Line= 'Gradient - a3 direction'
               do n1=0,nrbwf(1)-1
                  do n2=0,nrbwf(2)-1
                     do n3=0,nrbwf(3)-1
                        if(pwreal)then
                           if(ispin==1)then
                              read(io,*,end=20,err=30)avcgrad3(n1,n2,n3,band)
                           else
                              read(io,*,end=20,err=30)avcgrad32(n1,n2,n3,band)
                           endif
                        else
                           if(ispin==1)then
                              read(io,*,end=20,err=30)cavcgrad3(n1,n2,n3,band,k)
                           else
                              read(io,*,end=20,err=30)cavcgrad32(n1,n2,n3,band,k)
                           endif
                        endif
                     enddo
                  enddo
               enddo
             end select
          enddo ! bands
       enddo ! spin states
    enddo ! k
    if(.not.spin_polarized)eigenvalue(:,:,2)=eigenvalue(:,:,1)
    write(6,*)
    close(io) ; open_unit(io)=.false.

    write(6,'(t2,''Title: '',a)') trim(ltitle)
    write(6,'(t2,''Generating code                           : '',a)') trim(code)
    write(6,'(t2,''Method                                    : '',a)') trim(method)
    write(6,'(t2,''DFT functional                            : '',a)') trim(functional)
    write(6,'(t2,''Pseudopotential type                      : '',a)') trim(pseudo_type)
    tmpr=r2s(plane_wave_cutoff,'(f12.3)')
    write(6,'(t2,''Plane-wave cutoff (au)                    : '',a)') trim(tmpr)
    write(6,'(t2,''Smoothing B-spline grid:                  : '',3i5)') nrbwf(1:3)
    write(6,*)
    write(6,*)'Number of k points                        : ', trim(i2s(nkvec_bwfdet))
    write(6,*)'Max # bands per k point                   : ', trim(i2s(maxband))
    write(6,*)'Number of G vectors                       : ', trim(i2s(nwvec))
    write(6,*)
    write(6,*)'DFT ENERGY AND COMPONENTS (au per primitive cell):'
    tmpr=r2s(total_energy,'(f16.10)')
    write(6,'(t2,''Total energy                              : '',a)') trim(tmpr)
    tmpr=r2s(kinetic_energy,'(f16.10)')
    write(6,'(t2,''Kinetic energy                            : '',a)') trim(tmpr)
    tmpr=r2s(local_potential_energy,'(f16.10)')
    write(6,'(t2,''Local potential energy                    : '',a)') trim(tmpr)
    tmpr=r2s(non_local_potential_energy,'(f16.10)')
    write(6,'(t2,''Non-local potential energy                : '',a)') trim(tmpr)
    tmpr=r2s(electron_electron_energy,'(f16.10)')
    write(6,'(t2,''Electron-electron energy                  : '',a)') trim(tmpr)
    tmpr=r2s(teionion,'(f16.10)')
    write(6,'(t2,''Ion-ion energy                            : '',a)') trim(tmpr)
    write(6,*)
    if(pwreal)then
       write(6,*)'Real blip coefficients ==> GAMMA calculation'
    else
       write(6,*)'Complex blip coefficients ==> calculation with K-POINTS'
    endif
    eionion=teionion
 endif

#ifdef MPI
 CALL MPI_BCAST(eigenvalue, maxband*nkvec_bwfdet*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
 CALL MPI_BCAST(kvec_bwfdet, 3*nkvec_bwfdet, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
! Not needed, read from input later
! CALL MPI_BCAST(eionion, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

! Broadcast coefficients
 select case(igrad_lap)
 case (0)
   if(pwreal)then
      CALL MPI_BCAST(avc, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (spin_polarized) then
         CALL MPI_BCAST(avc2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      endif
   else
!!!   CALL MPI_BCAST(cavc(0,0,0,1,1), nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(cavc, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      if (spin_polarized) then
         write(6,'(''starting bcast of cavc2'')')
         CALL MPI_BCAST(cavc2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      endif
   endif
 case (1)
   if(pwreal)then
      CALL MPI_BCAST(avc, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(avclap, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (spin_polarized) then
         CALL MPI_BCAST(avc2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(avclap2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      endif
   else
!!!   CALL MPI_BCAST(cavc(0,0,0,1,1), nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(cavc, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(cavclap, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      if (spin_polarized) then
         write(6,'(''starting bcast of cavc2'')')
         CALL MPI_BCAST(cavc2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(cavclap2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      endif
   endif
 case (2)
   if(pwreal)then
      CALL MPI_BCAST(avc, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(avclap, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(avcgrad1, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(avcgrad2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(avcgrad3, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (spin_polarized) then
         CALL MPI_BCAST(avc2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(avclap2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(avcgrad12, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(avcgrad22, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(avcgrad32, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      endif
   else
!!!   CALL MPI_BCAST(cavc(0,0,0,1,1), nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(cavc, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(cavclap, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(cavclap, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(cavcgrad2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(cavcgrad3, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      if (spin_polarized) then
         write(6,'(''starting bcast of cavc2'')')
         CALL MPI_BCAST(cavc2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(cavclap2, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(cavcgrad12, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(cavcgrad22, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(cavcgrad32, nrbwf(1)*nrbwf(2)*nrbwf(3)*maxband*nkvec_bwfdet, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      endif
   endif
 end select
#endif
 write(6,'(''end of readbwf'')')

 return

10 call errstop('READBWF','Cannot find bwfn.data file.')
20 call errstop('READBWF','Read past end of bwfn.data file.')
30 call errstop('READBWF','Error reading bwfn.data file.')
 END SUBROUTINE readbwf


 SUBROUTINE bwfdet_setup(orb_norm)
!----------------------------------------------------------------------------!
! Setup for BLIP calculation                                                 !
! (adapted from PWFDET_SETUP MDT 3.2001)                                     !
!                                                                            !
! NB: This routine will shortly be replaced when the standard PWFDET_SETUP   !
! will be used for both PWs and blips. MDT 3.2002                            !
!----------------------------------------------------------------------------!
 use all_modules_mod
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: orb_norm
 INTEGER, PARAMETER :: mdet_max_mods=15,lsize=500,num_g=lsize+lsize-1, &
  &twonum_g=num_g*2,sixnum_g=num_g*6,fournum_gp2=twonum_g+twonum_g+2
 INTEGER ialloc,ik,jk,ig,ispin,eff_nele(2),ix,iy,iz,num_spins,band,k,        &
  &ngridpoints,i,nstates,ne,                 &
  &i1m,i2m,i3m,i1,i2,i3,j1,j2,j3,jmax,kmax,n,l,m,istart,&
  &nep,nseg,ngpri,j,nentry,iband

 INTEGER,ALLOCATABLE :: indx(:),kcheck(:),wf_np(:,:),wf_nm(:,:),&
  &istart_array(:),nentry_array(:),temp_i(:),vector_index(:),jm(:)
 REAL(dp) ksum(3),rvec(3),temp1,temp2,temp3,k_s(3),ktemp(1:3),             &
  &average_real_part,average_imaginary_part,pb1(3),pb2(3),pb3(3),              &
  &b1(3),b2(3),b3(3),           &
  &a1_bwfdet(3),a2_bwfdet(3),a3_bwfdet(3),           &
  &rkvec_shift_latt,glatt_inv_squared(3),r,z,low,high,current,           &
  &low2,current2,mean,glatt_sim_inv_squared(3)
 REAL(dp),ALLOCATABLE,DIMENSION(:) :: eigtemp,temp_r,low_array,high_array,d
 REAL(dp),ALLOCATABLE,DIMENSION(:,:) :: pr_lattice,sr_lattice
 COMPLEX(dp),ALLOCATABLE :: sum_orbs(:,:)
 LOGICAL :: metal=.false.
 LOGICAL,ALLOCATABLE :: ltemp(:)

 INTEGER :: irealimag
 REAL(dp) :: real_part,real_coefficient_sum,abs_real_coefficient_sum

 REAL(dp),ALLOCATABLE,DIMENSION(:,:) :: rpsi,lap
 REAL(dp) :: eps
 COMPLEX(dp),ALLOCATABLE,DIMENSION(:,:) :: rpsi_c,lap_c

 common /periodic2/ rkvec_shift_latt(3)

!Set local variables equal to common block equivalents

 painv=transpose(rlatt_inv)
 pb1=glatt(:,1)
 pb2=glatt(:,2)
 pb3=glatt(:,3)
 b1=glatt_sim(:,1)
 b2=glatt_sim(:,2)
 b3=glatt_sim(:,3)
 ndet_bwfdet=ndet

!generate lattice primitive reciprocal points (in Cartesian coordinates)

 allocate(istart_array(twonum_g),nentry_array(twonum_g),temp_i(sixnum_g),&
         &vector_index(twonum_g),temp_r(sixnum_g),low_array(twonum_g),&
         &high_array(twonum_g),d(sixnum_g),jm(twonum_g),pr_lattice(3,num_g),&
         &stat=ialloc)
 if(ialloc/=0)call errstop('READBWF','Lattice array allocation problem')

 z=0.d0
 ngpri=1
 glatt_inv_squared(3)=0.d0
 glatt_inv_squared(2)=0.d0
 glatt_inv_squared(1)=0.d0

 glatt_inv_squared(3)=sqrt(glatt_inv(1,3)**2+glatt_inv(2,3)**2+ &
                      &glatt_inv(3,3)**2)
 z=glatt_inv_squared(3)

 glatt_inv_squared(2)=sqrt(glatt_inv(1,2)**2+glatt_inv(2,2)**2+ &
                      &glatt_inv(3,2)**2)
 z=max(z,glatt_inv_squared(2))

 glatt_inv_squared(1)=sqrt(glatt_inv(1,1)**2+glatt_inv(2,1)**2+ &
                      &glatt_inv(3,1)**2)
 z=max(z,glatt_inv_squared(1))

 r=twonum_g
 i=int(r**(1.d0/3.d0))
 i=(i-1)/2
 r=i/z
 z=r+tolerance

4 i1m=int(z*glatt_inv_squared(1))
  i2m=int(z*glatt_inv_squared(2))
  i3m=int(z*glatt_inv_squared(3))
 if(((i1m+i1m+1)*(i2m+i2m+1)*(i3m+i3m+1))>=fournum_gp2)goto 5
 i1=i1m ; i2=i2m ; i3=i3m
 z=z+r
 goto 4

5 r=r*0.5d0
 z=z-r
 if(r>=(0.2d0))goto 4

 i2m=-i2
 i3m=-i3
 n=0
 jmax=i2

 do j1=-i1,0
   if(j1==0)jmax=0
   kmax=i3
   do j2=i2m,jmax
     if(abs(j2)==j1)kmax=-1
      do j3=i3m,kmax
        n=n+3
        temp_i(n-2)=j1
        temp_i(n-1)=j2
        temp_i(n)=j3
      enddo
   enddo
 enddo

 do l=1,n
   d(l)=temp_i(l)
   temp_r(l)=0.d0
 enddo
 n=n/3

 m=0
 do i=1,n
    m=m+3
    temp_r(m-2)=d(m-2)*glatt(1,1)+d(m-1)*glatt(1,2)+d(m)*glatt(1,3)
    temp_r(m-1)=d(m-2)*glatt(2,1)+d(m-1)*glatt(2,2)+d(m)*glatt(2,3)
    temp_r(m)  =d(m-2)*glatt(3,1)+d(m-1)*glatt(3,2)+d(m)*glatt(3,3)
 enddo
 m=0

 do l=1,n
   m=m+3
   d(l)=temp_r(m-2)**2+temp_r(m-1)**2+temp_r(m)**2
 enddo
 jmax=min(n,lsize-1)

!sort the temp_r vectors by their squared length d
 call gvsort(d,vector_index,low_array,high_array,istart_array,&
            &nentry_array,jm,n,jmax,tolerance)

 pr_lattice=0.d0
 do m=1,jmax
   j=vector_index(m)
   n=ngpri+1
   ngpri=ngpri+2
   j=j*3-3
   do i=1,3
     j=j+1
     pr_lattice(i,n)=temp_r(j)
     pr_lattice(i,ngpri)=-temp_r(j)
   enddo
 enddo

 deallocate(istart_array,nentry_array,temp_i,vector_index,temp_r,low_array,&
           &high_array,d,jm)
!end of pr_lattice setup

!generate lattice simulation reciprocal points (in Cartesian coordinates)

 allocate(istart_array(twonum_g),nentry_array(twonum_g),temp_i(sixnum_g),&
         &vector_index(twonum_g),temp_r(sixnum_g),low_array(twonum_g),&
         &high_array(twonum_g),d(sixnum_g),jm(twonum_g),sr_lattice(3,num_g),&
         &stat=ialloc)
 if(ialloc/=0)call errstop('READBWF','Lattice array allocation problem')

 ngpri=1
 glatt_sim_inv_squared(3)=0.d0
 glatt_sim_inv_squared(2)=0.d0
 glatt_sim_inv_squared(1)=0.d0

 glatt_sim_inv_squared(3)=sqrt(glatt_sim_inv(1,3)**2+glatt_sim_inv(2,3)**2+ &
                          &glatt_sim_inv(3,3)**2)
 z=glatt_sim_inv_squared(3)

 glatt_sim_inv_squared(2)=sqrt(glatt_sim_inv(1,2)**2+glatt_sim_inv(2,2)**2+ &
                          &glatt_sim_inv(3,2)**2)
 z=max(z,glatt_sim_inv_squared(2))

 glatt_sim_inv_squared(1)=sqrt(glatt_sim_inv(1,1)**2+glatt_sim_inv(2,1)**2+ &
                          &glatt_sim_inv(3,1)**2)
 z=max(z,glatt_sim_inv_squared(1))

 r=twonum_g
 i=int(r**(1.d0/3.d0))
 i=(i-1)/2
 r=i/z
 z=r+tolerance

8 i1m=int(z*glatt_sim_inv_squared(1))
  i2m=int(z*glatt_sim_inv_squared(2))
  i3m=int(z*glatt_sim_inv_squared(3))
 if(((i1m+i1m+1)*(i2m+i2m+1)*(i3m+i3m+1))>=fournum_gp2)goto 9
 i1=i1m ; i2=i2m ; i3=i3m
 z=z+r
 goto 8

9 r=r*0.5d0
 z=z-r
 if(r>=(0.2d0))goto 8

 i2m=-i2
 i3m=-i3
 n=0
 jmax=i2

 do j1=-i1,0
   if(j1==0)jmax=0
   kmax=i3
   do j2=i2m,jmax
     if(abs(j2)==j1)kmax=-1
      do j3=i3m,kmax
        n=n+3
        temp_i(n-2)=j1
        temp_i(n-1)=j2
        temp_i(n)=j3
      enddo
   enddo
 enddo


 do l=1,n
   d(l)=temp_i(l)
   temp_r(l)=0.d0
 enddo
 n=n/3

 m=0
 do i=1,n
    m=m+3
    temp_r(m-2)=d(m-2)*glatt_sim(1,1)+d(m-1)*glatt_sim(1,2)+d(m)*glatt_sim(1,3)
    temp_r(m-1)=d(m-2)*glatt_sim(2,1)+d(m-1)*glatt_sim(2,2)+d(m)*glatt_sim(2,3)
    temp_r(m)=d(m-2)*glatt_sim(3,1)+d(m-1)*glatt_sim(3,2)+d(m)*glatt_sim(3,3)
 enddo
 m=0

 do l=1,n
   m=m+3
   d(l)=temp_r(m-2)**2+temp_r(m-1)**2+temp_r(m)**2
 enddo
 jmax=min(n,lsize-1)

!sort the temp_r vectors by their squared length d
 call gvsort(d,vector_index,low_array,high_array,istart_array,&
            &nentry_array,jm,n,jmax,tolerance)

 sr_lattice=0.d0
 do m=1,jmax
   j=vector_index(m)
   n=ngpri+1
   ngpri=ngpri+2
   j=j*3-3
   do i=1,3
     j=j+1
     sr_lattice(i,n)=temp_r(j)
     sr_lattice(i,ngpri)=-temp_r(j)
   enddo
 enddo

 deallocate(istart_array,nentry_array,temp_i,vector_index,temp_r,low_array,&
           &high_array,d,jm)
!end of sr_lattice setup

!lattice vectors
 a1_bwfdet=rlatt_sim(:,1)
 a2_bwfdet=rlatt_sim(:,2)
 a3_bwfdet=rlatt_sim(:,3)

! Allocate workspace.
 allocate(gmap(3,nwvec),kdotr(nkvec_bwfdet),zdum(nkvec_bwfdet),lzdum(nkvec_bwfdet),gzdum(3,nkvec_bwfdet),  &
  &ztemp(nwvec),kvec2(3,nkvec_bwfdet),stat=ialloc)
 if(ialloc/=0)call errstop('BWFDET_SETUP','Allocation problem <1>')

 allocate (lkcalc(nkvec_bwfdet),lkpair(nkvec_bwfdet),iprom_repl_idx(mdet_max_mods,ndet_bwfdet,nspin)&
  &,iadd_idx(mdet_max_mods,ndet_bwfdet,nspin),isub_idx(2,mdet_max_mods,ndet_bwfdet,nspin)    &
  &,ltemp(nkvec_bwfdet),lkedge(nkvec_bwfdet),use_real_part(maxband,nkvec_bwfdet,nspin)              &
  &,sum_orbs(maxband,nkvec_bwfdet),stat=ialloc)
 if(ialloc/=0)call errstop('BWFDET_SETUP','Allocation problem <2>')

! Added allocation for CHAMP modifications
 allocate(wf_np(ndet_bwfdet,nspin),wf_nm(ndet_bwfdet,nspin))
 if(ialloc/=0)call errstop('BWFDET_SETUP','Allocation problem <3>')

 wf_np=0
 wf_nm=0

! Zero excited state stuff
 iprom_repl_idx=0 ; iadd_idx=0 ; isub_idx=0

! Do we have different orbitals for different spins (num_spins=2) or
! not (num_spins=1)?
 num_spins=1 ; if(spin_polarized)num_spins=2

! We need to impose constraints on the input k point set in order that
! we can make the orbitals real through linear combinations of Psi and Psi*
! whilst ensuring the linear combinations fit inside the supercell and thus
! give the right energy when integrated.

! Conditions:
! (1) When reduced into the supercell Brillouin zone, the primitive cell
!     k points map onto a unique single k_s vector which may be (0,0,0) or any
!     G_s/2 vector.
! (2) For non-equivalent k points (i.e. points not separated by a reciprocal
!     lattice vector) then it is allowed to have data both at +k and -k
!     in the pwfn.data file. However, one of these will be selected
!     at random and two states will be created from the real and imaginary
!     parts of the data for that k point, and its paired k will be ignored.
!     It is therefore encouraged that only 1 k point per pair is present.

! Flag the unnecessary k points we don't need because they are part
! of a k=k'+G pair
 lkcalc(:)=.true.
 lkpair(:)=.false.
 do ik=1,nkvec_bwfdet
  if(lkcalc(ik))then
aloop: do jk=ik+1,nkvec_bwfdet
    ksum(1:3)=kvec_bwfdet(1:3,ik)+kvec_bwfdet(1:3,jk)
    do ig=1,num_g
     if(abs(ksum(1)-pr_lattice(1,ig))<tolerance.and.&
      &abs(ksum(2)-pr_lattice(2,ig))<tolerance.and.&
      &abs(ksum(3)-pr_lattice(3,ig))<tolerance)then
      lkcalc(jk)=.false.
      lkpair(ik)=.true. ; lkpair(jk)=.true.
      exit aloop
     endif
    enddo ! G
   enddo aloop
  endif
 enddo ! ik

 if(any(.not.lkcalc))then
  call errwarn('BWFDET_SETUP','Redundant paired k points in bwfn.data file.')
 endif

! Flag the k points where +k and -k are equivalent (i.e. they differ by
! a primitive cell reciprocal lattive vector.)
!
! These are important because we can only create one state from combining
! Psi and Psi*, and we may thus choose the real or imaginary part according
! to our whim.
!
! An important complication is that such orbitals may be pure real everywhere,
! pure imaginary everywhere, or complex (unlike orbitals whose k is not
! equivalent to its -k, which cannot be pure real or pure imaginary
! everywhere). It is thus important that we don't e.g. take the
! real bit if the orbital is pure imaginary (see later)
!
 lkedge=.false.
 do ik=1,nkvec_bwfdet
  if(lkcalc(ik))then
   ksum(1:3)=2.d0*kvec_bwfdet(1:3,ik)
   do ig=1,num_g
    if(abs(ksum(1)-pr_lattice(1,ig))<tolerance.and.&
     &abs(ksum(2)-pr_lattice(2,ig))<tolerance.and.&
     &abs(ksum(3)-pr_lattice(3,ig))<tolerance)then
     lkedge(ik)=.true.
     exit
    endif
   enddo ! G
  endif
 enddo ! k

! Sort out which states to occupy with electrons (i.e. the nele(spin) states
! with the lowest eigenvalues for each spin.)

 i=max(sum(nband_bwfdet(:,1)),sum(nband_bwfdet(:,2)))
 allocate(eigtemp(i),kcheck(i),indx(i),stat=ialloc)
 if(ialloc/=0)call errstop('BWFDET_SETUP','EIGTEMP allocation.')

!we have only two possible spins in CHAMP
 allocate(nele(2))
 nele(1)=nup
 nele(2)=ndn
 nemax=maxval(nele(:))

 boccband(:,:)=0

 do ispin=1,num_spins
  if(nele(ispin)==0)cycle
  nstates=0
  do ik=1,nkvec_bwfdet
   if(lkcalc(ik))then
    do iband=1,nband_bwfdet(ik,ispin)
     nstates=nstates+1
     eigtemp(nstates)=eigenvalue(iband,ik,ispin)
     kcheck(nstates)=ik
    enddo
   endif
  enddo

  call indexx(eigtemp(1:nstates),indx(1:nstates)) ! creates index table in
                                                  ! indx array

! Add up the number of states at each k that are among the nele(spin) lowest
! energy states (taking into account that points away from the BZ edge create
! two states at +k and -k ) and stick the result in boccband(k,spin)
  ne=0
  do i=1,nstates
   ik=kcheck(indx(i))
   boccband(ik,ispin)=boccband(ik,ispin)+1
   ne=ne+1
   if(.not.lkedge(ik))ne=ne+1
   if(ne>=nele(ispin))exit
  enddo

  if(idtask == 0.and.ne>nele(ispin))then
   call errstop('BWFDET_SETUP','Problem in orbital counting. Not allowed to &
    &singly occupy +k/-k paired states - need more intelligent code?')
  endif

  if(idtask == 0.and.nele(ispin)<nstates)then
   if(abs(eigtemp(indx(nele(ispin)))-eigtemp(indx(nele(ispin)+1)))<1.d-6)then
    call errwarn('BWFDET_SETUP','Partially occupied degenerate states at the &
     &Fermi level. Multideterminant calculation probably advisable.')
   endif
  endif

 enddo ! num_spins

 if(.not.spin_polarized)boccband(:,2)=boccband(:,1)

 if(idtask == 0)then

  write(6,*)
  write(6,'( 1x,''Blip setup'',/1x,''=========='')')

s:do ispin=1,num_spins
   do ik=2,nkvec_bwfdet
    if(lkcalc(ik))then
     if(boccband(ik,ispin)/=boccband(1,ispin))then
      metal=.true.
      exit s
     endif
    endif
   enddo
  enddo s

  if(metal)then
   write(6,*)
   write(6,*)'METALLIC STATE DETECTED'
   if(.not.spin_polarized)then
    write(6,*)'Number of doubly-occupied bands filled at each k point : '
    write(6,*)'    k    nband_bwfdet'
    do ik=1,nkvec_bwfdet
     write(6,'(1x,a,1x,a)')trim(i2s(ik)),trim(i2s(boccband(ik,1)))
    enddo
   else
    write(6,*)'Number of singly-occupied bands filled at each k point : '
    write(6,*)'UP SPIN'
    write(6,*)'    k    nband_bwfdet'
    do ik=1,nkvec_bwfdet
     write(6,'(1x,a,1x,a)')trim(i2s(ik)),trim(i2s(boccband(ik,1)))
    enddo
    write(6,*)'DOWN SPIN'
    write(6,*)'    k    nband_bwfdet'
    do ik=1,nkvec_bwfdet
     write(6,'(1x,a,1x,a)')trim(i2s(ik)),trim(i2s(boccband(ik,2)))
    enddo
   endif
  else
   write(6,*)
   write(6,*)'INSULATING STATE DETECTED'
   if(.not.spin_polarized)then
    write(6,*)'Number of doubly-occupied bands filled at each k point : ', &
     &trim(i2s(boccband(1,1)))
   else
    write(6,*)'Number of singly-occupied bands filled at each k point : ', &
     &trim(i2s(sum(boccband(1,1:2))))
   endif
  endif
  write(6,*)

 endif

 deallocate(indx,eigtemp,eigenvalue,kcheck)
! An important complication with orbitals whose k is on the Brillouin zone
! edge is that such orbitals may be pure real everywhere, pure
! imaginary everywhere, or complex (unlike orbitals whose k is not
! equivalent to its -k, which cannot be pure real or pure imaginary
! everywhere). It is thus important that we don't e.g. take the
! real bit if the orbital is pure imaginary (see later)
!
! So, generate grid (over offset cube; avoid origin and include all octants).
! Accumulate real_squared and imaginary_squared values of orbitals at each
! grid point in the sum_orbs(band,k) vector, summing over grid points as
! we go.

! Changed by WDP to reflect CHAMP current method of choosing real or imaginary
! with plane waves (in read_orb_pw_pwscf and read_orb_pw_tm)
! Sums over real part of coefficients and divide by the sum of the absolute 
! values -- if absolute value greater than 10^-6, use real part
  use_real_part=.false.
  irealimag=0
! If we have only gamma then all this stuff can be skipped
  if(.not.pwreal)then
   do ispin=1,num_spins ! i.e. 2 for spin-polarized systems, 1 otherwise
    do k=1,nkvec_bwfdet
     if(lkedge(k))then
      irealimag=2
      do band=1,boccband(k,ispin)
        if(ispin==1)then
          real_coefficient_sum=sum(dble(cavc(:,:,:,band,k)))
          abs_real_coefficient_sum=sum(abs(dble(cavc(:,:,:,band,k))))
          real_part=abs(real_coefficient_sum/abs_real_coefficient_sum)
          if(real_part > 1.d-6) then
            use_real_part(band,k,ispin)=.true.
          endif
          if(use_real_part(band,k,ispin)) irealimag=1
          write(6,'(''ikv,iband,ireal_imag,sum,sum_abs='',3i4,9d12.4)')& 
             &k,band,irealimag,real_coefficient_sum,abs_real_coefficient_sum
        else
          real_coefficient_sum=sum(dble(cavc2(:,:,:,band,k)))
          abs_real_coefficient_sum=sum(abs(dble(cavc2(:,:,:,band,k))))
          real_part=real_coefficient_sum/abs_real_coefficient_sum
          if(real_part > 1.d-6) then
            use_real_part(band,k,ispin)=.true.
          endif
          if(use_real_part(band,k,ispin)) irealimag=2
          write(6,'(''ikv,iband,ireal_imag,sum,sum_abs='',3i4,9d12.4)')& 
             &k,band,irealimag,real_coefficient_sum,abs_real_coefficient_sum
        endif !ispin 1 or 2
      enddo ! occupied bands at k
     endif
    enddo ! k
   enddo ! spins
   if(.not.spin_polarized)use_real_part(:,:,2)=use_real_part(:,:,1)

  !next step in deciding between real and imaginary parts
!loop over ikv (kpoint index), set jorba
!call get_real_and_imaginary_parts_of_blips for orb_si and ddorb_si
!how do I choose r???
!it looks like it's just /0.1, 0.2, 0.3/

!  allocate(rpsi(nemax,nkvec_bwfdet),lap(nemax,nkvec_bwfdet))
!  allocate(rpsi_c(nemax,nkvec_bwfdet),lap_c(nemax,nkvec_bwfdet))
!
!  rvec(1)=0.1d0
!  rvec(2)=0.2d0
!  rvec(3)=0.3d0
!
!  eps=1.d-4
!
!  do ispin=1,num_spins ! i.e. 2 for spin-polarized systems, 1 otherwise
!   call get_real_and_imaginary_parts_of_blips(rvec,rpsi_c,lap_c,ispin)
!   do k=1,nkvec_bwfdet
!     do band=1,boccband(k,ispin)
!        write(6,'(''real orb_si,ddorb_si='',i5,9d12.4)') &
!             & band,dble(rpsi_c(band,k)),dble(lap_c(band,k))
!        write(6,'(''imag orb_si,ddorb_si='',i5,9d12.4)') &
!             & band,aimag(rpsi_c(band,k)),aimag(lap_c(band,k))
!     end do!band
!   end do!k
!   rpsi=dble(rpsi_c)-aimag(rpsi_c)
!   lap=dble(lap_c)-aimag(lap_c)
!   do k=1,nkvec_bwfdet
!    if(lkedge(k)) then
!      do band=1,boccband(k,ispin)
!       if(band <= boccband(k,ispin)) then
!         if(band >= 3) then
!           if(abs(rpsi(band,k)*lap(band-2,k)/&
!                &(rpsi(band-2,k)*lap(band,k))-1) < eps) then
!             continue
!           else if(abs(rpsi(band-1,k)*lap(band-2,k)/&
!                     &(rpsi(band-2,k)*lap(band-1,k))-1) < eps) then
!            use_real_part(band,k,ispin)=.true.
!            continue
!           endif
!         else if(band >= 4) then
!           if(abs(rpsi(band,k)*lap(band-3,k)/&
!                &(rpsi(band-3,k)*lap(band,k))-1) < eps) then
!             continue
!           else if(abs(rpsi(band-1,k)*lap(band-3,k)/&
!                     &(rpsi(band-3,k)*lap(band-1,k))-1) < eps) then
!             use_real_part(band,k,ispin)=.true.
!             continue
!           endif
!        endif ! band >= 2 or 3
!        if(band >=2) then
!          if(abs(rpsi(band-1,k)/rpsi(band,k)) < 1.d0) then
!            use_real_part(band,k,ispin)=.true.
!          endif
!        endif
!      endif ! band+1 <= boccbands(band,k,ispin)
!     enddo ! bands
!    endif ! lkedge
!   enddo ! k
!  enddo ! spins
!
!  deallocate(rpsi)
!  deallocate(lap)
!  deallocate(rpsi_c)
!  deallocate(lap_c)

   do ispin=1,num_spins
    do k=1,nkvec_bwfdet
     write(6,*)'For k-point ',k
     do band=1,boccband(k,ispin)
      if(use_real_part(band,k,ispin)) then
         write(6,*)'Using real part of band ',band
      else
         write(6,*)'Using imaginary part of band ',band
      end if
     enddo !band
    enddo ! k
   enddo ! spins
        
  endif ! pwreal


!  do ispin=1,num_spins ! i.e. 2 for spin-polarized systems, 1 otherwise
!   ngridpoints=0
!   sum_orbs=0.d0
!   do ix=-3,3
!    do iy=-3,3
!     do iz=-3,3
!      rvec(1)=dble(ix)*12.3456d0-0.01234d0
!      rvec(2)=dble(iy)*12.3456d0-0.04567d0
!      rvec(3)=dble(iz)*12.3456d0-0.08910d0
!      call sum_orbs_over_grid_with_blips(rvec,sum_orbs,ispin)
!      ngridpoints=ngridpoints+1
!     enddo
!    enddo
!   enddo
!   do k=1,nkvec_bwfdet
!    if(lkedge(k))then
!     do band=1,boccband(k,ispin)
!      average_real_part=sqrt(dble(sum_orbs(band,k))/(dble(ngridpoints)))
!      average_imaginary_part=sqrt(aimag(sum_orbs(band,k))/(dble(ngridpoints)))
!      if(average_real_part>average_imaginary_part)then
!       use_real_part(band,k,ispin)=.true.
!      endif
!     enddo ! occupied bands at k
!    endif
!   enddo ! k
!  enddo ! spins
!  if(.not.spin_polarized)use_real_part(:,:,2)=use_real_part(:,:,1)
! endif

! Generate potential k_s vector from each k point in pwfn.data
 kvec2=1.d5
 ltemp=.true.
 do jk=1,nkvec_bwfdet
bloop:do ik=1,num_g
   ktemp(1:3)=kvec_bwfdet(1:3,jk)-0.5d0*sr_lattice(1:3,ik)
   temp1=one_over_twopi*dot_product(ktemp,a1_bwfdet)
   temp2=one_over_twopi*dot_product(ktemp,a2_bwfdet)
   temp3=one_over_twopi*dot_product(ktemp,a3_bwfdet)
   if(abs(temp1-anint(temp1))<tolerance.and.abs(temp2-anint(temp2))<tolerance.and. &
    &abs(temp3-anint(temp3))<tolerance)then ! i.e. if projections are integer
    ltemp(jk)=.false.
    kvec2(1:3,jk)=0.5*sr_lattice(1:3,ik)
    exit bloop
   endif
  enddo bloop
 enddo

 k_s=kvec2(1:3,1)

 if(idtask == 0)then

! Complain if no unique k_s
  do ik=2,nkvec_bwfdet
   if(abs(k_s(1)-kvec2(1,ik))>tolerance.or.abs(k_s(2)-kvec2(2,ik))>tolerance.or. &
    &abs(k_s(3)-kvec2(3,ik))>tolerance)then
    write(6,*)'Primitive cell reciprocal lattice vectors: (au)'
    write(6,'(a,f16.10,1x,f16.10,1x,f16.10)')'pb1',pb1(1:3)
    write(6,'(a,f16.10,1x,f16.10,1x,f16.10)')'pb2',pb2(1:3)
    write(6,'(a,f16.10,1x,f16.10,1x,f16.10)')'pb3',pb3(1:3)
    write(6,*)'Simulation cell reciprocal lattice vectors: (au)'
    write(6,'(a,f16.10,1x,f16.10,1x,f16.10)')'b1 ',b1(1:3)
    write(6,'(a,f16.10,1x,f16.10,1x,f16.10)')'b2 ',b2(1:3)
    write(6,'(a,f16.10,1x,f16.10,1x,f16.10)')'b3 ',b3(1:3)
    write(6,*)'Input k points:'
    do jk=1,nkvec_bwfdet
     write(6,'(i3,f16.10,1x,f16.10,1x,f16.10)')jk,kvec_bwfdet(:,jk)
    enddo
    write(6,*)'Input k point set after reduction into 1st Brillouin &
     &zone of simulation cell:'
    do jk=1,nkvec_bwfdet
     write(6,'(i3,f16.10,1x,f16.10,1x,f16.10)')jk,kvec2(:,jk)
    enddo
    write(6,*)
    write(6,*)'The k point set in pwfn.data does not map to a unique simulation'
    write(6,*)'cell k_s vector. Are you sure the npcell block in input is'
    write(6,*)'correct?'
    write(6,*)
    call errstop('BWFDET_SETUP','Quitting')
   endif
  enddo

! Print k point analysis to output
  write(6,*)'K POINT ANALYSIS'
  write(6,'(1x,&
  &''  k    kx         ky         kz       use pair edge'')')
  do ik=1,nkvec_bwfdet
   write(6,'(1x,i3,1x,f10.6,1x,f10.6,1x,f10.6,3x,l1,4x,l1,4x,l1)') &
    &ik,kvec_bwfdet(1,ik),kvec_bwfdet(2,ik),kvec_bwfdet(3,ik),lkcalc(ik),lkpair(ik),lkedge(ik)
  enddo ! k

  if(any(ltemp))then
   call errstop('BWFDET_SETUP','Supercell k_s vector not the Gamma point or &
   &half a supercell reciprocal lattice vector. Unable to make orbitals into &
   &real Bloch functions.')
  else
   write(6,*)
   write(6,*)'MAPPING ONTO UNIQUE K_S VECTOR:'
   write(6,'(3f20.8,a)')k_s(1:3),' (Cartesian a.u.)'
   if(any(abs(k_s(1:3))>tolerance))then
    temp1=one_over_twopi*dot_product(k_s,a1_bwfdet)
    temp2=one_over_twopi*dot_product(k_s,a2_bwfdet)
    temp3=one_over_twopi*dot_product(k_s,a3_bwfdet)
    write(6,'(3f20.8,a)')temp1,temp2,temp3,' (frac supercell reciprocal &
     &lattice vectors)'
   endif
  endif

! REMOVED by WDP to allow for spin-polarized calculations with
! equivalent orbitals for up and down spin
!
! Check number of orbitals and supposed number of electrons match
!
! eff_nele=0
! do ispin=1,nspin
!  do ik=1,nkvec_bwfdet
!   if(lkcalc(ik))then
!    if(lkedge(ik))then
!     eff_nele(ispin)=eff_nele(ispin)+boccband(ik,ispin)
!    else
!     eff_nele(ispin)=eff_nele(ispin)+boccband(ik,ispin)*2
!    endif
!   endif
!  enddo ! k
! enddo ! spin
! eff_nele(1)=eff_nele(1)+wf_np(1,1)-wf_nm(1,1)
! eff_nele(2)=eff_nele(2)+wf_np(1,2)-wf_nm(1,2)
! write(6,*)
! write(6,'(1x,''Computed no of orbitals                   = '',i4,1x,i4,/1x,&
!  &''Declared no of up and down spin electrons = '',i4,1x,i4)')eff_nele,nele
! if(any(eff_nele/=nele))call errstop('BWFDET_SETUP',&
!  &'Electron/orbital number mismatch.')

! Calculate kinetic energy and compare with DFT kinetic energy to make sure
! we're occupying the same set of orbitals.
!
! [Actually not so easy to do this for blips. MDT]

!  ke_check=0.d0
!  do spin=1,num_spins ! 2 for spin_polarized, 1 otherwise
!   do k=1,nkvec_bwfdet
!    if(lkcalc(k))then
!     if(lkedge(k))then
!      scale=1.d0
!     else
!      scale=2.d0
!     endif
!     do band=1,boccband(k,spin)
!      do g=1,nwvec
!       vec(:)=kvec_bwfdet(:,k)+gvecwf(:,g)
!       kg=vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
!       if(pwreal)then
!        if(spin==1)then
!         ke_check=ke_check+ckg(g,band,k)**2*kg*scale
!        else
!         ke_check=ke_check+ckg2(g,band,k)**2*kg*scale
!        endif
!       else
!        if(spin==1)then
!         ke_check=ke_check+dble(cckg(g,band,k)*conjg(cckg(g,band,k)))*kg*scale
!        else
!        ke_check=ke_check+dble(cckg2(g,band,k)*conjg(cckg2(g,band,k)))*kg*scale
!        endif
!       endif
!      enddo ! G
!     enddo ! band
!    endif
!   enddo ! k
!  enddo ! spin
!  ke_check=0.5d0*ke_check/npcells
!  if(.not.spin_polarized)ke_check=2.d0*ke_check
!
!  if(kinetic_energy/=0.d0)then
!   write(6,*)
!   write(6,*)'Kinetic energy from pwfn.data   : ',kinetic_energy
!   write(6,*)'Calculated kinetic energy       : ',ke_check
!   if(abs(kinetic_energy-ke_check)<1.d-5)then
!    write(6,*)'Tick.'
!   else
!    call errstop('BWFDET_SETUP','Kinetic energy mismatch. Probably some kind&
!     & of state occupation problem.')
!   endif
!  else
!   write(6,*)'Not possible to do KE check as info not provided in bwfn.data.'
!   write(6,*)
!  endif

 endif ! (idtask == 0
 call qmc_barrier

!norm=0.d0
!do i=2,nemax
! norm=norm+log(dble(i))
!enddo
!norm=exp(-norm/(2.d0*dble(nemax)))
 norm=1.d0
 norm=norm*orb_norm

 END SUBROUTINE bwfdet_setup


 SUBROUTINE bwfdet_main(rvec,iw,igl,ispin,rpsi,grad,lap)
!--------------------------------------------------------------------------!
! Evaluate orbitals, gradients and Laplacians for periodic system          !
! expanded in blip basis                                                   !
! DA 3.2001                                                                !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER idet,ik
 INTEGER,INTENT(in) :: iw,igl,ispin
 REAL(dp) r(3),kd
 REAL(dp),INTENT(in) :: rvec(3)
 REAL(dp),INTENT(out) :: rpsi(nemax,ndet_bwfdet),grad(3,nemax,ndet_bwfdet),lap(nemax,ndet_bwfdet)
 LOGICAL spin1

 spin1=.true. ; if(ispin==2.and.spin_polarized)spin1=.false.
 do ik=1,nkvec_bwfdet
  kd=rvec(1)*kvec_bwfdet(1,ik)+rvec(2)*kvec_bwfdet(2,ik)+rvec(3)*kvec_bwfdet(3,ik)
  zdum(ik)=norm*exp(zi*kd) ! normalization times phase factor used to multiply wavefn in blip3d
  if(igl==1)then
   lzdum(ik)=-dot_product(kvec_bwfdet(:,ik),kvec_bwfdet(:,ik))*norm*exp(zi*kd)
   gzdum(1,ik)=kvec_bwfdet(1,ik)*zi*norm*exp(zi*kd)
   gzdum(2,ik)=kvec_bwfdet(2,ik)*zi*norm*exp(zi*kd)
   gzdum(3,ik)=kvec_bwfdet(3,ik)*zi*norm*exp(zi*kd)
  endif
 enddo

! Position must be in units of crystal lattice
!r=matmul(transpose(painv),rvec)
 r=matmul(rvec,painv)

! Warning: Just as William split blip3dgamma into blip3dgamma_w and blip3dgamma_gl (no derivs and derivs) to save time
! one could do the same for blip3d
 if(spin1)then ! spin 1 (and spin 2 in spin restricted case)
  if(pwreal)then
   if(igl==0) then
    call blip3dgamma_w(rpsi,lap,grad,r,avc,nrbwf,painv,igl,boccband(1,ispin), &
     &nemax,maxband)
   else
    call blip3dgamma_gl(rpsi,lap,grad,r,avc,avclap,avcgrad1,avcgrad2,avcgrad3, &
     &nrbwf,painv,igl,boccband(1,ispin),nemax,maxband)
   endif
  else
   call blip3d(rpsi,lap,grad,r,cavc,cavclap,cavcgrad1,cavcgrad2,cavcgrad3,nrbwf,&
    &painv,iw,igl,boccband(1,ispin),nemax,nkvec_bwfdet,maxband,&
    &use_real_part(1,1,ispin))
  endif
 else          ! spin 2 (in spin unrbwfestricted case)
  if(pwreal)then
   if(igl==0) then
    call blip3dgamma_w(rpsi,lap,grad,r,avc2,nrbwf,painv,igl,boccband(1,ispin), &
     &nemax,maxband)
   else
    call blip3dgamma_gl(rpsi,lap,grad,r,avc2,avclap2,avcgrad12,avcgrad22,avcgrad32, &
     &nrbwf,painv,igl,boccband(1,ispin),nemax,maxband)
   endif
  else
   call blip3d(rpsi,lap,grad,r,cavc2,cavclap2,cavcgrad12,cavcgrad22,cavcgrad32,&
    &nrbwf,painv,iw,igl,boccband(1,ispin),nemax,nkvec_bwfdet,maxband,          &
    &use_real_part(1,1,ispin))
  endif
 endif

! William says there is no possibility of using more than 1 det in this module.
! It just sets them all equal.
 do idet=2,ndet_bwfdet
  if(iw==1)rpsi(:,idet)=rpsi(:,1)
  if(igl==1)then
   grad(:,:,idet)=grad(:,:,1) ; lap(:,idet)=lap(:,1)
  endif
 enddo ! dets

 END SUBROUTINE bwfdet_main


 SUBROUTINE blip3d(rpsi,lap,grad,r,avc,avclap,avcgrad1,avcgrad2,avcgrad3,nrbwf,&
  &bg,iw,igl,boccband,nemax,nkvec_bwfdet,maxband,use_real_part)
!---------------------------------------------------------------------------!
! This subroutine evaluates the value of a function, its gradient and its   !
! Laplacian at a vector point r, using the overlapping of blip functions.   !
! The blip grid is defined on a cubic cell, so r should always be given in  !
! units of the crystal lattice vectors.                                     !
!                                                                           !
! Input:                                                                    !
!  r(3)                 position in units of lattice vectors                !
!  avc(0:nrbwf,0:nrbwf,0:nrbwf,boccband,nkvec_bwfdet) blip coefficients     !
!  nrbwf(3)                num of divisions for each side of the box        !
!                       (defines the blip grid)                             !
!  bg(3,3)              the reciprocal lattice vectors (in a.u./tpi)        !
!                                                                           !
! DA 3.2001                                                                 !
!---------------------------------------------------------------------------!

 IMPLICIT NONE
 INTEGER,INTENT(in) :: iw,igl,nemax,nkvec_bwfdet,maxband
 INTEGER,INTENT(in) :: boccband(nkvec_bwfdet)
 INTEGER nrbwf(3),ix,iy,iz,ixp,ixpp,ixm,iyp,iypp,iym,izp,izpp,izm,idum,ib,ik
 REAL(dp) r(3),t(3),txm,tx,txp,txpp,tym,ty,typ,typp,tzm,tz,tzp,                &
  &tzpp,bg(3,3),rpsi(nemax),lap(nemax),grad(3,nemax),x,y,z,b(6)
 REAL(dp) :: dtxm=0.d0,d2txm=0.d0,dtx=0.d0,d2tx=0.d0,dtxp=0.d0,d2txp=0.d0,     &
  &dtxpp=0.d0,d2txpp=0.d0,dtym=0.d0,d2tym=0.d0,dty=0.d0,d2ty=0.d0,dtyp=0.d0,   &
  &d2typ=0.d0,dtypp=0.d0,d2typp=0.d0,dtzm=0.d0,d2tzm=0.d0,dtz=0.d0,d2tz=0.d0,  &
  &dtzp=0.d0,d2tzp=0.d0,dtzpp=0.d0,d2tzpp=0.d0
 COMPLEX(dp) :: avc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),d1(16),rdum,  &
  &ldum,gdum(3)
 COMPLEX(dp) :: avclap(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),l1(16)
 COMPLEX(dp) :: avcgrad1(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),g1(16)
 COMPLEX(dp) :: avcgrad2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),g2(16)
 COMPLEX(dp) :: avcgrad3(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband,nkvec_bwfdet),g3(16)
 LOGICAL use_real_part(maxband,nkvec_bwfdet)

 ix=int(mod(r(1)+abs(int(r(1)))+1,1.d0)*nrbwf(1))
 iy=int(mod(r(2)+abs(int(r(2)))+1,1.d0)*nrbwf(2))
 iz=int(mod(r(3)+abs(int(r(3)))+1,1.d0)*nrbwf(3))

 if(ix<0.or.iy<0.or.iz<0)call errstop('BLIP3D','Negative index found.')

! The blips are defined as the product of one-dimensional cubic splines, These
! are different from zero only on four adjacent grid points.  It follows that
! for any vector r there are only 64 overlapping three-dimensional blips, which
! are the product of all possible combinations of the three one-dimensional
! splines.

! These are the extra 3 coefficients for each dimension needed
 ix=mod(ix,nrbwf(1))
 iy=mod(iy,nrbwf(2))
 iz=mod(iz,nrbwf(3))
 ixp=mod(ix+1,nrbwf(1))
 ixpp=mod(ix+2,nrbwf(1))
 ixm=mod(ix-1+nrbwf(1),nrbwf(1))
 iyp=mod(iy+1,nrbwf(2))
 iypp=mod(iy+2,nrbwf(2))
 iym=mod(iy-1+nrbwf(2),nrbwf(2))
 izp=mod(iz+1,nrbwf(3))
 izpp=mod(iz+2,nrbwf(3))
 izm=mod(iz-1+nrbwf(3),nrbwf(3))

! Now calculate the 12 monodimensional blip functions
 t=mod(r+abs(int(r))+1,1.d0)*nrbwf
 t=mod(t,dble(nrbwf))

!Alternate method with only one definition of x
!x=t(1)-ix
!txm   = 0.25d0-0.75d0*x+0.75d0*x*x-0.25d0*x*x*x
!tx    = 1.d0           -1.5d0 *x*x+0.75d0*x*x*x
!txp   =-1.25d0+5.25d0*x-3.75d0*x*x+0.75d0*x*x*x
!txpp  =                            0.75d0*x*x*x
!y=t(1)-iy
!tym   = 0.25d0-0.75d0*y+0.75d0*y*y-0.25d0*y*y*y
!ty    = 1.d0           -1.5d0 *y*y+0.75d0*y*y*y
!typ   =-1.25d0+5.25d0*y-3.75d0*y*y+0.75d0*y*y*y
!typp  =                            0.75d0*y*y*y
!z=t(1)-iz
!tzm   = 0.25d0-0.75d0*z+0.75d0*z*z-0.25d0*z*z*z
!tz    = 1.d0           -1.5d0 *z*z+0.75d0*z*z*z
!tzp   =-1.25d0+5.25d0*z-3.75d0*z*z+0.75d0*z*z*z
!tzpp  =                            0.75d0*z*z*z
!if(igl==1)then
!  dtxm  =(       -0.75d0  +1.5d0*x   -0.75d0*x*x)*nrbwf(1)
!  dtx   =(                -3.d0*x    +2.25d0*x*x)*nrbwf(1)
!  dtxp  =(        5.25d0  -7.5d0*x   +2.25d0*x*x)*nrbwf(1)
!  dtxpp =(                            2.25d0*x*x)*nrbwf(1)
!  d2txm =(                 1.5d0     -1.5d0*x   )*nrbwf(1)**2
!  d2tx  =(                -3.d0      +4.5d0*x   )*nrbwf(1)**2
!  d2txp =(                -7.5d0     +4.5d0*x   )*nrbwf(1)**2
!  d2txpp=(                            4.5d0*x   )*nrbwf(1)**2
!  dtym  =(       -0.75d0  +1.5d0*y   -0.75d0*y*y)*nrbwf(2)
!  dty   =(                -3.d0*y    +2.25d0*y*y)*nrbwf(2)
!  dtyp  =(        5.25d0  -7.5d0*y   +2.25d0*y*y)*nrbwf(2)
!  dtypp =(                            2.25d0*y*y)*nrbwf(2)
!  d2tym =(                 1.5d0     -1.5d0*y   )*nrbwf(2)**2
!  d2ty  =(                -3.d0      +4.5d0*y   )*nrbwf(2)**2
!  d2typ =(                -7.5d0     +4.5d0*y   )*nrbwf(2)**2
!  d2typp=(                            4.5d0*y   )*nrbwf(2)**2
!  dtzm  =(       -0.75d0  +1.5d0*z   -0.75d0*z*z)*nrbwf(3)
!  dtz   =(                -3.d0*z    +2.25d0*z*z)*nrbwf(3)
!  dtzp  =(        5.25d0  -7.5d0*z   +2.25d0*z*z)*nrbwf(3)
!  dtzpp =(                            2.25d0*z*z)*nrbwf(3)
!  d2tzm =(                 1.5d0     -1.5d0*z   )*nrbwf(3)**2
!  d2tz  =(                -3.d0      +4.5d0*z   )*nrbwf(3)**2
!  d2tzp =(                -7.5d0     +4.5d0*z   )*nrbwf(3)**2
!  d2tzpp=(                            4.5d0*z   )*nrbwf(3)**2
!endif

 select case (igrad_lap)
 case(0)
   x=t(1)-ix+1.d0
   txm=2.d0-3*x+1.5d0*x*x-0.25d0*x*x*x
   if(igl==1)then
    dtxm=(-3.d0+3*x-0.75d0*x*x)*nrbwf(1)
    d2txm=(3.d0-1.5d0*x)*nrbwf(1)*nrbwf(1)
   endif
   x=t(1)-ix
   tx=1.d0-1.5d0*x*x+0.75d0*x*x*x
   if(igl==1)then
    dtx=(-3.d0*x+2.25d0*x*x)*nrbwf(1)
    d2tx=(-3.d0+4.5d0*x)*nrbwf(1)*nrbwf(1)
   endif
   x=t(1)-ix-1.d0
   txp=1.d0-1.5d0*x*x-0.75d0*x*x*x
   if(igl==1)then
    dtxp=(-3.d0*x-2.25d0*x*x)*nrbwf(1)
    d2txp=(-3.d0-4.5d0*x)*nrbwf(1)*nrbwf(1)
   endif
   x=t(1)-ix-2.d0
   txpp=2.d0+3*x+1.5d0*x*x+0.25d0*x*x*x
   if(igl==1)then
    dtxpp=(3.d0+3*x+0.75d0*x*x)*nrbwf(1)
    d2txpp=(3.d0+1.5d0*x)*nrbwf(1)*nrbwf(1)
   endif

   y=t(2)-iy+1.d0
   tym=2.d0-3*y+1.5d0*y*y-0.25d0*y*y*y
   if(igl==1)then
    dtym=(-3.d0+3*y-0.75d0*y*y)*nrbwf(2)
    d2tym=(3.d0-1.5d0*y)*nrbwf(2)*nrbwf(2)
   endif
   y=t(2)-iy
   ty=1.d0-1.5d0*y*y+0.75d0*y*y*y
   if(igl==1)then
    dty=(-3.d0*y+2.25d0*y*y)*nrbwf(2)
    d2ty=(-3.d0+4.5d0*y)*nrbwf(2)*nrbwf(2)
   endif
   y=t(2)-iy-1.d0
   typ=1.d0-1.5d0*y*y-0.75d0*y*y*y
   if(igl==1)then
    dtyp=(-3.d0*y-2.25d0*y*y)*nrbwf(2)
    d2typ=(-3.d0-4.5d0*y)*nrbwf(2)*nrbwf(2)
   endif
   y=t(2)-iy-2.d0
   typp=2.d0+3*y+1.5d0*y*y+0.25d0*y*y*y
   if(igl==1) then
    dtypp=(3.d0+3*y+0.75d0*y*y)*nrbwf(2)
    d2typp=(3.d0+1.5d0*y)*nrbwf(2)*nrbwf(2)
   endif

   z=t(3)-iz+1.d0
   tzm=2.d0-3*z+1.5d0*z*z-0.25d0*z*z*z
   if(igl==1)then
    dtzm=(-3.d0+3*z-0.75d0*z*z)*nrbwf(3)
    d2tzm=(3.d0-1.5d0*z)*nrbwf(3)*nrbwf(3)
   endif
   z=t(3)-iz
   tz=1.d0-1.5d0*z*z+0.75d0*z*z*z
   if(igl==1)then
    dtz=(-3.d0*z+2.25d0*z*z)*nrbwf(3)
    d2tz=(-3.d0+4.5d0*z)*nrbwf(3)*nrbwf(3)
   endif
   z=t(3)-iz-1.d0
   tzp=1.d0-1.5d0*z*z-0.75d0*z*z*z
   if(igl==1)then
    dtzp=(-3.d0*z-2.25d0*z*z)*nrbwf(3)
    d2tzp=(-3.d0-4.5d0*z)*nrbwf(3)*nrbwf(3)
   endif
   z=t(3)-iz-2.d0
   tzpp=2.d0+3*z+1.5d0*z*z+0.25d0*z*z*z
   if(igl==1)then
    dtzpp=(3.d0+3*z+0.75d0*z*z)*nrbwf(3)
    d2tzpp=(3.d0+1.5d0*z)*nrbwf(3)*nrbwf(3)
   endif
 case(1)
   x=t(1)-ix+1.d0
   txm=2.d0-3*x+1.5d0*x*x-0.25d0*x*x*x
   if(igl==1)then
    dtxm=(-3.d0+3*x-0.75d0*x*x)*nrbwf(1)
   endif
   x=t(1)-ix
   tx=1.d0-1.5d0*x*x+0.75d0*x*x*x
   if(igl==1)then
    dtx=(-3.d0*x+2.25d0*x*x)*nrbwf(1)
   endif
   x=t(1)-ix-1.d0
   txp=1.d0-1.5d0*x*x-0.75d0*x*x*x
   if(igl==1)then
    dtxp=(-3.d0*x-2.25d0*x*x)*nrbwf(1)
   endif
   x=t(1)-ix-2.d0
   txpp=2.d0+3*x+1.5d0*x*x+0.25d0*x*x*x
   if(igl==1)then
    dtxpp=(3.d0+3*x+0.75d0*x*x)*nrbwf(1)
   endif

   y=t(2)-iy+1.d0
   tym=2.d0-3*y+1.5d0*y*y-0.25d0*y*y*y
   if(igl==1)then
    dtym=(-3.d0+3*y-0.75d0*y*y)*nrbwf(2)
   endif
   y=t(2)-iy
   ty=1.d0-1.5d0*y*y+0.75d0*y*y*y
   if(igl==1)then
    dty=(-3.d0*y+2.25d0*y*y)*nrbwf(2)
   endif
   y=t(2)-iy-1.d0
   typ=1.d0-1.5d0*y*y-0.75d0*y*y*y
   if(igl==1)then
    dtyp=(-3.d0*y-2.25d0*y*y)*nrbwf(2)
   endif
   y=t(2)-iy-2.d0
   typp=2.d0+3*y+1.5d0*y*y+0.25d0*y*y*y
   if(igl==1) then
    dtypp=(3.d0+3*y+0.75d0*y*y)*nrbwf(2)
   endif

   z=t(3)-iz+1.d0
   tzm=2.d0-3*z+1.5d0*z*z-0.25d0*z*z*z
   if(igl==1)then
    dtzm=(-3.d0+3*z-0.75d0*z*z)*nrbwf(3)
   endif
   z=t(3)-iz
   tz=1.d0-1.5d0*z*z+0.75d0*z*z*z
   if(igl==1)then
    dtz=(-3.d0*z+2.25d0*z*z)*nrbwf(3)
   endif
   z=t(3)-iz-1.d0
   tzp=1.d0-1.5d0*z*z-0.75d0*z*z*z
   if(igl==1)then
    dtzp=(-3.d0*z-2.25d0*z*z)*nrbwf(3)
   endif
   z=t(3)-iz-2.d0
   tzpp=2.d0+3*z+1.5d0*z*z+0.25d0*z*z*z
   if(igl==1)then
    dtzpp=(3.d0+3*z+0.75d0*z*z)*nrbwf(3)
   endif
 case(2)
   x=t(1)-ix+1.d0
   txm=2.d0-3*x+1.5d0*x*x-0.25d0*x*x*x
   x=t(1)-ix
   tx=1.d0-1.5d0*x*x+0.75d0*x*x*x
   x=t(1)-ix-1.d0
   txp=1.d0-1.5d0*x*x-0.75d0*x*x*x
   x=t(1)-ix-2.d0
   txpp=2.d0+3*x+1.5d0*x*x+0.25d0*x*x*x

   y=t(2)-iy+1.d0
   tym=2.d0-3*y+1.5d0*y*y-0.25d0*y*y*y
   y=t(2)-iy
   ty=1.d0-1.5d0*y*y+0.75d0*y*y*y
   y=t(2)-iy-1.d0
   typ=1.d0-1.5d0*y*y-0.75d0*y*y*y
   y=t(2)-iy-2.d0
   typp=2.d0+3*y+1.5d0*y*y+0.25d0*y*y*y

   z=t(3)-iz+1.d0
   tzm=2.d0-3*z+1.5d0*z*z-0.25d0*z*z*z
   z=t(3)-iz
   tz=1.d0-1.5d0*z*z+0.75d0*z*z*z
   z=t(3)-iz-1.d0
   tzp=1.d0-1.5d0*z*z-0.75d0*z*z*z
   z=t(3)-iz-2.d0
   tzpp=2.d0+3*z+1.5d0*z*z+0.25d0*z*z*z
 end select

! bg are the primitive? cell reciprocal lattice vectors
 b(1)=(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2)
 b(2)=(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2)
 b(3)=(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2)
 b(4)=2*(bg(1,1)*bg(1,2)+bg(2,1)*bg(2,2)+bg(3,1)*bg(3,2))
 b(5)=2*(bg(1,2)*bg(1,3)+bg(2,2)*bg(2,3)+bg(3,2)*bg(3,3))
 b(6)=2*(bg(1,3)*bg(1,1)+bg(2,3)*bg(2,1)+bg(3,3)*bg(3,1))

 idum=0
 do ik=1,nkvec_bwfdet
  if(lkcalc(ik))then
   do ib=1,boccband(ik)
    d1(1)=avc(ix,iy,iz,ib,ik)*tz+avc(ix,iy,izp,ib,ik)*tzp+                     &
         avc(ix,iy,izpp,ib,ik)*tzpp+avc(ix,iy,izm,ib,ik)*tzm
    d1(2)=avc(ix,iyp,iz,ib,ik)*tz+avc(ix,iyp,izp,ib,ik)*tzp+                   &
         avc(ix,iyp,izpp,ib,ik)*tzpp+avc(ix,iyp,izm,ib,ik)*tzm
    d1(3)=avc(ix,iypp,iz,ib,ik)*tz+avc(ix,iypp,izp,ib,ik)*tzp+                 &
         avc(ix,iypp,izpp,ib,ik)*tzpp+avc(ix,iypp,izm,ib,ik)*tzm
    d1(4)=avc(ix,iym,iz,ib,ik)*tz+avc(ix,iym,izp,ib,ik)*tzp+                   &
         avc(ix,iym,izpp,ib,ik)*tzpp+avc(ix,iym,izm,ib,ik)*tzm

    d1(5)=avc(ixp,iy,iz,ib,ik)*tz+avc(ixp,iy,izp,ib,ik)*tzp                    &
         +avc(ixp,iy,izpp,ib,ik)*tzpp+avc(ixp,iy,izm,ib,ik)*tzm
    d1(6)=avc(ixp,iyp,iz,ib,ik)*tz+avc(ixp,iyp,izp,ib,ik)*tzp                  &
         +avc(ixp,iyp,izpp,ib,ik)*tzpp+avc(ixp,iyp,izm,ib,ik)*tzm
    d1(7)=avc(ixp,iypp,iz,ib,ik)*tz+avc(ixp,iypp,izp,ib,ik)*tzp                &
         +avc(ixp,iypp,izpp,ib,ik)*tzpp+avc(ixp,iypp,izm,ib,ik)*tzm
    d1(8)=avc(ixp,iym,iz,ib,ik)*tz+avc(ixp,iym,izp,ib,ik)*tzp                  &
         +avc(ixp,iym,izpp,ib,ik)*tzpp+avc(ixp,iym,izm,ib,ik)*tzm

    d1(9)=avc(ixpp,iy,iz,ib,ik)*tz+avc(ixpp,iy,izp,ib,ik)*tzp                  &
         +avc(ixpp,iy,izpp,ib,ik)*tzpp+avc(ixpp,iy,izm,ib,ik)*tzm
    d1(10)=avc(ixpp,iyp,iz,ib,ik)*tz+avc(ixpp,iyp,izp,ib,ik)*tzp               &
         +avc(ixpp,iyp,izpp,ib,ik)*tzpp+avc(ixpp,iyp,izm,ib,ik)*tzm
    d1(11)=avc(ixpp,iypp,iz,ib,ik)*tz+avc(ixpp,iypp,izp,ib,ik)*tzp             &
         +avc(ixpp,iypp,izpp,ib,ik)*tzpp+avc(ixpp,iypp,izm,ib,ik)*tzm
    d1(12)=avc(ixpp,iym,iz,ib,ik)*tz+avc(ixpp,iym,izp,ib,ik)*tzp               &
         +avc(ixpp,iym,izpp,ib,ik)*tzpp+avc(ixpp,iym,izm,ib,ik)*tzm

    d1(13)=avc(ixm,iy,iz,ib,ik)*tz+avc(ixm,iy,izp,ib,ik)*tzp                   &
         +avc(ixm,iy,izpp,ib,ik)*tzpp+avc(ixm,iy,izm,ib,ik)*tzm
    d1(14)=avc(ixm,iyp,iz,ib,ik)*tz+avc(ixm,iyp,izp,ib,ik)*tzp                 &
         +avc(ixm,iyp,izpp,ib,ik)*tzpp+avc(ixm,iyp,izm,ib,ik)*tzm
    d1(15)=avc(ixm,iypp,iz,ib,ik)*tz+avc(ixm,iypp,izp,ib,ik)*tzp               &
         +avc(ixm,iypp,izpp,ib,ik)*tzpp+avc(ixm,iypp,izm,ib,ik)*tzm
    d1(16)=avc(ixm,iym,iz,ib,ik)*tz+avc(ixm,iym,izp,ib,ik)*tzp                 &
         +avc(ixm,iym,izpp,ib,ik)*tzpp+avc(ixm,iym,izm,ib,ik)*tzm

! The function
    rdum=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*typ      &
         +d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*   &
         tym)*txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm
    if(iw==1) then
     if(lkedge(ik))then
      if(use_real_part(ib,ik))then
       rpsi(idum+ib)=dble(zdum(ik)*rdum) ! multiply wavefn by phase factor
      else
       rpsi(idum+ib)=aimag(zdum(ik)*rdum)
      endif
     else
      rpsi(idum+1)=dble(zdum(ik)*rdum)
      rpsi(idum+2)=aimag(zdum(ik)*rdum)
     endif
    endif

    if(igl==1)then
     select case (igrad_lap)
     case (0)
! The Laplacian: first term involving Theta''(x)
     ldum=((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*d2tx+(d1(5)*ty+d1(6)*typ+ &
          d1(7)*typp+d1(8)*tym)*d2txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)* &
          tym)*d2txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*d2txm)*b(1)

! The Laplacian: term involving Theta'(x)Theta'(y)
     ldum=ldum+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*dtx+(d1(5)*dty+  &
          d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*dtxp+(d1(9)*dty+d1(10)*dtyp+      &
          d1(11)*dtypp+d1(12)*dtym)*dtxpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*      &
          dtypp+d1(16)*dtym)*dtxm)*b(4)

! The gradient, first term involving Theta'(x)
     gdum(1)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+(d1(5)*ty+d1(6)*    &
          typ+d1(7)*typp+d1(8)*tym)*dtxp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+     &
          d1(12)*tym)*dtxpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm

! Second term of the laplacian involving Theta''(y)
     ldum=ldum+((d1(1)*d2ty+d1(2)*d2typ+d1(3)*d2typp+d1(4)*d2tym)*tx+          &
          (d1(5)*d2ty+d1(6)*d2typ+d1(7)*d2typp+d1(8)*d2tym)*txp+(d1(9)*d2ty+   &
          d1(10)*d2typ+d1(11)*d2typp+d1(12)*d2tym)*txpp+(d1(13)*d2ty+d1(14)*   &
          d2typ+d1(15)*d2typp+d1(16)*d2tym)*txm)*b(2)

! Second term of the gradient involving Theta'(y)
     gdum(2)=(d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+(d1(5)*dty+      &
          d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+(d1(9)*dty+d1(10)*dtyp+       &
          d1(11)*dtypp+d1(12)*dtym)*txpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+ &
          d1(16)*dtym)*txm

! And now the third term of the laplacian involving Theta''(z)
     d1(1)=avc(ix,iy,iz,ib,ik)*d2tz+avc(ix,iy,izp,ib,ik)*d2tzp+                &
          avc(ix,iy,izpp,ib,ik)*d2tzpp+avc(ix,iy,izm,ib,ik)*d2tzm
     d1(2)=avc(ix,iyp,iz,ib,ik)*d2tz+avc(ix,iyp,izp,ib,ik)*d2tzp+              &
          avc(ix,iyp,izpp,ib,ik)*d2tzpp+avc(ix,iyp,izm,ib,ik)*d2tzm
     d1(3)=avc(ix,iypp,iz,ib,ik)*d2tz+avc(ix,iypp,izp,ib,ik)*d2tzp+            &
          avc(ix,iypp,izpp,ib,ik)*d2tzpp+avc(ix,iypp,izm,ib,ik)*d2tzm
     d1(4)=avc(ix,iym,iz,ib,ik)*d2tz+avc(ix,iym,izp,ib,ik)*d2tzp+              &
          avc(ix,iym,izpp,ib,ik)*d2tzpp+avc(ix,iym,izm,ib,ik)*d2tzm

     d1(5)=avc(ixp,iy,iz,ib,ik)*d2tz+avc(ixp,iy,izp,ib,ik)*d2tzp+              &
          avc(ixp,iy,izpp,ib,ik)*d2tzpp+avc(ixp,iy,izm,ib,ik)*d2tzm
     d1(6)=avc(ixp,iyp,iz,ib,ik)*d2tz+avc(ixp,iyp,izp,ib,ik)*d2tzp+            &
          avc(ixp,iyp,izpp,ib,ik)*d2tzpp+avc(ixp,iyp,izm,ib,ik)*d2tzm
     d1(7)=avc(ixp,iypp,iz,ib,ik)*d2tz+avc(ixp,iypp,izp,ib,ik)*d2tzp+          &
          avc(ixp,iypp,izpp,ib,ik)*d2tzpp+avc(ixp, iypp, izm,ib,ik)*d2tzm
     d1(8)=avc(ixp,iym,iz,ib,ik)*d2tz+avc(ixp,iym,izp,ib,ik)*d2tzp+            &
          avc(ixp,iym,izpp,ib,ik)*d2tzpp+avc(ixp,iym,izm,ib,ik)*d2tzm

     d1(9)=avc(ixpp,iy,iz,ib,ik)*d2tz+avc(ixpp,iy,izp,ib,ik)*d2tzp+            &
          avc(ixpp,iy,izpp,ib,ik)*d2tzpp+avc(ixpp,iy,izm,ib,ik)*d2tzm
     d1(10)=avc(ixpp,iyp,iz,ib,ik)*d2tz+avc(ixpp,iyp,izp,ib,ik)*d2tzp+         &
          avc(ixpp,iyp,izpp,ib,ik)*d2tzpp+avc(ixpp,iyp,izm,ib,ik)*d2tzm
     d1(11)=avc(ixpp,iypp,iz,ib,ik)*d2tz+avc(ixpp,iypp,izp,ib,ik)*d2tzp+       &
          avc(ixpp,iypp,izpp,ib,ik)*d2tzpp+avc(ixpp,iypp,izm,ib,ik)*d2tzm
     d1(12)=avc(ixpp,iym,iz,ib,ik)*d2tz+avc(ixpp,iym,izp,ib,ik)*d2tzp+         &
          avc(ixpp,iym,izpp,ib,ik)*d2tzpp+avc(ixpp,iym,izm,ib,ik)*d2tzm

     d1(13)=avc(ixm,iy,iz,ib,ik)*d2tz+avc(ixm,iy,izp,ib,ik)*d2tzp+             &
          avc(ixm,iy,izpp,ib,ik)*d2tzpp+avc(ixm,iy,izm,ib,ik)*d2tzm
     d1(14)=avc(ixm,iyp,iz,ib,ik)*d2tz+avc(ixm,iyp,izp,ib,ik)*d2tzp+           &
          avc(ixm,iyp,izpp,ib,ik)*d2tzpp+avc(ixm,iyp,izm,ib,ik)*d2tzm
     d1(15)=avc(ixm,iypp,iz,ib,ik)*d2tz+avc(ixm,iypp,izp,ib,ik)*d2tzp+         &
          avc(ixm,iypp,izpp,ib,ik)*d2tzpp+avc(ixm,iypp,izm,ib,ik)*d2tzm
     d1(16)=avc(ixm,iym,iz,ib,ik)*d2tz+avc(ixm,iym,izp,ib,ik)*d2tzp+           &
          avc(ixm,iym,izpp,ib,ik)*d2tzpp+avc(ixm,iym,izm,ib,ik)*d2tzm

! Theta''(z)
     ldum=ldum+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*  &
          typ+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+      &
          d1(12)* tym)*txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm) &
          *b(3)

! And the third term of the gradient involving Theta'(z)
     d1(1)=avc(ix,iy,iz,ib,ik)*dtz+avc(ix,iy,izp,ib,ik)*dtzp+                  &
          avc(ix,iy,izpp,ib,ik)*dtzpp+avc(ix, iy, izm,ib,ik)*dtzm
     d1(2)=avc(ix,iyp,iz,ib,ik)*dtz+avc(ix,iyp,izp,ib,ik)*dtzp+                &
          avc(ix,iyp,izpp,ib,ik)*dtzpp+avc(ix,iyp,izm,ib,ik)*dtzm
     d1(3)=avc(ix,iypp,iz,ib,ik)*dtz+avc(ix,iypp,izp,ib,ik)*dtzp+              &
          avc(ix,iypp,izpp,ib,ik)*dtzpp+avc(ix,iypp,izm,ib,ik)*dtzm
     d1(4)=avc(ix,iym,iz,ib,ik)*dtz+avc(ix,iym,izp,ib,ik)*dtzp+                &
          avc(ix,iym,izpp,ib,ik)*dtzpp+avc(ix,iym,izm,ib,ik)*dtzm

     d1(5)=avc(ixp,iy,iz,ib,ik)*dtz+avc(ixp,iy,izp,ib,ik)*dtzp+                &
          avc(ixp,iy,izpp,ib,ik)*dtzpp+avc(ixp,iy,izm,ib,ik)*dtzm
     d1(6)=avc(ixp,iyp,iz,ib,ik)*dtz+avc(ixp,iyp,izp,ib,ik)*dtzp+              &
          avc(ixp,iyp,izpp,ib,ik)*dtzpp+avc(ixp,iyp,izm,ib,ik)*dtzm
     d1(7)=avc(ixp,iypp,iz,ib,ik)*dtz+avc(ixp,iypp,izp,ib,ik)*dtzp+            &
          avc(ixp,iypp,izpp,ib,ik)*dtzpp+avc(ixp,iypp,izm,ib,ik)*dtzm
     d1(8)=avc(ixp,iym,iz,ib,ik)*dtz+avc(ixp,iym,izp,ib,ik)*dtzp+              &
          avc(ixp,iym,izpp,ib,ik)*dtzpp+avc(ixp,iym,izm,ib,ik)*dtzm

     d1(9)=avc(ixpp,iy,iz,ib,ik)*dtz+avc(ixpp,iy,izp,ib,ik)*dtzp+              &
          avc(ixpp,iy,izpp,ib,ik)*dtzpp+avc(ixpp,iy,izm,ib,ik)*dtzm
     d1(10)=avc(ixpp,iyp,iz,ib,ik)*dtz+avc(ixpp,iyp,izp,ib,ik)*dtzp+           &
          avc(ixpp,iyp,izpp,ib,ik)*dtzpp+avc(ixpp,iyp,izm,ib,ik)*dtzm
     d1(11)=avc(ixpp,iypp,iz,ib,ik)*dtz+avc(ixpp,iypp,izp,ib,ik)*dtzp+         &
          avc(ixpp,iypp,izpp,ib,ik)*dtzpp+avc(ixpp,iypp,izm,ib,ik)*dtzm
     d1(12)=avc(ixpp,iym,iz,ib,ik)*dtz+avc(ixpp,iym,izp,ib,ik)*dtzp+           &
          avc(ixpp,iym,izpp,ib,ik)*dtzpp+avc(ixpp,iym,izm,ib,ik)*dtzm

     d1(13)=avc(ixm,iy,iz,ib,ik)*dtz+avc(ixm,iy,izp,ib,ik)*dtzp+               &
          avc(ixm,iy,izpp,ib,ik)*dtzpp+avc(ixm,iy,izm,ib,ik)*dtzm
     d1(14)=avc(ixm,iyp,iz,ib,ik)*dtz+avc(ixm,iyp,izp,ib,ik)*dtzp+             &
          avc(ixm,iyp,izpp,ib,ik)*dtzpp+avc(ixm,iyp,izm,ib,ik)*dtzm
     d1(15)=avc(ixm,iypp,iz,ib,ik)*dtz+avc(ixm,iypp,izp,ib,ik)*dtzp+           &
          avc(ixm,iypp,izpp,ib,ik)*dtzpp+avc(ixm,iypp,izm,ib,ik)*dtzm
     d1(16)=avc(ixm,iym,iz,ib,ik)*dtz+avc(ixm,iym,izp,ib,ik)*dtzp+             &
          avc(ixm,iym,izpp,ib,ik)*dtzpp+avc(ixm,iym,izm,ib,ik)*dtzm

! The Laplacian: term involving Theta'(x)Theta'(z)
     ldum=ldum+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+(d1(5)*ty+       &
          d1(6)*typ+d1(7)*typp+d1(8)*tym)*dtxp+(d1(9)*ty+d1(10)*typ+d1(11)*    &
          typp+d1(12)*tym)*dtxpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym) &
          *dtxm)*b(6)

! the Laplacian: term involving Theta'(y)Theta'(z)
     ldum=ldum+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+(d1(5)        &
          *dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+(d1(9)*dty+d1(10)*dtyp+  &
          d1(11)*dtypp+d1(12)*dtym)*txpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+ &
          d1(16)*dtym)*txm)*b(5)

     gdum(3)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*     &
          typ+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+      &
          d1(12)*tym)*txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

! go to the Cartesian grid
     gdum=matmul(bg,gdum)

     if(lkedge(ik))then
      if(use_real_part(ib,ik))then
       lap(idum+ib)=dble(lzdum(ik)*rdum+                                       &
            2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+   &
            zdum(ik)*ldum)
       grad(1,idum+ib)=dble(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
       grad(2,idum+ib)=dble(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
       grad(3,idum+ib)=dble(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
      else
         lap(idum+ib)=aimag(lzdum(ik)*rdum+                                    &
              2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+ &
              zdum(ik)*ldum)
         grad(1,idum+ib)=aimag(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
         grad(2,idum+ib)=aimag(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
         grad(3,idum+ib)=aimag(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
      endif
     else
      lap(idum+1)=dble(lzdum(ik)*rdum+                                         &
           2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+    &
           zdum(ik)*ldum)
      lap(idum+2)=aimag(lzdum(ik)*rdum+                                        &
           2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+    &
           zdum(ik)*ldum)
      grad(1,idum+1)=dble(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
      grad(2,idum+1)=dble(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
      grad(3,idum+1)=dble(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
      grad(1,idum+2)=aimag(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
      grad(2,idum+2)=aimag(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
      grad(3,idum+2)=aimag(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
     endif
    case(1)
    l1(1)=avclap(ix,iy,iz,ib,ik)*tz+avclap(ix,iy,izp,ib,ik)*tzp+        &
         avclap(ix,iy,izpp,ib,ik)*tzpp+avclap(ix,iy,izm,ib,ik)*tzm
    l1(2)=avclap(ix,iyp,iz,ib,ik)*tz+avclap(ix,iyp,izp,ib,ik)*tzp+      &
         avclap(ix,iyp,izpp,ib,ik)*tzpp+avclap(ix,iyp,izm,ib,ik)*tzm
    l1(3)=avclap(ix,iypp,iz,ib,ik)*tz+avclap(ix,iypp,izp,ib,ik)*tzp+    &
         avclap(ix,iypp,izpp,ib,ik)*tzpp+avclap(ix,iypp,izm,ib,ik)*tzm
    l1(4)=avclap(ix,iym,iz,ib,ik)*tz+avclap(ix,iym,izp,ib,ik)*tzp+      &
         avclap(ix,iym,izpp,ib,ik)*tzpp+avclap(ix,iym,izm,ib,ik)*tzm

    l1(5)=avclap(ixp,iy,iz,ib,ik)*tz+avclap(ixp,iy,izp,ib,ik)*tzp       &
         +avclap(ixp,iy,izpp,ib,ik)*tzpp+avclap(ixp,iy,izm,ib,ik)*tzm
    l1(6)=avclap(ixp,iyp,iz,ib,ik)*tz+avclap(ixp,iyp,izp,ib,ik)*tzp     &
         +avclap(ixp,iyp,izpp,ib,ik)*tzpp+avclap(ixp,iyp,izm,ib,ik)*tzm
    l1(7)=avclap(ixp,iypp,iz,ib,ik)*tz+avclap(ixp,iypp,izp,ib,ik)*tzp   &
         +avclap(ixp,iypp,izpp,ib,ik)*tzpp+avclap(ixp,iypp,izm,ib,ik)*tzm
    l1(8)=avclap(ixp,iym,iz,ib,ik)*tz+avclap(ixp,iym,izp,ib,ik)*tzp     &
         +avclap(ixp,iym,izpp,ib,ik)*tzpp+avclap(ixp,iym,izm,ib,ik)*tzm

    l1(9)=avclap(ixpp,iy,iz,ib,ik)*tz+avclap(ixpp,iy,izp,ib,ik)*tzp     &
         +avclap(ixpp,iy,izpp,ib,ik)*tzpp+avclap(ixpp,iy,izm,ib,ik)*tzm
    l1(10)=avclap(ixpp,iyp,iz,ib,ik)*tz+avclap(ixpp,iyp,izp,ib,ik)*tzp  &
         +avclap(ixpp,iyp,izpp,ib,ik)*tzpp+avclap(ixpp,iyp,izm,ib,ik)*tzm
    l1(11)=avclap(ixpp,iypp,iz,ib,ik)*tz+avclap(ixpp,iypp,izp,ib,ik)*tzp&
         +avclap(ixpp,iypp,izpp,ib,ik)*tzpp+avclap(ixpp,iypp,izm,ib,ik)*tzm
    l1(12)=avclap(ixpp,iym,iz,ib,ik)*tz+avclap(ixpp,iym,izp,ib,ik)*tzp  &
         +avclap(ixpp,iym,izpp,ib,ik)*tzpp+avclap(ixpp,iym,izm,ib,ik)*tzm

    l1(13)=avclap(ixm,iy,iz,ib,ik)*tz+avclap(ixm,iy,izp,ib,ik)*tzp      &
         +avclap(ixm,iy,izpp,ib,ik)*tzpp+avclap(ixm,iy,izm,ib,ik)*tzm
    l1(14)=avclap(ixm,iyp,iz,ib,ik)*tz+avclap(ixm,iyp,izp,ib,ik)*tzp    &
         +avclap(ixm,iyp,izpp,ib,ik)*tzpp+avclap(ixm,iyp,izm,ib,ik)*tzm
    l1(15)=avclap(ixm,iypp,iz,ib,ik)*tz+avclap(ixm,iypp,izp,ib,ik)*tzp  &
         +avclap(ixm,iypp,izpp,ib,ik)*tzpp+avclap(ixm,iypp,izm,ib,ik)*tzm
    l1(16)=avclap(ixm,iym,iz,ib,ik)*tz+avclap(ixm,iym,izp,ib,ik)*tzp    &
         +avclap(ixm,iym,izpp,ib,ik)*tzpp+avclap(ixm,iym,izm,ib,ik)*tzm

    ldum=(l1(1)*ty+l1(2)*typ+l1(3)*typp+l1(4)*tym)*tx+(l1(5)*ty+l1(6)*typ      &
         +l1(7)*typp+l1(8)*tym)*txp+(l1(9)*ty+l1(10)*typ+l1(11)*typp+l1(12)*   &
         tym)*txpp+(l1(13)*ty+l1(14)*typ+l1(15)*typp+l1(16)*tym)*txm

! The gradient, first term involving Theta'(x)
     gdum(1)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+(d1(5)*ty+d1(6)*    &
          typ+d1(7)*typp+d1(8)*tym)*dtxp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+     &
          d1(12)*tym)*dtxpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm

! Second term of the gradient involving Theta'(y)
     gdum(2)=(d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+(d1(5)*dty+      &
          d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+(d1(9)*dty+d1(10)*dtyp+       &
          d1(11)*dtypp+d1(12)*dtym)*txpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+ &
          d1(16)*dtym)*txm

! And the third term of the gradient involving Theta'(z)
     d1(1)=avc(ix,iy,iz,ib,ik)*dtz+avc(ix,iy,izp,ib,ik)*dtzp+                  &
          avc(ix,iy,izpp,ib,ik)*dtzpp+avc(ix, iy, izm,ib,ik)*dtzm
     d1(2)=avc(ix,iyp,iz,ib,ik)*dtz+avc(ix,iyp,izp,ib,ik)*dtzp+                &
          avc(ix,iyp,izpp,ib,ik)*dtzpp+avc(ix,iyp,izm,ib,ik)*dtzm
     d1(3)=avc(ix,iypp,iz,ib,ik)*dtz+avc(ix,iypp,izp,ib,ik)*dtzp+              &
          avc(ix,iypp,izpp,ib,ik)*dtzpp+avc(ix,iypp,izm,ib,ik)*dtzm
     d1(4)=avc(ix,iym,iz,ib,ik)*dtz+avc(ix,iym,izp,ib,ik)*dtzp+                &
          avc(ix,iym,izpp,ib,ik)*dtzpp+avc(ix,iym,izm,ib,ik)*dtzm

     d1(5)=avc(ixp,iy,iz,ib,ik)*dtz+avc(ixp,iy,izp,ib,ik)*dtzp+                &
          avc(ixp,iy,izpp,ib,ik)*dtzpp+avc(ixp,iy,izm,ib,ik)*dtzm
     d1(6)=avc(ixp,iyp,iz,ib,ik)*dtz+avc(ixp,iyp,izp,ib,ik)*dtzp+              &
          avc(ixp,iyp,izpp,ib,ik)*dtzpp+avc(ixp,iyp,izm,ib,ik)*dtzm
     d1(7)=avc(ixp,iypp,iz,ib,ik)*dtz+avc(ixp,iypp,izp,ib,ik)*dtzp+            &
          avc(ixp,iypp,izpp,ib,ik)*dtzpp+avc(ixp,iypp,izm,ib,ik)*dtzm
     d1(8)=avc(ixp,iym,iz,ib,ik)*dtz+avc(ixp,iym,izp,ib,ik)*dtzp+              &
          avc(ixp,iym,izpp,ib,ik)*dtzpp+avc(ixp,iym,izm,ib,ik)*dtzm

     d1(9)=avc(ixpp,iy,iz,ib,ik)*dtz+avc(ixpp,iy,izp,ib,ik)*dtzp+              &
          avc(ixpp,iy,izpp,ib,ik)*dtzpp+avc(ixpp,iy,izm,ib,ik)*dtzm
     d1(10)=avc(ixpp,iyp,iz,ib,ik)*dtz+avc(ixpp,iyp,izp,ib,ik)*dtzp+           &
          avc(ixpp,iyp,izpp,ib,ik)*dtzpp+avc(ixpp,iyp,izm,ib,ik)*dtzm
     d1(11)=avc(ixpp,iypp,iz,ib,ik)*dtz+avc(ixpp,iypp,izp,ib,ik)*dtzp+         &
          avc(ixpp,iypp,izpp,ib,ik)*dtzpp+avc(ixpp,iypp,izm,ib,ik)*dtzm
     d1(12)=avc(ixpp,iym,iz,ib,ik)*dtz+avc(ixpp,iym,izp,ib,ik)*dtzp+           &
          avc(ixpp,iym,izpp,ib,ik)*dtzpp+avc(ixpp,iym,izm,ib,ik)*dtzm

     d1(13)=avc(ixm,iy,iz,ib,ik)*dtz+avc(ixm,iy,izp,ib,ik)*dtzp+               &
          avc(ixm,iy,izpp,ib,ik)*dtzpp+avc(ixm,iy,izm,ib,ik)*dtzm
     d1(14)=avc(ixm,iyp,iz,ib,ik)*dtz+avc(ixm,iyp,izp,ib,ik)*dtzp+             &
          avc(ixm,iyp,izpp,ib,ik)*dtzpp+avc(ixm,iyp,izm,ib,ik)*dtzm
     d1(15)=avc(ixm,iypp,iz,ib,ik)*dtz+avc(ixm,iypp,izp,ib,ik)*dtzp+           &
          avc(ixm,iypp,izpp,ib,ik)*dtzpp+avc(ixm,iypp,izm,ib,ik)*dtzm
     d1(16)=avc(ixm,iym,iz,ib,ik)*dtz+avc(ixm,iym,izp,ib,ik)*dtzp+             &
          avc(ixm,iym,izpp,ib,ik)*dtzpp+avc(ixm,iym,izm,ib,ik)*dtzm

     gdum(3)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*     &
          typ+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+      &
          d1(12)*tym)*txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

! go to the Cartesian grid
     gdum=matmul(bg,gdum)

     if(lkedge(ik))then
      if(use_real_part(ib,ik))then
       lap(idum+ib)=dble(lzdum(ik)*rdum+                                       &
            2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+   &
                          zdum(ik)*ldum)
       grad(1,idum+ib)=dble(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
       grad(2,idum+ib)=dble(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
       grad(3,idum+ib)=dble(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
      else
         lap(idum+ib)=aimag(lzdum(ik)*rdum+                                    &
              2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+ &
                             zdum(ik)*ldum)
         grad(1,idum+ib)=aimag(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
         grad(2,idum+ib)=aimag(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
         grad(3,idum+ib)=aimag(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
      endif
     else
      lap(idum+1)=dble(lzdum(ik)*rdum+                                         &
            2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+   &
                        zdum(ik)*ldum)
      lap(idum+2)=aimag(lzdum(ik)*rdum+                                        &
             2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+  &
                         zdum(ik)*ldum)
      grad(1,idum+1)=dble(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
      grad(2,idum+1)=dble(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
      grad(3,idum+1)=dble(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
      grad(1,idum+2)=aimag(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
      grad(2,idum+2)=aimag(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
      grad(3,idum+2)=aimag(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
     endif
    case(2)
    l1(1)=avclap(ix,iy,iz,ib,ik)*tz+avclap(ix,iy,izp,ib,ik)*tzp+        &
         avclap(ix,iy,izpp,ib,ik)*tzpp+avclap(ix,iy,izm,ib,ik)*tzm
    l1(2)=avclap(ix,iyp,iz,ib,ik)*tz+avclap(ix,iyp,izp,ib,ik)*tzp+      &
         avclap(ix,iyp,izpp,ib,ik)*tzpp+avclap(ix,iyp,izm,ib,ik)*tzm
    l1(3)=avclap(ix,iypp,iz,ib,ik)*tz+avclap(ix,iypp,izp,ib,ik)*tzp+    &
         avclap(ix,iypp,izpp,ib,ik)*tzpp+avclap(ix,iypp,izm,ib,ik)*tzm
    l1(4)=avclap(ix,iym,iz,ib,ik)*tz+avclap(ix,iym,izp,ib,ik)*tzp+      &
         avclap(ix,iym,izpp,ib,ik)*tzpp+avclap(ix,iym,izm,ib,ik)*tzm

    l1(5)=avclap(ixp,iy,iz,ib,ik)*tz+avclap(ixp,iy,izp,ib,ik)*tzp       &
         +avclap(ixp,iy,izpp,ib,ik)*tzpp+avclap(ixp,iy,izm,ib,ik)*tzm
    l1(6)=avclap(ixp,iyp,iz,ib,ik)*tz+avclap(ixp,iyp,izp,ib,ik)*tzp     &
         +avclap(ixp,iyp,izpp,ib,ik)*tzpp+avclap(ixp,iyp,izm,ib,ik)*tzm
    l1(7)=avclap(ixp,iypp,iz,ib,ik)*tz+avclap(ixp,iypp,izp,ib,ik)*tzp   &
         +avclap(ixp,iypp,izpp,ib,ik)*tzpp+avclap(ixp,iypp,izm,ib,ik)*tzm
    l1(8)=avclap(ixp,iym,iz,ib,ik)*tz+avclap(ixp,iym,izp,ib,ik)*tzp     &
         +avclap(ixp,iym,izpp,ib,ik)*tzpp+avclap(ixp,iym,izm,ib,ik)*tzm

    l1(9)=avclap(ixpp,iy,iz,ib,ik)*tz+avclap(ixpp,iy,izp,ib,ik)*tzp     &
         +avclap(ixpp,iy,izpp,ib,ik)*tzpp+avclap(ixpp,iy,izm,ib,ik)*tzm
    l1(10)=avclap(ixpp,iyp,iz,ib,ik)*tz+avclap(ixpp,iyp,izp,ib,ik)*tzp  &
         +avclap(ixpp,iyp,izpp,ib,ik)*tzpp+avclap(ixpp,iyp,izm,ib,ik)*tzm
    l1(11)=avclap(ixpp,iypp,iz,ib,ik)*tz+avclap(ixpp,iypp,izp,ib,ik)*tzp&
         +avclap(ixpp,iypp,izpp,ib,ik)*tzpp+avclap(ixpp,iypp,izm,ib,ik)*tzm
    l1(12)=avclap(ixpp,iym,iz,ib,ik)*tz+avclap(ixpp,iym,izp,ib,ik)*tzp  &
         +avclap(ixpp,iym,izpp,ib,ik)*tzpp+avclap(ixpp,iym,izm,ib,ik)*tzm

    l1(13)=avclap(ixm,iy,iz,ib,ik)*tz+avclap(ixm,iy,izp,ib,ik)*tzp      &
         +avclap(ixm,iy,izpp,ib,ik)*tzpp+avclap(ixm,iy,izm,ib,ik)*tzm
    l1(14)=avclap(ixm,iyp,iz,ib,ik)*tz+avclap(ixm,iyp,izp,ib,ik)*tzp    &
         +avclap(ixm,iyp,izpp,ib,ik)*tzpp+avclap(ixm,iyp,izm,ib,ik)*tzm
    l1(15)=avclap(ixm,iypp,iz,ib,ik)*tz+avclap(ixm,iypp,izp,ib,ik)*tzp  &
         +avclap(ixm,iypp,izpp,ib,ik)*tzpp+avclap(ixm,iypp,izm,ib,ik)*tzm
    l1(16)=avclap(ixm,iym,iz,ib,ik)*tz+avclap(ixm,iym,izp,ib,ik)*tzp    &
         +avclap(ixm,iym,izpp,ib,ik)*tzpp+avclap(ixm,iym,izm,ib,ik)*tzm

    ldum=(l1(1)*ty+l1(2)*typ+l1(3)*typp+l1(4)*tym)*tx+(l1(5)*ty+l1(6)*typ      &
         +l1(7)*typp+l1(8)*tym)*txp+(l1(9)*ty+l1(10)*typ+l1(11)*typp+l1(12)*   &
         tym)*txpp+(l1(13)*ty+l1(14)*typ+l1(15)*typp+l1(16)*tym)*txm

    g1(1)=avcgrad1(ix,iy,iz,ib,ik)*tz+avcgrad1(ix,iy,izp,ib,ik)*tzp+        &
         avcgrad1(ix,iy,izpp,ib,ik)*tzpp+avcgrad1(ix,iy,izm,ib,ik)*tzm
    g1(2)=avcgrad1(ix,iyp,iz,ib,ik)*tz+avcgrad1(ix,iyp,izp,ib,ik)*tzp+      &
         avcgrad1(ix,iyp,izpp,ib,ik)*tzpp+avcgrad1(ix,iyp,izm,ib,ik)*tzm
    g1(3)=avcgrad1(ix,iypp,iz,ib,ik)*tz+avcgrad1(ix,iypp,izp,ib,ik)*tzp+    &
         avcgrad1(ix,iypp,izpp,ib,ik)*tzpp+avcgrad1(ix,iypp,izm,ib,ik)*tzm
    g1(4)=avcgrad1(ix,iym,iz,ib,ik)*tz+avcgrad1(ix,iym,izp,ib,ik)*tzp+      &
         avcgrad1(ix,iym,izpp,ib,ik)*tzpp+avcgrad1(ix,iym,izm,ib,ik)*tzm

    g1(5)=avcgrad1(ixp,iy,iz,ib,ik)*tz+avcgrad1(ixp,iy,izp,ib,ik)*tzp       &
         +avcgrad1(ixp,iy,izpp,ib,ik)*tzpp+avcgrad1(ixp,iy,izm,ib,ik)*tzm
    g1(6)=avcgrad1(ixp,iyp,iz,ib,ik)*tz+avcgrad1(ixp,iyp,izp,ib,ik)*tzp     &
         +avcgrad1(ixp,iyp,izpp,ib,ik)*tzpp+avcgrad1(ixp,iyp,izm,ib,ik)*tzm
    g1(7)=avcgrad1(ixp,iypp,iz,ib,ik)*tz+avcgrad1(ixp,iypp,izp,ib,ik)*tzp   &
         +avcgrad1(ixp,iypp,izpp,ib,ik)*tzpp+avcgrad1(ixp,iypp,izm,ib,ik)*tzm
    g1(8)=avcgrad1(ixp,iym,iz,ib,ik)*tz+avcgrad1(ixp,iym,izp,ib,ik)*tzp     &
         +avcgrad1(ixp,iym,izpp,ib,ik)*tzpp+avcgrad1(ixp,iym,izm,ib,ik)*tzm

    g1(9)=avcgrad1(ixpp,iy,iz,ib,ik)*tz+avcgrad1(ixpp,iy,izp,ib,ik)*tzp     &
         +avcgrad1(ixpp,iy,izpp,ib,ik)*tzpp+avcgrad1(ixpp,iy,izm,ib,ik)*tzm
    g1(10)=avcgrad1(ixpp,iyp,iz,ib,ik)*tz+avcgrad1(ixpp,iyp,izp,ib,ik)*tzp  &
         +avcgrad1(ixpp,iyp,izpp,ib,ik)*tzpp+avcgrad1(ixpp,iyp,izm,ib,ik)*tzm
    g1(11)=avcgrad1(ixpp,iypp,iz,ib,ik)*tz+avcgrad1(ixpp,iypp,izp,ib,ik)*tzp&
         +avcgrad1(ixpp,iypp,izpp,ib,ik)*tzpp+avcgrad1(ixpp,iypp,izm,ib,ik)*tzm
    g1(12)=avcgrad1(ixpp,iym,iz,ib,ik)*tz+avcgrad1(ixpp,iym,izp,ib,ik)*tzp  &
         +avcgrad1(ixpp,iym,izpp,ib,ik)*tzpp+avcgrad1(ixpp,iym,izm,ib,ik)*tzm

    g1(13)=avcgrad1(ixm,iy,iz,ib,ik)*tz+avcgrad1(ixm,iy,izp,ib,ik)*tzp      &
         +avcgrad1(ixm,iy,izpp,ib,ik)*tzpp+avcgrad1(ixm,iy,izm,ib,ik)*tzm
    g1(14)=avcgrad1(ixm,iyp,iz,ib,ik)*tz+avcgrad1(ixm,iyp,izp,ib,ik)*tzp    &
         +avcgrad1(ixm,iyp,izpp,ib,ik)*tzpp+avcgrad1(ixm,iyp,izm,ib,ik)*tzm
    g1(15)=avcgrad1(ixm,iypp,iz,ib,ik)*tz+avcgrad1(ixm,iypp,izp,ib,ik)*tzp  &
         +avcgrad1(ixm,iypp,izpp,ib,ik)*tzpp+avcgrad1(ixm,iypp,izm,ib,ik)*tzm
    g1(16)=avcgrad1(ixm,iym,iz,ib,ik)*tz+avcgrad1(ixm,iym,izp,ib,ik)*tzp    &
         +avcgrad1(ixm,iym,izpp,ib,ik)*tzpp+avcgrad1(ixm,iym,izm,ib,ik)*tzm

    gdum(1)=(g1(1)*ty+g1(2)*typ+g1(3)*typp+g1(4)*tym)*tx+ &
            (g1(5)*ty+g1(6)*typ+g1(7)*typp+g1(8)*tym)*txp+&
            (g1(9)*ty+g1(10)*typ+g1(11)*typp+g1(12)*tym)*txpp+&
            (g1(13)*ty+g1(14)*typ+g1(15)*typp+g1(16)*tym)*txm

    g2(1)=avcgrad2(ix,iy,iz,ib,ik)*tz+avcgrad2(ix,iy,izp,ib,ik)*tzp+        &
         avcgrad2(ix,iy,izpp,ib,ik)*tzpp+avcgrad2(ix,iy,izm,ib,ik)*tzm
    g2(2)=avcgrad2(ix,iyp,iz,ib,ik)*tz+avcgrad2(ix,iyp,izp,ib,ik)*tzp+      &
         avcgrad2(ix,iyp,izpp,ib,ik)*tzpp+avcgrad2(ix,iyp,izm,ib,ik)*tzm
    g2(3)=avcgrad2(ix,iypp,iz,ib,ik)*tz+avcgrad2(ix,iypp,izp,ib,ik)*tzp+    &
         avcgrad2(ix,iypp,izpp,ib,ik)*tzpp+avcgrad2(ix,iypp,izm,ib,ik)*tzm
    g2(4)=avcgrad2(ix,iym,iz,ib,ik)*tz+avcgrad2(ix,iym,izp,ib,ik)*tzp+      &
         avcgrad2(ix,iym,izpp,ib,ik)*tzpp+avcgrad2(ix,iym,izm,ib,ik)*tzm

    g2(5)=avcgrad2(ixp,iy,iz,ib,ik)*tz+avcgrad2(ixp,iy,izp,ib,ik)*tzp       &
         +avcgrad2(ixp,iy,izpp,ib,ik)*tzpp+avcgrad2(ixp,iy,izm,ib,ik)*tzm
    g2(6)=avcgrad2(ixp,iyp,iz,ib,ik)*tz+avcgrad2(ixp,iyp,izp,ib,ik)*tzp     &
         +avcgrad2(ixp,iyp,izpp,ib,ik)*tzpp+avcgrad2(ixp,iyp,izm,ib,ik)*tzm
    g2(7)=avcgrad2(ixp,iypp,iz,ib,ik)*tz+avcgrad2(ixp,iypp,izp,ib,ik)*tzp   &
         +avcgrad2(ixp,iypp,izpp,ib,ik)*tzpp+avcgrad2(ixp,iypp,izm,ib,ik)*tzm
    g2(8)=avcgrad2(ixp,iym,iz,ib,ik)*tz+avcgrad2(ixp,iym,izp,ib,ik)*tzp     &
         +avcgrad2(ixp,iym,izpp,ib,ik)*tzpp+avcgrad2(ixp,iym,izm,ib,ik)*tzm

    g2(9)=avcgrad2(ixpp,iy,iz,ib,ik)*tz+avcgrad2(ixpp,iy,izp,ib,ik)*tzp     &
         +avcgrad2(ixpp,iy,izpp,ib,ik)*tzpp+avcgrad2(ixpp,iy,izm,ib,ik)*tzm
    g2(10)=avcgrad2(ixpp,iyp,iz,ib,ik)*tz+avcgrad2(ixpp,iyp,izp,ib,ik)*tzp  &
         +avcgrad2(ixpp,iyp,izpp,ib,ik)*tzpp+avcgrad2(ixpp,iyp,izm,ib,ik)*tzm
    g2(11)=avcgrad2(ixpp,iypp,iz,ib,ik)*tz+avcgrad2(ixpp,iypp,izp,ib,ik)*tzp&
         +avcgrad2(ixpp,iypp,izpp,ib,ik)*tzpp+avcgrad2(ixpp,iypp,izm,ib,ik)*tzm
    g2(12)=avcgrad2(ixpp,iym,iz,ib,ik)*tz+avcgrad2(ixpp,iym,izp,ib,ik)*tzp  &
         +avcgrad2(ixpp,iym,izpp,ib,ik)*tzpp+avcgrad2(ixpp,iym,izm,ib,ik)*tzm

    g2(13)=avcgrad2(ixm,iy,iz,ib,ik)*tz+avcgrad2(ixm,iy,izp,ib,ik)*tzp      &
         +avcgrad2(ixm,iy,izpp,ib,ik)*tzpp+avcgrad2(ixm,iy,izm,ib,ik)*tzm
    g2(14)=avcgrad2(ixm,iyp,iz,ib,ik)*tz+avcgrad2(ixm,iyp,izp,ib,ik)*tzp    &
         +avcgrad2(ixm,iyp,izpp,ib,ik)*tzpp+avcgrad2(ixm,iyp,izm,ib,ik)*tzm
    g2(15)=avcgrad2(ixm,iypp,iz,ib,ik)*tz+avcgrad2(ixm,iypp,izp,ib,ik)*tzp  &
         +avcgrad2(ixm,iypp,izpp,ib,ik)*tzpp+avcgrad2(ixm,iypp,izm,ib,ik)*tzm
    g2(16)=avcgrad2(ixm,iym,iz,ib,ik)*tz+avcgrad2(ixm,iym,izp,ib,ik)*tzp    &
         +avcgrad2(ixm,iym,izpp,ib,ik)*tzpp+avcgrad2(ixm,iym,izm,ib,ik)*tzm

    gdum(2)=(g2(1)*ty+g2(2)*typ+g2(3)*typp+g2(4)*tym)*tx+ &
            (g2(5)*ty+g2(6)*typ+g2(7)*typp+g2(8)*tym)*txp+&
            (g2(9)*ty+g2(10)*typ+g2(11)*typp+g2(12)*tym)*txpp+&
            (g2(13)*ty+g2(14)*typ+g2(15)*typp+g2(16)*tym)*txm

    g3(1)=avcgrad3(ix,iy,iz,ib,ik)*tz+avcgrad3(ix,iy,izp,ib,ik)*tzp+        &
         avcgrad3(ix,iy,izpp,ib,ik)*tzpp+avcgrad3(ix,iy,izm,ib,ik)*tzm
    g3(2)=avcgrad3(ix,iyp,iz,ib,ik)*tz+avcgrad3(ix,iyp,izp,ib,ik)*tzp+      &
         avcgrad3(ix,iyp,izpp,ib,ik)*tzpp+avcgrad3(ix,iyp,izm,ib,ik)*tzm
    g3(3)=avcgrad3(ix,iypp,iz,ib,ik)*tz+avcgrad3(ix,iypp,izp,ib,ik)*tzp+    &
         avcgrad3(ix,iypp,izpp,ib,ik)*tzpp+avcgrad3(ix,iypp,izm,ib,ik)*tzm
    g3(4)=avcgrad3(ix,iym,iz,ib,ik)*tz+avcgrad3(ix,iym,izp,ib,ik)*tzp+      &
         avcgrad3(ix,iym,izpp,ib,ik)*tzpp+avcgrad3(ix,iym,izm,ib,ik)*tzm

    g3(5)=avcgrad3(ixp,iy,iz,ib,ik)*tz+avcgrad3(ixp,iy,izp,ib,ik)*tzp       &
         +avcgrad3(ixp,iy,izpp,ib,ik)*tzpp+avcgrad3(ixp,iy,izm,ib,ik)*tzm
    g3(6)=avcgrad3(ixp,iyp,iz,ib,ik)*tz+avcgrad3(ixp,iyp,izp,ib,ik)*tzp     &
         +avcgrad3(ixp,iyp,izpp,ib,ik)*tzpp+avcgrad3(ixp,iyp,izm,ib,ik)*tzm
    g3(7)=avcgrad3(ixp,iypp,iz,ib,ik)*tz+avcgrad3(ixp,iypp,izp,ib,ik)*tzp   &
         +avcgrad3(ixp,iypp,izpp,ib,ik)*tzpp+avcgrad3(ixp,iypp,izm,ib,ik)*tzm
    g3(8)=avcgrad3(ixp,iym,iz,ib,ik)*tz+avcgrad3(ixp,iym,izp,ib,ik)*tzp     &
         +avcgrad3(ixp,iym,izpp,ib,ik)*tzpp+avcgrad3(ixp,iym,izm,ib,ik)*tzm

    g3(9)=avcgrad3(ixpp,iy,iz,ib,ik)*tz+avcgrad3(ixpp,iy,izp,ib,ik)*tzp     &
         +avcgrad3(ixpp,iy,izpp,ib,ik)*tzpp+avcgrad3(ixpp,iy,izm,ib,ik)*tzm
    g3(10)=avcgrad3(ixpp,iyp,iz,ib,ik)*tz+avcgrad3(ixpp,iyp,izp,ib,ik)*tzp  &
         +avcgrad3(ixpp,iyp,izpp,ib,ik)*tzpp+avcgrad3(ixpp,iyp,izm,ib,ik)*tzm
    g3(11)=avcgrad3(ixpp,iypp,iz,ib,ik)*tz+avcgrad3(ixpp,iypp,izp,ib,ik)*tzp&
         +avcgrad3(ixpp,iypp,izpp,ib,ik)*tzpp+avcgrad3(ixpp,iypp,izm,ib,ik)*tzm
    g3(12)=avcgrad3(ixpp,iym,iz,ib,ik)*tz+avcgrad3(ixpp,iym,izp,ib,ik)*tzp  &
         +avcgrad3(ixpp,iym,izpp,ib,ik)*tzpp+avcgrad3(ixpp,iym,izm,ib,ik)*tzm

    g3(13)=avcgrad3(ixm,iy,iz,ib,ik)*tz+avcgrad3(ixm,iy,izp,ib,ik)*tzp      &
         +avcgrad3(ixm,iy,izpp,ib,ik)*tzpp+avcgrad3(ixm,iy,izm,ib,ik)*tzm
    g3(14)=avcgrad3(ixm,iyp,iz,ib,ik)*tz+avcgrad3(ixm,iyp,izp,ib,ik)*tzp    &
         +avcgrad3(ixm,iyp,izpp,ib,ik)*tzpp+avcgrad3(ixm,iyp,izm,ib,ik)*tzm
    g3(15)=avcgrad3(ixm,iypp,iz,ib,ik)*tz+avcgrad3(ixm,iypp,izp,ib,ik)*tzp  &
         +avcgrad3(ixm,iypp,izpp,ib,ik)*tzpp+avcgrad3(ixm,iypp,izm,ib,ik)*tzm
    g3(16)=avcgrad3(ixm,iym,iz,ib,ik)*tz+avcgrad3(ixm,iym,izp,ib,ik)*tzp    &
         +avcgrad3(ixm,iym,izpp,ib,ik)*tzpp+avcgrad3(ixm,iym,izm,ib,ik)*tzm

    gdum(3)=(g3(1)*ty+g3(2)*typ+g3(3)*typp+g3(4)*tym)*tx+ &
            (g3(5)*ty+g3(6)*typ+g3(7)*typp+g3(8)*tym)*txp+&
            (g3(9)*ty+g3(10)*typ+g3(11)*typp+g3(12)*tym)*txpp+&
            (g3(13)*ty+g3(14)*typ+g3(15)*typp+g3(16)*tym)*txm

! go to the Cartesian grid
!    gdum=matmul(bg,gdum)

     if(lkedge(ik))then
      if(use_real_part(ib,ik))then
       lap(idum+ib)=dble(lzdum(ik)*rdum+                                       &
            2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+   &
                        zdum(ik)*ldum)
       grad(1,idum+ib)=dble(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
       grad(2,idum+ib)=dble(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
       grad(3,idum+ib)=dble(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
      else
        lap(idum+ib)=aimag(lzdum(ik)*rdum+                                     &
            2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+   &
                        zdum(ik)*ldum)
         grad(1,idum+ib)=aimag(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
         grad(2,idum+ib)=aimag(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
         grad(3,idum+ib)=aimag(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
      endif
     else
       lap(idum+1)=dble(lzdum(ik)*rdum+                                        &
            2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+   &
                         zdum(ik)*ldum)
       lap(idum+2)=aimag(lzdum(ik)*rdum+                                       &
            2*(gdum(1)*gzdum(1,ik)+gdum(2)*gzdum(2,ik)+gdum(3)*gzdum(3,ik))+   &
                          zdum(ik)*ldum)
      grad(1,idum+1)=dble(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
      grad(2,idum+1)=dble(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
      grad(3,idum+1)=dble(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
      grad(1,idum+2)=aimag(gzdum(1,ik)*rdum+zdum(ik)*gdum(1))
      grad(2,idum+2)=aimag(gzdum(2,ik)*rdum+zdum(ik)*gdum(2))
      grad(3,idum+2)=aimag(gzdum(3,ik)*rdum+zdum(ik)*gdum(3))
     endif
    end select
    endif
! Increment idum by 2 if not at BZ edge and by occup of band if at edge set around line 763? Don't follow logic.
    if(.not.lkedge(ik))idum=idum+2
   enddo
   if(lkedge(ik))idum=idum+boccband(ik)
  endif
 enddo

 END SUBROUTINE blip3d


 SUBROUTINE blip3dgamma_w(rpsi,lap,grad,r,avc,nrbwf,bg,igl,boccband,nemax,maxband)

!---------------------------------------------------------------------------!
! This subroutine evaluates the value of a function                         !
! at a vector point r, using the overlapping of blip functions.             !
! The blip grid is defined on a cubic cell, so r should always be given in  !
! units of the crystal lattice vectors.                                     !
!                                                                           !
! Input:                                                                    !
!  r(3)                 position in units of lattice vectors                !
!  avc(0:nrbwf,0:nrbwf,0:nrbwf,maxband)             blip coefficients                !
!  nrbwf(3)                num of divisions for each side of the box           !
!                       (defines the blip grid)                             !
!  bg(3,3)              the reciprocal lattice vectors (in a.u./tpi)        !
!                                                                           !
! GAMMA ONLY                                                                !
! DA 3.2001                                                                 !
! WP 8.2006                                                                 !
!---------------------------------------------------------------------------!

 IMPLICIT NONE
 INTEGER,INTENT(in) :: igl,boccband,maxband,nemax
 INTEGER nrbwf(3),ix,iy,iz,ixp,ixpp,ixm,iyp,iypp,iym,izp,izpp,izm,idum,ib
 REAL(dp) avc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),r(3),t(3),txm,tx,txp,txpp,&
  &tym,ty,typ,typp,tzm,tz,tzp,tzpp,bg(3,3),d1(16),rpsi(nemax),lap(nemax),      &
  &grad(3,nemax),x,y,z,b(6)

 ix=int(mod(r(1)+abs(int(r(1)))+1,1.d0)*nrbwf(1))
 iy=int(mod(r(2)+abs(int(r(2)))+1,1.d0)*nrbwf(2))
 iz=int(mod(r(3)+abs(int(r(3)))+1,1.d0)*nrbwf(3))

 if(ix<0.or.iy<0.or.iz<0)call errstop('BLIP3D','Negative index found.')

! The blips are defined as the product of one-dimensional cubic splines, These
! are different from zero only on four adjacent grid points.  It follows that
! for any vector r there are only 64 overlapping three-dimensional blips, which
! are the product of all possible combinations of the three one-dimensional
! splines.

! These are the extra 3 coefficients for each dimension needed
 ix=mod(ix,nrbwf(1))
 iy=mod(iy,nrbwf(2))
 iz=mod(iz,nrbwf(3))
 ixp=mod(ix+1,nrbwf(1))
 ixpp=mod(ix+2,nrbwf(1))
 ixm=mod(ix-1+nrbwf(1),nrbwf(1))
 iyp=mod(iy+1,nrbwf(2))
 iypp=mod(iy+2,nrbwf(2))
 iym=mod(iy-1+nrbwf(2),nrbwf(2))
 izp=mod(iz+1,nrbwf(3))
 izpp=mod(iz+2,nrbwf(3))
 izm=mod(iz-1+nrbwf(3),nrbwf(3))

! Now calculate the 12 monodimensional blip functions
 t=mod(r+abs(int(r))+1,1.d0)*nrbwf
 t=mod(t,dble(nrbwf))

 x=t(1)-ix+1.d0
 txm=2.d0-3*x+1.5d0*x*x-0.25d0*x*x*x
 x=t(1)-ix
 tx=1.d0-1.5d0*x*x+0.75d0*x*x*x
 x=t(1)-ix-1.d0
 txp=1.d0-1.5d0*x*x-0.75d0*x*x*x
 x=t(1)-ix-2.d0
 txpp=2.d0+3*x+1.5d0*x*x+0.25d0*x*x*x

 y=t(2)-iy+1.d0
 tym=2.d0-3*y+1.5d0*y*y-0.25d0*y*y*y
 y=t(2)-iy
 ty=1.d0-1.5d0*y*y+0.75d0*y*y*y
 y=t(2)-iy-1.d0
 typ=1.d0-1.5d0*y*y-0.75d0*y*y*y
 y=t(2)-iy-2.d0
 typp=2.d0+3*y+1.5d0*y*y+0.25d0*y*y*y

 z=t(3)-iz+1.d0
 tzm=2.d0-3*z+1.5d0*z*z-0.25d0*z*z*z
 z=t(3)-iz
 tz=1.d0-1.5d0*z*z+0.75d0*z*z*z
 z=t(3)-iz-1.d0
 tzp=1.d0-1.5d0*z*z-0.75d0*z*z*z
 z=t(3)-iz-2.d0
 tzpp=2.d0+3*z+1.5d0*z*z+0.25d0*z*z*z

 b(1)=(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2)
 b(2)=(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2)
 b(3)=(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2)
 b(4)=2*(bg(1,1)*bg(1,2)+bg(2,1)*bg(2,2)+bg(3,1)*bg(3,2))
 b(5)=2*(bg(1,2)*bg(1,3)+bg(2,2)*bg(2,3)+bg(3,2)*bg(3,3))
 b(6)=2*(bg(1,3)*bg(1,1)+bg(2,3)*bg(2,1)+bg(3,3)*bg(3,1))

 idum=0
 do ib=1,boccband
  idum=idum+1
  d1(1)=avc(ix,iy,iz,ib)*tz+avc(ix,iy,izp,ib)*tzp+                       &
   &avc(ix,iy,izpp,ib)*tzpp+avc(ix,iy,izm,ib)*tzm
  d1(2)=avc(ix,iyp,iz,ib)*tz+avc(ix,iyp,izp,ib)*tzp+                     &
   &avc(ix,iyp,izpp,ib)*tzpp+avc(ix,iyp,izm,ib)*tzm
  d1(3)=avc(ix,iypp,iz,ib)*tz+avc(ix,iypp,izp,ib)*tzp+                   &
   &avc(ix,iypp,izpp,ib)*tzpp+avc(ix,iypp,izm,ib)*tzm
  d1(4)=avc(ix,iym,iz,ib)*tz+avc(ix,iym,izp,ib)*tzp+                     &
   &avc(ix,iym,izpp,ib)*tzpp+avc(ix,iym,izm,ib)*tzm

  d1(5)=avc(ixp,iy,iz,ib)*tz+avc(ixp,iy,izp,ib)*tzp                      &
   &+avc(ixp,iy,izpp,ib)*tzpp+avc(ixp,iy,izm,ib)*tzm
  d1(6)=avc(ixp,iyp,iz,ib)*tz+avc(ixp,iyp,izp,ib)*tzp                    &
   &+avc(ixp,iyp,izpp,ib)*tzpp+avc(ixp,iyp,izm,ib)*tzm
  d1(7)=avc(ixp,iypp,iz,ib)*tz+avc(ixp,iypp,izp,ib)*tzp                  &
   &+avc(ixp,iypp,izpp,ib)*tzpp+avc(ixp,iypp,izm,ib)*tzm
  d1(8)=avc(ixp,iym,iz,ib)*tz+avc(ixp,iym,izp,ib)*tzp                    &
   &+avc(ixp,iym,izpp,ib)*tzpp+avc(ixp,iym,izm,ib)*tzm

  d1(9)=avc(ixpp,iy,iz,ib)*tz+avc(ixpp,iy,izp,ib)*tzp                    &
   &+avc(ixpp,iy,izpp,ib)*tzpp+avc(ixpp,iy,izm,ib)*tzm
  d1(10)=avc(ixpp,iyp,iz,ib)*tz+avc(ixpp,iyp,izp,ib)*tzp                 &
   &+avc(ixpp,iyp,izpp,ib)*tzpp+avc(ixpp,iyp,izm,ib)*tzm
  d1(11)=avc(ixpp,iypp,iz,ib)*tz+avc(ixpp,iypp,izp,ib)*tzp               &
   &+avc(ixpp,iypp,izpp,ib)*tzpp+avc(ixpp,iypp,izm,ib)*tzm
  d1(12)=avc(ixpp,iym,iz,ib)*tz+avc(ixpp,iym,izp,ib)*tzp                 &
   &+avc(ixpp,iym,izpp,ib)*tzpp+avc(ixpp,iym,izm,ib)*tzm

  d1(13)=avc(ixm,iy,iz,ib)*tz+avc(ixm,iy,izp,ib)*tzp                     &
   &+avc(ixm,iy,izpp,ib)*tzpp+avc(ixm,iy,izm,ib)*tzm
  d1(14)=avc(ixm,iyp,iz,ib)*tz+avc(ixm,iyp,izp,ib)*tzp                   &
   &+avc(ixm,iyp,izpp,ib)*tzpp+avc(ixm,iyp,izm,ib)*tzm
  d1(15)=avc(ixm,iypp,iz,ib)*tz+avc(ixm,iypp,izp,ib)*tzp                 &
   &+avc(ixm,iypp,izpp,ib)*tzpp+avc(ixm,iypp,izm,ib)*tzm
  d1(16)=avc(ixm,iym,iz,ib)*tz+avc(ixm,iym,izp,ib)*tzp                   &
   &+avc(ixm,iym,izpp,ib)*tzpp+avc(ixm,iym,izm,ib)*tzm

! The function
  rpsi(idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*typ &
   &+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*  &
   &txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

 enddo

 END SUBROUTINE blip3dgamma_w

 SUBROUTINE blip3dgamma_gl(rpsi,lap,grad,r,avc,avclap,avcgrad1,avcgrad2,avcgrad3,&
 &nrbwf,bg,igl,boccband,nemax,maxband)

!---------------------------------------------------------------------------!
! This subroutine evaluates the value of a function,its gradient and        !
! Laplacian at a vector point r, using the overlapping of blip functions.   !
! The blip grid is defined on a cubic cell, so r should always be given in  !
! units of the crystal lattice vectors.                                     !
!                                                                           !
! Input:                                                                    !
!  r(3)                 position in units of lattice vectors                !
!  avc(0:nrbwf,0:nrbwf,0:nrbwf,maxband)             blip coefficients                !
!  nrbwf(3)                num of divisions for each side of the box           !
!                       (defines the blip grid)                             !
!  bg(3,3)              the reciprocal lattice vectors (in a.u./pi)        !
!                                                                           !
! GAMMA ONLY                                                                !
! DA 3.2001                                                                 !
! WP 8.2006                                                                 !
!---------------------------------------------------------------------------!

 IMPLICIT NONE
 INTEGER,INTENT(in) :: igl,boccband,maxband,nemax
 INTEGER nrbwf(3),ix,iy,iz,ixp,ixpp,ixm,iyp,iypp,iym,izp,izpp,izm,idum,ib
 REAL(dp) avc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),r(3),t(3),txm,tx,txp,txpp,&
  &tym,ty,typ,typp,tzm,tz,tzp,tzpp,bg(3,3),d1(16),rpsi(nemax),lap(nemax),      &
  &grad(3,nemax),x,y,z,b(6)
 REAL(dp) avclap(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband)
 REAL(dp) avcgrad1(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband), &
  &avcgrad2(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband),        &
  &avcgrad3(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1,maxband)
 REAL(dp) :: dtxm=0.d0,d2txm=0.d0,dtx=0.d0,d2tx=0.d0,dtxp=0.d0,d2txp=0.d0,     &
  &dtxpp=0.d0,d2txpp=0.d0,dtym=0.d0,d2tym=0.d0,dty=0.d0,d2ty=0.d0,dtyp=0.d0,   &
  &d2typ=0.d0,dtypp=0.d0,d2typp=0.d0,dtzm=0.d0,d2tzm=0.d0,dtz=0.d0,d2tz=0.d0,  &
  &dtzp=0.d0,d2tzp=0.d0,dtzpp=0.d0,d2tzpp=0.d0

 ix=int(mod(r(1)+abs(int(r(1)))+1,1.d0)*nrbwf(1))
 iy=int(mod(r(2)+abs(int(r(2)))+1,1.d0)*nrbwf(2))
 iz=int(mod(r(3)+abs(int(r(3)))+1,1.d0)*nrbwf(3))

 if(ix<0.or.iy<0.or.iz<0)call errstop('BLIP3D','Negative index found.')

! The blips are defined as the product of one-dimensional cubic splines, These
! are different from zero only on four adjacent grid points.  It follows that
! for any vector r there are only 64 overlapping three-dimensional blips, which
! are the product of all possible combinations of the three one-dimensional
! splines.

! These are the extra 3 coefficients for each dimension needed
 ix=mod(ix,nrbwf(1))
 iy=mod(iy,nrbwf(2))
 iz=mod(iz,nrbwf(3))
 ixp=mod(ix+1,nrbwf(1))
 ixpp=mod(ix+2,nrbwf(1))
 ixm=mod(ix-1+nrbwf(1),nrbwf(1))
 iyp=mod(iy+1,nrbwf(2))
 iypp=mod(iy+2,nrbwf(2))
 iym=mod(iy-1+nrbwf(2),nrbwf(2))
 izp=mod(iz+1,nrbwf(3))
 izpp=mod(iz+2,nrbwf(3))
 izm=mod(iz-1+nrbwf(3),nrbwf(3))

! Now calculate the 12 monodimensional blip functions
 t=mod(r+abs(int(r))+1,1.d0)*nrbwf
 t=mod(t,dble(nrbwf))

! write(6,*)"In bwfdet_mod, igrad_lap=",igrad_lap

 select case(igrad_lap)
 case (0)
   x=t(1)-ix+1.d0
   txm=2.d0-3*x+1.5d0*x*x-0.25d0*x*x*x
    dtxm=(-3.d0+3*x-0.75d0*x*x)*nrbwf(1)
    d2txm=(3.d0-1.5d0*x)*nrbwf(1)*nrbwf(1)
   x=t(1)-ix
   tx=1.d0-1.5d0*x*x+0.75d0*x*x*x
    dtx=(-3.d0*x+2.25d0*x*x)*nrbwf(1)
    d2tx=(-3.d0+4.5d0*x)*nrbwf(1)*nrbwf(1)
   x=t(1)-ix-1.d0
   txp=1.d0-1.5d0*x*x-0.75d0*x*x*x
    dtxp=(-3.d0*x-2.25d0*x*x)*nrbwf(1)
    d2txp=(-3.d0-4.5d0*x)*nrbwf(1)*nrbwf(1)
   x=t(1)-ix-2.d0
   txpp=2.d0+3*x+1.5d0*x*x+0.25d0*x*x*x
    dtxpp=(3.d0+3*x+0.75d0*x*x)*nrbwf(1)
    d2txpp=(3.d0+1.5d0*x)*nrbwf(1)*nrbwf(1)

   y=t(2)-iy+1.d0
   tym=2.d0-3*y+1.5d0*y*y-0.25d0*y*y*y
    dtym=(-3.d0+3*y-0.75d0*y*y)*nrbwf(2)
    d2tym=(3.d0-1.5d0*y)*nrbwf(2)*nrbwf(2)
   y=t(2)-iy
   ty=1.d0-1.5d0*y*y+0.75d0*y*y*y
    dty=(-3.d0*y+2.25d0*y*y)*nrbwf(2)
    d2ty=(-3.d0+4.5d0*y)*nrbwf(2)*nrbwf(2)
   y=t(2)-iy-1.d0
   typ=1.d0-1.5d0*y*y-0.75d0*y*y*y
    dtyp=(-3.d0*y-2.25d0*y*y)*nrbwf(2)
    d2typ=(-3.d0-4.5d0*y)*nrbwf(2)*nrbwf(2)
   y=t(2)-iy-2.d0
   typp=2.d0+3*y+1.5d0*y*y+0.25d0*y*y*y
    dtypp=(3.d0+3*y+0.75d0*y*y)*nrbwf(2)
    d2typp=(3.d0+1.5d0*y)*nrbwf(2)*nrbwf(2)

   z=t(3)-iz+1.d0
   tzm=2.d0-3*z+1.5d0*z*z-0.25d0*z*z*z
    dtzm=(-3.d0+3*z-0.75d0*z*z)*nrbwf(3)
    d2tzm=(3.d0-1.5d0*z)*nrbwf(3)*nrbwf(3)
   z=t(3)-iz
   tz=1.d0-1.5d0*z*z+0.75d0*z*z*z
    dtz=(-3.d0*z+2.25d0*z*z)*nrbwf(3)
    d2tz=(-3.d0+4.5d0*z)*nrbwf(3)*nrbwf(3)
   z=t(3)-iz-1.d0
   tzp=1.d0-1.5d0*z*z-0.75d0*z*z*z
    dtzp=(-3.d0*z-2.25d0*z*z)*nrbwf(3)
    d2tzp=(-3.d0-4.5d0*z)*nrbwf(3)*nrbwf(3)
   z=t(3)-iz-2.d0
   tzpp=2.d0+3*z+1.5d0*z*z+0.25d0*z*z*z
    dtzpp=(3.d0+3*z+0.75d0*z*z)*nrbwf(3)
    d2tzpp=(3.d0+1.5d0*z)*nrbwf(3)*nrbwf(3)
 case (1)
   x=t(1)-ix+1.d0
   txm=2.d0-3*x+1.5d0*x*x-0.25d0*x*x*x
    dtxm=(-3.d0+3*x-0.75d0*x*x)*nrbwf(1)
   x=t(1)-ix
   tx=1.d0-1.5d0*x*x+0.75d0*x*x*x
    dtx=(-3.d0*x+2.25d0*x*x)*nrbwf(1)
   x=t(1)-ix-1.d0
   txp=1.d0-1.5d0*x*x-0.75d0*x*x*x
    dtxp=(-3.d0*x-2.25d0*x*x)*nrbwf(1)
   x=t(1)-ix-2.d0
   txpp=2.d0+3*x+1.5d0*x*x+0.25d0*x*x*x
    dtxpp=(3.d0+3*x+0.75d0*x*x)*nrbwf(1)

   y=t(2)-iy+1.d0
   tym=2.d0-3*y+1.5d0*y*y-0.25d0*y*y*y
    dtym=(-3.d0+3*y-0.75d0*y*y)*nrbwf(2)
   y=t(2)-iy
   ty=1.d0-1.5d0*y*y+0.75d0*y*y*y
    dty=(-3.d0*y+2.25d0*y*y)*nrbwf(2)
   y=t(2)-iy-1.d0
   typ=1.d0-1.5d0*y*y-0.75d0*y*y*y
    dtyp=(-3.d0*y-2.25d0*y*y)*nrbwf(2)
   y=t(2)-iy-2.d0
   typp=2.d0+3*y+1.5d0*y*y+0.25d0*y*y*y
    dtypp=(3.d0+3*y+0.75d0*y*y)*nrbwf(2)

   z=t(3)-iz+1.d0
   tzm=2.d0-3*z+1.5d0*z*z-0.25d0*z*z*z
    dtzm=(-3.d0+3*z-0.75d0*z*z)*nrbwf(3)
   z=t(3)-iz
   tz=1.d0-1.5d0*z*z+0.75d0*z*z*z
    dtz=(-3.d0*z+2.25d0*z*z)*nrbwf(3)
   z=t(3)-iz-1.d0
   tzp=1.d0-1.5d0*z*z-0.75d0*z*z*z
    dtzp=(-3.d0*z-2.25d0*z*z)*nrbwf(3)
   z=t(3)-iz-2.d0
   tzpp=2.d0+3*z+1.5d0*z*z+0.25d0*z*z*z
    dtzpp=(3.d0+3*z+0.75d0*z*z)*nrbwf(3)
 case (2)
   x=t(1)-ix+1.d0
   txm=2.d0-3*x+1.5d0*x*x-0.25d0*x*x*x
   x=t(1)-ix
   tx=1.d0-1.5d0*x*x+0.75d0*x*x*x
   x=t(1)-ix-1.d0
   txp=1.d0-1.5d0*x*x-0.75d0*x*x*x
   x=t(1)-ix-2.d0
   txpp=2.d0+3*x+1.5d0*x*x+0.25d0*x*x*x

   y=t(2)-iy+1.d0
   tym=2.d0-3*y+1.5d0*y*y-0.25d0*y*y*y
   y=t(2)-iy
   ty=1.d0-1.5d0*y*y+0.75d0*y*y*y
   y=t(2)-iy-1.d0
   typ=1.d0-1.5d0*y*y-0.75d0*y*y*y
   y=t(2)-iy-2.d0
   typp=2.d0+3*y+1.5d0*y*y+0.25d0*y*y*y

   z=t(3)-iz+1.d0
   tzm=2.d0-3*z+1.5d0*z*z-0.25d0*z*z*z
   z=t(3)-iz
   tz=1.d0-1.5d0*z*z+0.75d0*z*z*z
   z=t(3)-iz-1.d0
   tzp=1.d0-1.5d0*z*z-0.75d0*z*z*z
   z=t(3)-iz-2.d0
   tzpp=2.d0+3*z+1.5d0*z*z+0.25d0*z*z*z
 end select

 b(1)=(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2)
 b(2)=(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2)
 b(3)=(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2)
 b(4)=2*(bg(1,1)*bg(1,2)+bg(2,1)*bg(2,2)+bg(3,1)*bg(3,2))
 b(5)=2*(bg(1,2)*bg(1,3)+bg(2,2)*bg(2,3)+bg(3,2)*bg(3,3))
 b(6)=2*(bg(1,3)*bg(1,1)+bg(2,3)*bg(2,1)+bg(3,3)*bg(3,1))

 idum=0
 do ib=1,boccband
  idum=idum+1
  d1(1)=avc(ix,iy,iz,ib)*tz+avc(ix,iy,izp,ib)*tzp+                       &
   &avc(ix,iy,izpp,ib)*tzpp+avc(ix,iy,izm,ib)*tzm
  d1(2)=avc(ix,iyp,iz,ib)*tz+avc(ix,iyp,izp,ib)*tzp+                     &
   &avc(ix,iyp,izpp,ib)*tzpp+avc(ix,iyp,izm,ib)*tzm
  d1(3)=avc(ix,iypp,iz,ib)*tz+avc(ix,iypp,izp,ib)*tzp+                   &
   &avc(ix,iypp,izpp,ib)*tzpp+avc(ix,iypp,izm,ib)*tzm
  d1(4)=avc(ix,iym,iz,ib)*tz+avc(ix,iym,izp,ib)*tzp+                     &
   &avc(ix,iym,izpp,ib)*tzpp+avc(ix,iym,izm,ib)*tzm

  d1(5)=avc(ixp,iy,iz,ib)*tz+avc(ixp,iy,izp,ib)*tzp                      &
   &+avc(ixp,iy,izpp,ib)*tzpp+avc(ixp,iy,izm,ib)*tzm
  d1(6)=avc(ixp,iyp,iz,ib)*tz+avc(ixp,iyp,izp,ib)*tzp                    &
   &+avc(ixp,iyp,izpp,ib)*tzpp+avc(ixp,iyp,izm,ib)*tzm
  d1(7)=avc(ixp,iypp,iz,ib)*tz+avc(ixp,iypp,izp,ib)*tzp                  &
   &+avc(ixp,iypp,izpp,ib)*tzpp+avc(ixp,iypp,izm,ib)*tzm
  d1(8)=avc(ixp,iym,iz,ib)*tz+avc(ixp,iym,izp,ib)*tzp                    &
   &+avc(ixp,iym,izpp,ib)*tzpp+avc(ixp,iym,izm,ib)*tzm

  d1(9)=avc(ixpp,iy,iz,ib)*tz+avc(ixpp,iy,izp,ib)*tzp                    &
   &+avc(ixpp,iy,izpp,ib)*tzpp+avc(ixpp,iy,izm,ib)*tzm
  d1(10)=avc(ixpp,iyp,iz,ib)*tz+avc(ixpp,iyp,izp,ib)*tzp                 &
   &+avc(ixpp,iyp,izpp,ib)*tzpp+avc(ixpp,iyp,izm,ib)*tzm
  d1(11)=avc(ixpp,iypp,iz,ib)*tz+avc(ixpp,iypp,izp,ib)*tzp               &
   &+avc(ixpp,iypp,izpp,ib)*tzpp+avc(ixpp,iypp,izm,ib)*tzm
  d1(12)=avc(ixpp,iym,iz,ib)*tz+avc(ixpp,iym,izp,ib)*tzp                 &
   &+avc(ixpp,iym,izpp,ib)*tzpp+avc(ixpp,iym,izm,ib)*tzm

  d1(13)=avc(ixm,iy,iz,ib)*tz+avc(ixm,iy,izp,ib)*tzp                     &
   &+avc(ixm,iy,izpp,ib)*tzpp+avc(ixm,iy,izm,ib)*tzm
  d1(14)=avc(ixm,iyp,iz,ib)*tz+avc(ixm,iyp,izp,ib)*tzp                   &
   &+avc(ixm,iyp,izpp,ib)*tzpp+avc(ixm,iyp,izm,ib)*tzm
  d1(15)=avc(ixm,iypp,iz,ib)*tz+avc(ixm,iypp,izp,ib)*tzp                 &
   &+avc(ixm,iypp,izpp,ib)*tzpp+avc(ixm,iypp,izm,ib)*tzm
  d1(16)=avc(ixm,iym,iz,ib)*tz+avc(ixm,iym,izp,ib)*tzp                   &
   &+avc(ixm,iym,izpp,ib)*tzpp+avc(ixm,iym,izm,ib)*tzm

! The function
  rpsi(idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*typ &
   &+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*  &
   &txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

  if(igl==1)then
   select case (igrad_lap)
   case(0)
! The Laplacian: first term involving Theta''(x)
   lap(idum)=((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*d2tx+(d1(5)*ty+d1(6)* &
    &typ+d1(7)*typp+d1(8)*tym)*d2txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)  &
    &*tym)*d2txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*d2txm)*b(1)

! The Laplacian: term involving Theta'(x)Theta'(y)
   lap(idum)=lap(idum)+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*dtx+    &
    &(d1(5)*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*dtxp+(d1(9)*dty+d1(10)*    &
    &dtyp+d1(11)*dtypp+d1(12)*dtym)*dtxpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*     &
    &dtypp+d1(16)*dtym)*dtxm)*b(4)

! The gradient, first term involving Theta'(x)
   grad(1,idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+(d1(5)*ty+d1(6)*&
    &typ+d1(7)*typp+d1(8)*tym)*dtxp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*  &
    &tym)*dtxpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm

! Second term of the laplacian involving Theta''(y)
   lap(idum)=lap(idum)+((d1(1)*d2ty+d1(2)*d2typ+d1(3)*d2typp+d1(4)*d2tym)*tx+ &
    &(d1(5)*d2ty+d1(6)*d2typ+d1(7)*d2typp+d1(8)*d2tym)*txp+(d1(9)*d2ty+       &
    &d1(10)*d2typ+d1(11)*d2typp+d1(12)*d2tym)*txpp+(d1(13)*d2ty+d1(14)*d2typ  &
    &+d1(15)*d2typp+d1(16)*d2tym)*txm)*b(2)

! Second term of the gradient involving Theta'(y)
   grad(2,idum)=(d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+(d1(5)*dty+  &
    &d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+(d1(9)*dty+d1(10)*dtyp+d1(11)*    &
    &dtypp+d1(12)*dtym)*txpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*     &
    &dtym)*txm

! And now the third term of the laplacian involving Theta''(z)
   d1(1)=avc(ix,iy,iz,ib)*d2tz+avc(ix,iy,izp,ib)*d2tzp+                  &
    &avc(ix,iy,izpp,ib)*d2tzpp+avc(ix,iy,izm,ib)*d2tzm
   d1(2)=avc(ix,iyp,iz,ib)*d2tz+avc(ix,iyp,izp,ib)*d2tzp+                &
    &avc(ix,iyp,izpp,ib)*d2tzpp+avc(ix,iyp,izm,ib)*d2tzm
   d1(3)=avc(ix,iypp,iz,ib)*d2tz+avc(ix,iypp,izp,ib)*d2tzp+              &
    &avc(ix,iypp,izpp,ib)*d2tzpp+avc(ix,iypp,izm,ib)*d2tzm
   d1(4)=avc(ix,iym,iz,ib)*d2tz+avc(ix,iym,izp,ib)*d2tzp+                &
    &avc(ix,iym,izpp,ib)*d2tzpp+avc(ix,iym,izm,ib)*d2tzm

   d1(5)=avc(ixp,iy,iz,ib)*d2tz+avc(ixp,iy,izp,ib)*d2tzp+                &
    &avc(ixp,iy,izpp,ib)*d2tzpp+avc(ixp,iy,izm,ib)*d2tzm
   d1(6)=avc(ixp,iyp,iz,ib)*d2tz+avc(ixp,iyp,izp,ib)*d2tzp+              &
    &avc(ixp,iyp,izpp,ib)*d2tzpp+avc(ixp,iyp,izm,ib)*d2tzm
   d1(7)=avc(ixp,iypp,iz,ib)*d2tz+avc(ixp,iypp,izp,ib)*d2tzp+            &
    &avc(ixp,iypp,izpp,ib)*d2tzpp+avc(ixp, iypp, izm,ib)*d2tzm
   d1(8)=avc(ixp,iym,iz,ib)*d2tz+avc(ixp,iym,izp,ib)*d2tzp+              &
    &avc(ixp,iym,izpp,ib)*d2tzpp+avc(ixp,iym,izm,ib)*d2tzm

   d1(9)=avc(ixpp,iy,iz,ib)*d2tz+avc(ixpp,iy,izp,ib)*d2tzp+              &
    &avc(ixpp,iy,izpp,ib)*d2tzpp+avc(ixpp,iy,izm,ib)*d2tzm
   d1(10)=avc(ixpp,iyp,iz,ib)*d2tz+avc(ixpp,iyp,izp,ib)*d2tzp+           &
    &avc(ixpp,iyp,izpp,ib)*d2tzpp+avc(ixpp,iyp,izm,ib)*d2tzm
   d1(11)=avc(ixpp,iypp,iz,ib)*d2tz+avc(ixpp,iypp,izp,ib)*d2tzp+         &
    &avc(ixpp,iypp,izpp,ib)*d2tzpp+avc(ixpp,iypp,izm,ib)*d2tzm
   d1(12)=avc(ixpp,iym,iz,ib)*d2tz+avc(ixpp,iym,izp,ib)*d2tzp+           &
    &avc(ixpp,iym,izpp,ib)*d2tzpp+avc(ixpp,iym,izm,ib)*d2tzm

   d1(13)=avc(ixm,iy,iz,ib)*d2tz+avc(ixm,iy,izp,ib)*d2tzp+               &
    &avc(ixm,iy,izpp,ib)*d2tzpp+avc(ixm,iy,izm,ib)*d2tzm
   d1(14)=avc(ixm,iyp,iz,ib)*d2tz+avc(ixm,iyp,izp,ib)*d2tzp+             &
    &avc(ixm,iyp,izpp,ib)*d2tzpp+avc(ixm,iyp,izm,ib)*d2tzm
   d1(15)=avc(ixm,iypp,iz,ib)*d2tz+avc(ixm,iypp,izp,ib)*d2tzp+           &
    &avc(ixm,iypp,izpp,ib)*d2tzpp+avc(ixm,iypp,izm,ib)*d2tzm
   d1(16)=avc(ixm,iym,iz,ib)*d2tz+avc(ixm,iym,izp,ib)*d2tzp+             &
    &avc(ixm,iym,izpp,ib)*d2tzpp+avc(ixm,iym,izm,ib)*d2tzm

! Theta''(z)
   lap(idum)=lap(idum)+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+&
    &d1(6)*typ+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+     &
    &d1(12)*tym)*txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm)*b(3)

! And the third term of the gradient involving Theta'(z)
   d1(1)=avc(ix,iy,iz,ib)*dtz+avc(ix,iy,izp,ib)*dtzp+                    &
    &avc(ix,iy,izpp,ib)*dtzpp+avc(ix,iy,izm,ib)*dtzm
   d1(2)=avc(ix,iyp,iz,ib)*dtz+avc(ix,iyp,izp,ib)*dtzp+                  &
    &avc(ix,iyp,izpp,ib)*dtzpp+avc(ix,iyp,izm,ib)*dtzm
   d1(3)=avc(ix,iypp,iz,ib)*dtz+avc(ix,iypp,izp,ib)*dtzp+                &
    &avc(ix,iypp,izpp,ib)*dtzpp+avc(ix,iypp,izm,ib)*dtzm
   d1(4)=avc(ix,iym,iz,ib)*dtz+avc(ix,iym,izp,ib)*dtzp+                  &
    &avc(ix,iym,izpp,ib)*dtzpp+avc(ix,iym,izm,ib)*dtzm

   d1(5)=avc(ixp,iy,iz,ib)*dtz+avc(ixp,iy,izp,ib)*dtzp+                  &
    &avc(ixp,iy,izpp,ib)*dtzpp+avc(ixp,iy,izm,ib)*dtzm
   d1(6)=avc(ixp,iyp,iz,ib)*dtz+avc(ixp,iyp,izp,ib)*dtzp+                &
    &avc(ixp,iyp,izpp,ib)*dtzpp+avc(ixp,iyp,izm,ib)*dtzm
   d1(7)=avc(ixp,iypp,iz,ib)*dtz+avc(ixp,iypp,izp,ib)*dtzp+              &
    &avc(ixp,iypp,izpp,ib)*dtzpp+avc(ixp,iypp,izm,ib)*dtzm
   d1(8)=avc(ixp,iym,iz,ib)*dtz+avc(ixp,iym,izp,ib)*dtzp+                &
    &avc(ixp,iym,izpp,ib)*dtzpp+avc(ixp,iym,izm,ib)*dtzm

   d1(9)=avc(ixpp,iy,iz,ib)*dtz+avc(ixpp,iy,izp,ib)*dtzp+                &
    &avc(ixpp,iy,izpp,ib)*dtzpp+avc(ixpp,iy,izm,ib)*dtzm
   d1(10)=avc(ixpp,iyp,iz,ib)*dtz+avc(ixpp,iyp,izp,ib)*dtzp+             &
    &avc(ixpp,iyp,izpp,ib)*dtzpp+avc(ixpp,iyp,izm,ib)*dtzm
   d1(11)=avc(ixpp,iypp,iz,ib)*dtz+avc(ixpp,iypp,izp,ib)*dtzp+           &
    &avc(ixpp,iypp,izpp,ib)*dtzpp+avc(ixpp,iypp,izm,ib)*dtzm
   d1(12)=avc(ixpp,iym,iz,ib)*dtz+avc(ixpp,iym,izp,ib)*dtzp+             &
    &avc(ixpp,iym,izpp,ib)*dtzpp+avc(ixpp,iym,izm,ib)*dtzm

   d1(13)=avc(ixm,iy,iz,ib)*dtz+avc(ixm,iy,izp,ib)*dtzp+                 &
    &avc(ixm,iy,izpp,ib)*dtzpp+avc(ixm,iy,izm,ib)*dtzm
   d1(14)=avc(ixm,iyp,iz,ib)*dtz+avc(ixm,iyp,izp,ib)*dtzp+               &
    &avc(ixm,iyp,izpp,ib)*dtzpp+avc(ixm,iyp,izm,ib)*dtzm
   d1(15)=avc(ixm,iypp,iz,ib)*dtz+avc(ixm,iypp,izp,ib)*dtzp+             &
    &avc(ixm,iypp,izpp,ib)*dtzpp+avc(ixm,iypp,izm,ib)*dtzm
   d1(16)=avc(ixm,iym,iz,ib)*dtz+avc(ixm,iym,izp,ib)*dtzp+               &
    &avc(ixm,iym,izpp,ib)*dtzpp+avc(ixm,iym,izm,ib)*dtzm

! The Laplacian: term involving Theta'(x)Theta'(z)
   lap(idum)=lap(idum)+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+(d1(5)*ty&
    &+d1(6)*typ+d1(7)*typp+d1(8)*tym)*dtxp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+   &
    &d1(12)*tym)*dtxpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm)*b(6)

! the Laplacian: term involving Theta'(y)Theta'(z)
   lap(idum)=lap(idum)+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+(d1(5)&
    &*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+(d1(9)*dty+d1(10)*dtyp+    &
    &d1(11)*dtypp+d1(12)*dtym)*txpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+   &
    &d1(16)*dtym)*txm)*b(5)

   grad(3,idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*  &
    &typ+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)* &
    &tym)*txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

! Go to the Cartesian grid
   grad(:,idum)=matmul(bg,grad(:,idum))

   case(1)

! The gradient, first term involving Theta'(x)
   grad(1,idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+(d1(5)*ty+d1(6)*&
    &typ+d1(7)*typp+d1(8)*tym)*dtxp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*  &
    &tym)*dtxpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm

! Second term of the gradient involving Theta'(y)
   grad(2,idum)=(d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+(d1(5)*dty+  &
    &d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+(d1(9)*dty+d1(10)*dtyp+d1(11)*    &
    &dtypp+d1(12)*dtym)*txpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*     &
    &dtym)*txm

! And the third term of the gradient involving Theta'(z)
   d1(1)=avc(ix,iy,iz,ib)*dtz+avc(ix,iy,izp,ib)*dtzp+                    &
    &avc(ix,iy,izpp,ib)*dtzpp+avc(ix,iy,izm,ib)*dtzm
   d1(2)=avc(ix,iyp,iz,ib)*dtz+avc(ix,iyp,izp,ib)*dtzp+                  &
    &avc(ix,iyp,izpp,ib)*dtzpp+avc(ix,iyp,izm,ib)*dtzm
   d1(3)=avc(ix,iypp,iz,ib)*dtz+avc(ix,iypp,izp,ib)*dtzp+                &
    &avc(ix,iypp,izpp,ib)*dtzpp+avc(ix,iypp,izm,ib)*dtzm
   d1(4)=avc(ix,iym,iz,ib)*dtz+avc(ix,iym,izp,ib)*dtzp+                  &
    &avc(ix,iym,izpp,ib)*dtzpp+avc(ix,iym,izm,ib)*dtzm

   d1(5)=avc(ixp,iy,iz,ib)*dtz+avc(ixp,iy,izp,ib)*dtzp+                  &
    &avc(ixp,iy,izpp,ib)*dtzpp+avc(ixp,iy,izm,ib)*dtzm
   d1(6)=avc(ixp,iyp,iz,ib)*dtz+avc(ixp,iyp,izp,ib)*dtzp+                &
    &avc(ixp,iyp,izpp,ib)*dtzpp+avc(ixp,iyp,izm,ib)*dtzm
   d1(7)=avc(ixp,iypp,iz,ib)*dtz+avc(ixp,iypp,izp,ib)*dtzp+              &
    &avc(ixp,iypp,izpp,ib)*dtzpp+avc(ixp,iypp,izm,ib)*dtzm
   d1(8)=avc(ixp,iym,iz,ib)*dtz+avc(ixp,iym,izp,ib)*dtzp+                &
    &avc(ixp,iym,izpp,ib)*dtzpp+avc(ixp,iym,izm,ib)*dtzm

   d1(9)=avc(ixpp,iy,iz,ib)*dtz+avc(ixpp,iy,izp,ib)*dtzp+                &
    &avc(ixpp,iy,izpp,ib)*dtzpp+avc(ixpp,iy,izm,ib)*dtzm
   d1(10)=avc(ixpp,iyp,iz,ib)*dtz+avc(ixpp,iyp,izp,ib)*dtzp+             &
    &avc(ixpp,iyp,izpp,ib)*dtzpp+avc(ixpp,iyp,izm,ib)*dtzm
   d1(11)=avc(ixpp,iypp,iz,ib)*dtz+avc(ixpp,iypp,izp,ib)*dtzp+           &
    &avc(ixpp,iypp,izpp,ib)*dtzpp+avc(ixpp,iypp,izm,ib)*dtzm
   d1(12)=avc(ixpp,iym,iz,ib)*dtz+avc(ixpp,iym,izp,ib)*dtzp+             &
    &avc(ixpp,iym,izpp,ib)*dtzpp+avc(ixpp,iym,izm,ib)*dtzm

   d1(13)=avc(ixm,iy,iz,ib)*dtz+avc(ixm,iy,izp,ib)*dtzp+                 &
    &avc(ixm,iy,izpp,ib)*dtzpp+avc(ixm,iy,izm,ib)*dtzm
   d1(14)=avc(ixm,iyp,iz,ib)*dtz+avc(ixm,iyp,izp,ib)*dtzp+               &
    &avc(ixm,iyp,izpp,ib)*dtzpp+avc(ixm,iyp,izm,ib)*dtzm
   d1(15)=avc(ixm,iypp,iz,ib)*dtz+avc(ixm,iypp,izp,ib)*dtzp+             &
    &avc(ixm,iypp,izpp,ib)*dtzpp+avc(ixm,iypp,izm,ib)*dtzm
   d1(16)=avc(ixm,iym,iz,ib)*dtz+avc(ixm,iym,izp,ib)*dtzp+               &
    &avc(ixm,iym,izpp,ib)*dtzpp+avc(ixm,iym,izm,ib)*dtzm

   grad(3,idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*  &
    &typ+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)* &
    &tym)*txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

! Go to the Cartesian grid
   grad(:,idum)=matmul(bg,grad(:,idum))

   d1(1)=avclap(ix,iy,iz,ib)*tz+avclap(ix,iy,izp,ib)*tzp+                       &
    &avclap(ix,iy,izpp,ib)*tzpp+avclap(ix,iy,izm,ib)*tzm
   d1(2)=avclap(ix,iyp,iz,ib)*tz+avclap(ix,iyp,izp,ib)*tzp+                     &
    &avclap(ix,iyp,izpp,ib)*tzpp+avclap(ix,iyp,izm,ib)*tzm
   d1(3)=avclap(ix,iypp,iz,ib)*tz+avclap(ix,iypp,izp,ib)*tzp+                   &
    &avclap(ix,iypp,izpp,ib)*tzpp+avclap(ix,iypp,izm,ib)*tzm
   d1(4)=avclap(ix,iym,iz,ib)*tz+avclap(ix,iym,izp,ib)*tzp+                     &
    &avclap(ix,iym,izpp,ib)*tzpp+avclap(ix,iym,izm,ib)*tzm
  
   d1(5)=avclap(ixp,iy,iz,ib)*tz+avclap(ixp,iy,izp,ib)*tzp                      &
    &+avclap(ixp,iy,izpp,ib)*tzpp+avclap(ixp,iy,izm,ib)*tzm
   d1(6)=avclap(ixp,iyp,iz,ib)*tz+avclap(ixp,iyp,izp,ib)*tzp                    &
    &+avclap(ixp,iyp,izpp,ib)*tzpp+avclap(ixp,iyp,izm,ib)*tzm
   d1(7)=avclap(ixp,iypp,iz,ib)*tz+avclap(ixp,iypp,izp,ib)*tzp                  &
    &+avclap(ixp,iypp,izpp,ib)*tzpp+avclap(ixp,iypp,izm,ib)*tzm
   d1(8)=avclap(ixp,iym,iz,ib)*tz+avclap(ixp,iym,izp,ib)*tzp                    &
    &+avclap(ixp,iym,izpp,ib)*tzpp+avclap(ixp,iym,izm,ib)*tzm
  
   d1(9)=avclap(ixpp,iy,iz,ib)*tz+avclap(ixpp,iy,izp,ib)*tzp                    &
    &+avclap(ixpp,iy,izpp,ib)*tzpp+avclap(ixpp,iy,izm,ib)*tzm
   d1(10)=avclap(ixpp,iyp,iz,ib)*tz+avclap(ixpp,iyp,izp,ib)*tzp                 &
    &+avclap(ixpp,iyp,izpp,ib)*tzpp+avclap(ixpp,iyp,izm,ib)*tzm
   d1(11)=avclap(ixpp,iypp,iz,ib)*tz+avclap(ixpp,iypp,izp,ib)*tzp               &
    &+avclap(ixpp,iypp,izpp,ib)*tzpp+avclap(ixpp,iypp,izm,ib)*tzm
   d1(12)=avclap(ixpp,iym,iz,ib)*tz+avclap(ixpp,iym,izp,ib)*tzp                 &
    &+avclap(ixpp,iym,izpp,ib)*tzpp+avclap(ixpp,iym,izm,ib)*tzm
  
   d1(13)=avclap(ixm,iy,iz,ib)*tz+avclap(ixm,iy,izp,ib)*tzp                     &
    &+avclap(ixm,iy,izpp,ib)*tzpp+avclap(ixm,iy,izm,ib)*tzm
   d1(14)=avclap(ixm,iyp,iz,ib)*tz+avclap(ixm,iyp,izp,ib)*tzp                   &
    &+avclap(ixm,iyp,izpp,ib)*tzpp+avclap(ixm,iyp,izm,ib)*tzm
   d1(15)=avclap(ixm,iypp,iz,ib)*tz+avclap(ixm,iypp,izp,ib)*tzp                 &
    &+avclap(ixm,iypp,izpp,ib)*tzpp+avclap(ixm,iypp,izm,ib)*tzm
   d1(16)=avclap(ixm,iym,iz,ib)*tz+avclap(ixm,iym,izp,ib)*tzp                   &
    &+avclap(ixm,iym,izpp,ib)*tzpp+avclap(ixm,iym,izm,ib)*tzm
  
!  The Laplacian
   lap(idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*typ &
    &+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*  &
    &txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

   case(2)

   d1(1)=avclap(ix,iy,iz,ib)*tz+avclap(ix,iy,izp,ib)*tzp+                       &
    &avclap(ix,iy,izpp,ib)*tzpp+avclap(ix,iy,izm,ib)*tzm
   d1(2)=avclap(ix,iyp,iz,ib)*tz+avclap(ix,iyp,izp,ib)*tzp+                     &
    &avclap(ix,iyp,izpp,ib)*tzpp+avclap(ix,iyp,izm,ib)*tzm
   d1(3)=avclap(ix,iypp,iz,ib)*tz+avclap(ix,iypp,izp,ib)*tzp+                   &
    &avclap(ix,iypp,izpp,ib)*tzpp+avclap(ix,iypp,izm,ib)*tzm
   d1(4)=avclap(ix,iym,iz,ib)*tz+avclap(ix,iym,izp,ib)*tzp+                     &
    &avclap(ix,iym,izpp,ib)*tzpp+avclap(ix,iym,izm,ib)*tzm
  
   d1(5)=avclap(ixp,iy,iz,ib)*tz+avclap(ixp,iy,izp,ib)*tzp                      &
    &+avclap(ixp,iy,izpp,ib)*tzpp+avclap(ixp,iy,izm,ib)*tzm
   d1(6)=avclap(ixp,iyp,iz,ib)*tz+avclap(ixp,iyp,izp,ib)*tzp                    &
    &+avclap(ixp,iyp,izpp,ib)*tzpp+avclap(ixp,iyp,izm,ib)*tzm
   d1(7)=avclap(ixp,iypp,iz,ib)*tz+avclap(ixp,iypp,izp,ib)*tzp                  &
    &+avclap(ixp,iypp,izpp,ib)*tzpp+avclap(ixp,iypp,izm,ib)*tzm
   d1(8)=avclap(ixp,iym,iz,ib)*tz+avclap(ixp,iym,izp,ib)*tzp                    &
    &+avclap(ixp,iym,izpp,ib)*tzpp+avclap(ixp,iym,izm,ib)*tzm
  
   d1(9)=avclap(ixpp,iy,iz,ib)*tz+avclap(ixpp,iy,izp,ib)*tzp                    &
    &+avclap(ixpp,iy,izpp,ib)*tzpp+avclap(ixpp,iy,izm,ib)*tzm
   d1(10)=avclap(ixpp,iyp,iz,ib)*tz+avclap(ixpp,iyp,izp,ib)*tzp                 &
    &+avclap(ixpp,iyp,izpp,ib)*tzpp+avclap(ixpp,iyp,izm,ib)*tzm
   d1(11)=avclap(ixpp,iypp,iz,ib)*tz+avclap(ixpp,iypp,izp,ib)*tzp               &
    &+avclap(ixpp,iypp,izpp,ib)*tzpp+avclap(ixpp,iypp,izm,ib)*tzm
   d1(12)=avclap(ixpp,iym,iz,ib)*tz+avclap(ixpp,iym,izp,ib)*tzp                 &
    &+avclap(ixpp,iym,izpp,ib)*tzpp+avclap(ixpp,iym,izm,ib)*tzm
  
   d1(13)=avclap(ixm,iy,iz,ib)*tz+avclap(ixm,iy,izp,ib)*tzp                     &
    &+avclap(ixm,iy,izpp,ib)*tzpp+avclap(ixm,iy,izm,ib)*tzm
   d1(14)=avclap(ixm,iyp,iz,ib)*tz+avclap(ixm,iyp,izp,ib)*tzp                   &
    &+avclap(ixm,iyp,izpp,ib)*tzpp+avclap(ixm,iyp,izm,ib)*tzm
   d1(15)=avclap(ixm,iypp,iz,ib)*tz+avclap(ixm,iypp,izp,ib)*tzp                 &
    &+avclap(ixm,iypp,izpp,ib)*tzpp+avclap(ixm,iypp,izm,ib)*tzm
   d1(16)=avclap(ixm,iym,iz,ib)*tz+avclap(ixm,iym,izp,ib)*tzp                   &
    &+avclap(ixm,iym,izpp,ib)*tzpp+avclap(ixm,iym,izm,ib)*tzm
  
!  The Laplacian
   lap(idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*typ &
    &+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*  &
    &txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

   d1(1)=avcgrad1(ix,iy,iz,ib)*tz+avcgrad1(ix,iy,izp,ib)*tzp+                       &
    &avcgrad1(ix,iy,izpp,ib)*tzpp+avcgrad1(ix,iy,izm,ib)*tzm
   d1(2)=avcgrad1(ix,iyp,iz,ib)*tz+avcgrad1(ix,iyp,izp,ib)*tzp+                     &
    &avcgrad1(ix,iyp,izpp,ib)*tzpp+avcgrad1(ix,iyp,izm,ib)*tzm
   d1(3)=avcgrad1(ix,iypp,iz,ib)*tz+avcgrad1(ix,iypp,izp,ib)*tzp+                   &
    &avcgrad1(ix,iypp,izpp,ib)*tzpp+avcgrad1(ix,iypp,izm,ib)*tzm
   d1(4)=avcgrad1(ix,iym,iz,ib)*tz+avcgrad1(ix,iym,izp,ib)*tzp+                     &
    &avcgrad1(ix,iym,izpp,ib)*tzpp+avcgrad1(ix,iym,izm,ib)*tzm
  
   d1(5)=avcgrad1(ixp,iy,iz,ib)*tz+avcgrad1(ixp,iy,izp,ib)*tzp                      &
    &+avcgrad1(ixp,iy,izpp,ib)*tzpp+avcgrad1(ixp,iy,izm,ib)*tzm
   d1(6)=avcgrad1(ixp,iyp,iz,ib)*tz+avcgrad1(ixp,iyp,izp,ib)*tzp                    &
    &+avcgrad1(ixp,iyp,izpp,ib)*tzpp+avcgrad1(ixp,iyp,izm,ib)*tzm
   d1(7)=avcgrad1(ixp,iypp,iz,ib)*tz+avcgrad1(ixp,iypp,izp,ib)*tzp                  &
    &+avcgrad1(ixp,iypp,izpp,ib)*tzpp+avcgrad1(ixp,iypp,izm,ib)*tzm
   d1(8)=avcgrad1(ixp,iym,iz,ib)*tz+avcgrad1(ixp,iym,izp,ib)*tzp                    &
    &+avcgrad1(ixp,iym,izpp,ib)*tzpp+avcgrad1(ixp,iym,izm,ib)*tzm
  
   d1(9)=avcgrad1(ixpp,iy,iz,ib)*tz+avcgrad1(ixpp,iy,izp,ib)*tzp                    &
    &+avcgrad1(ixpp,iy,izpp,ib)*tzpp+avcgrad1(ixpp,iy,izm,ib)*tzm
   d1(10)=avcgrad1(ixpp,iyp,iz,ib)*tz+avcgrad1(ixpp,iyp,izp,ib)*tzp                 &
    &+avcgrad1(ixpp,iyp,izpp,ib)*tzpp+avcgrad1(ixpp,iyp,izm,ib)*tzm
   d1(11)=avcgrad1(ixpp,iypp,iz,ib)*tz+avcgrad1(ixpp,iypp,izp,ib)*tzp               &
    &+avcgrad1(ixpp,iypp,izpp,ib)*tzpp+avcgrad1(ixpp,iypp,izm,ib)*tzm
   d1(12)=avcgrad1(ixpp,iym,iz,ib)*tz+avcgrad1(ixpp,iym,izp,ib)*tzp                 &
    &+avcgrad1(ixpp,iym,izpp,ib)*tzpp+avcgrad1(ixpp,iym,izm,ib)*tzm
  
   d1(13)=avcgrad1(ixm,iy,iz,ib)*tz+avcgrad1(ixm,iy,izp,ib)*tzp                     &
    &+avcgrad1(ixm,iy,izpp,ib)*tzpp+avcgrad1(ixm,iy,izm,ib)*tzm
   d1(14)=avcgrad1(ixm,iyp,iz,ib)*tz+avcgrad1(ixm,iyp,izp,ib)*tzp                   &
    &+avcgrad1(ixm,iyp,izpp,ib)*tzpp+avcgrad1(ixm,iyp,izm,ib)*tzm
   d1(15)=avcgrad1(ixm,iypp,iz,ib)*tz+avcgrad1(ixm,iypp,izp,ib)*tzp                 &
    &+avcgrad1(ixm,iypp,izpp,ib)*tzpp+avcgrad1(ixm,iypp,izm,ib)*tzm
   d1(16)=avcgrad1(ixm,iym,iz,ib)*tz+avcgrad1(ixm,iym,izp,ib)*tzp                   &
    &+avcgrad1(ixm,iym,izpp,ib)*tzpp+avcgrad1(ixm,iym,izm,ib)*tzm
  
!  The Gradient in direction of a1
   grad(1,idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*typ &
    &+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*  &
    &txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

   d1(1)=avcgrad2(ix,iy,iz,ib)*tz+avcgrad2(ix,iy,izp,ib)*tzp+                       &
    &avcgrad2(ix,iy,izpp,ib)*tzpp+avcgrad2(ix,iy,izm,ib)*tzm
   d1(2)=avcgrad2(ix,iyp,iz,ib)*tz+avcgrad2(ix,iyp,izp,ib)*tzp+                     &
    &avcgrad2(ix,iyp,izpp,ib)*tzpp+avcgrad2(ix,iyp,izm,ib)*tzm
   d1(3)=avcgrad2(ix,iypp,iz,ib)*tz+avcgrad2(ix,iypp,izp,ib)*tzp+                   &
    &avcgrad2(ix,iypp,izpp,ib)*tzpp+avcgrad2(ix,iypp,izm,ib)*tzm
   d1(4)=avcgrad2(ix,iym,iz,ib)*tz+avcgrad2(ix,iym,izp,ib)*tzp+                     &
    &avcgrad2(ix,iym,izpp,ib)*tzpp+avcgrad2(ix,iym,izm,ib)*tzm
  
   d1(5)=avcgrad2(ixp,iy,iz,ib)*tz+avcgrad2(ixp,iy,izp,ib)*tzp                      &
    &+avcgrad2(ixp,iy,izpp,ib)*tzpp+avcgrad2(ixp,iy,izm,ib)*tzm
   d1(6)=avcgrad2(ixp,iyp,iz,ib)*tz+avcgrad2(ixp,iyp,izp,ib)*tzp                    &
    &+avcgrad2(ixp,iyp,izpp,ib)*tzpp+avcgrad2(ixp,iyp,izm,ib)*tzm
   d1(7)=avcgrad2(ixp,iypp,iz,ib)*tz+avcgrad2(ixp,iypp,izp,ib)*tzp                  &
    &+avcgrad2(ixp,iypp,izpp,ib)*tzpp+avcgrad2(ixp,iypp,izm,ib)*tzm
   d1(8)=avcgrad2(ixp,iym,iz,ib)*tz+avcgrad2(ixp,iym,izp,ib)*tzp                    &
    &+avcgrad2(ixp,iym,izpp,ib)*tzpp+avcgrad2(ixp,iym,izm,ib)*tzm
  
   d1(9)=avcgrad2(ixpp,iy,iz,ib)*tz+avcgrad2(ixpp,iy,izp,ib)*tzp                    &
    &+avcgrad2(ixpp,iy,izpp,ib)*tzpp+avcgrad2(ixpp,iy,izm,ib)*tzm
   d1(10)=avcgrad2(ixpp,iyp,iz,ib)*tz+avcgrad2(ixpp,iyp,izp,ib)*tzp                 &
    &+avcgrad2(ixpp,iyp,izpp,ib)*tzpp+avcgrad2(ixpp,iyp,izm,ib)*tzm
   d1(11)=avcgrad2(ixpp,iypp,iz,ib)*tz+avcgrad2(ixpp,iypp,izp,ib)*tzp               &
    &+avcgrad2(ixpp,iypp,izpp,ib)*tzpp+avcgrad2(ixpp,iypp,izm,ib)*tzm
   d1(12)=avcgrad2(ixpp,iym,iz,ib)*tz+avcgrad2(ixpp,iym,izp,ib)*tzp                 &
    &+avcgrad2(ixpp,iym,izpp,ib)*tzpp+avcgrad2(ixpp,iym,izm,ib)*tzm
  
   d1(13)=avcgrad2(ixm,iy,iz,ib)*tz+avcgrad2(ixm,iy,izp,ib)*tzp                     &
    &+avcgrad2(ixm,iy,izpp,ib)*tzpp+avcgrad2(ixm,iy,izm,ib)*tzm
   d1(14)=avcgrad2(ixm,iyp,iz,ib)*tz+avcgrad2(ixm,iyp,izp,ib)*tzp                   &
    &+avcgrad2(ixm,iyp,izpp,ib)*tzpp+avcgrad2(ixm,iyp,izm,ib)*tzm
   d1(15)=avcgrad2(ixm,iypp,iz,ib)*tz+avcgrad2(ixm,iypp,izp,ib)*tzp                 &
    &+avcgrad2(ixm,iypp,izpp,ib)*tzpp+avcgrad2(ixm,iypp,izm,ib)*tzm
   d1(16)=avcgrad2(ixm,iym,iz,ib)*tz+avcgrad2(ixm,iym,izp,ib)*tzp                   &
    &+avcgrad2(ixm,iym,izpp,ib)*tzpp+avcgrad2(ixm,iym,izm,ib)*tzm
  
!  The Gradient in direction of a2
   grad(2,idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*typ &
    &+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*  &
    &txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

   d1(1)=avcgrad3(ix,iy,iz,ib)*tz+avcgrad3(ix,iy,izp,ib)*tzp+                       &
    &avcgrad3(ix,iy,izpp,ib)*tzpp+avcgrad3(ix,iy,izm,ib)*tzm
   d1(2)=avcgrad3(ix,iyp,iz,ib)*tz+avcgrad3(ix,iyp,izp,ib)*tzp+                     &
    &avcgrad3(ix,iyp,izpp,ib)*tzpp+avcgrad3(ix,iyp,izm,ib)*tzm
   d1(3)=avcgrad3(ix,iypp,iz,ib)*tz+avcgrad3(ix,iypp,izp,ib)*tzp+                   &
    &avcgrad3(ix,iypp,izpp,ib)*tzpp+avcgrad3(ix,iypp,izm,ib)*tzm
   d1(4)=avcgrad3(ix,iym,iz,ib)*tz+avcgrad3(ix,iym,izp,ib)*tzp+                     &
    &avcgrad3(ix,iym,izpp,ib)*tzpp+avcgrad3(ix,iym,izm,ib)*tzm
  
   d1(5)=avcgrad3(ixp,iy,iz,ib)*tz+avcgrad3(ixp,iy,izp,ib)*tzp                      &
    &+avcgrad3(ixp,iy,izpp,ib)*tzpp+avcgrad3(ixp,iy,izm,ib)*tzm
   d1(6)=avcgrad3(ixp,iyp,iz,ib)*tz+avcgrad3(ixp,iyp,izp,ib)*tzp                    &
    &+avcgrad3(ixp,iyp,izpp,ib)*tzpp+avcgrad3(ixp,iyp,izm,ib)*tzm
   d1(7)=avcgrad3(ixp,iypp,iz,ib)*tz+avcgrad3(ixp,iypp,izp,ib)*tzp                  &
    &+avcgrad3(ixp,iypp,izpp,ib)*tzpp+avcgrad3(ixp,iypp,izm,ib)*tzm
   d1(8)=avcgrad3(ixp,iym,iz,ib)*tz+avcgrad3(ixp,iym,izp,ib)*tzp                    &
    &+avcgrad3(ixp,iym,izpp,ib)*tzpp+avcgrad3(ixp,iym,izm,ib)*tzm
  
   d1(9)=avcgrad3(ixpp,iy,iz,ib)*tz+avcgrad3(ixpp,iy,izp,ib)*tzp                    &
    &+avcgrad3(ixpp,iy,izpp,ib)*tzpp+avcgrad3(ixpp,iy,izm,ib)*tzm
   d1(10)=avcgrad3(ixpp,iyp,iz,ib)*tz+avcgrad3(ixpp,iyp,izp,ib)*tzp                 &
    &+avcgrad3(ixpp,iyp,izpp,ib)*tzpp+avcgrad3(ixpp,iyp,izm,ib)*tzm
   d1(11)=avcgrad3(ixpp,iypp,iz,ib)*tz+avcgrad3(ixpp,iypp,izp,ib)*tzp               &
    &+avcgrad3(ixpp,iypp,izpp,ib)*tzpp+avcgrad3(ixpp,iypp,izm,ib)*tzm
   d1(12)=avcgrad3(ixpp,iym,iz,ib)*tz+avcgrad3(ixpp,iym,izp,ib)*tzp                 &
    &+avcgrad3(ixpp,iym,izpp,ib)*tzpp+avcgrad3(ixpp,iym,izm,ib)*tzm
  
   d1(13)=avcgrad3(ixm,iy,iz,ib)*tz+avcgrad3(ixm,iy,izp,ib)*tzp                     &
    &+avcgrad3(ixm,iy,izpp,ib)*tzpp+avcgrad3(ixm,iy,izm,ib)*tzm
   d1(14)=avcgrad3(ixm,iyp,iz,ib)*tz+avcgrad3(ixm,iyp,izp,ib)*tzp                   &
    &+avcgrad3(ixm,iyp,izpp,ib)*tzpp+avcgrad3(ixm,iyp,izm,ib)*tzm
   d1(15)=avcgrad3(ixm,iypp,iz,ib)*tz+avcgrad3(ixm,iypp,izp,ib)*tzp                 &
    &+avcgrad3(ixm,iypp,izpp,ib)*tzpp+avcgrad3(ixm,iypp,izm,ib)*tzm
   d1(16)=avcgrad3(ixm,iym,iz,ib)*tz+avcgrad3(ixm,iym,izp,ib)*tzp                   &
    &+avcgrad3(ixm,iym,izpp,ib)*tzpp+avcgrad3(ixm,iym,izm,ib)*tzm
  
!  The Gradient in direction of a3
   grad(3,idum)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*typ &
    &+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*  &
    &txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

   end select

  endif
 enddo

 END SUBROUTINE blip3dgamma_gl

! Called by sum_orbs_over_grid_with_blips to check if orbs are pure real or imag
 SUBROUTINE blip_one_band(rpsi,grad,lap,r,avc,nrbwf,bg,igl)

 IMPLICIT NONE

 INTEGER,INTENT(in) :: igl
 INTEGER,INTENT(in) :: nrbwf(3)
 REAL(dp),INTENT(in) :: r(3),bg(3,3)
 COMPLEX(dp),INTENT(in) :: avc(0:nrbwf(1)-1,0:nrbwf(2)-1,0:nrbwf(3)-1)

 COMPLEX(dp),INTENT(out) :: rpsi,grad(3),lap
 
 INTEGER :: ix,iy,iz,ixp,ixpp,ixm,iyp,iypp,iym,izp,izpp,izm
 REAL(dp) :: t(3),txm,tx,txp,txpp,tym,ty,typ,typp,tzm,tz,tzp,tzpp,  &
  &x,y,z,b(6)
 REAL(dp) :: dtxm=0.d0,d2txm=0.d0,dtx=0.d0,d2tx=0.d0,dtxp=0.d0,d2txp=0.d0,     &
  &dtxpp=0.d0,d2txpp=0.d0,dtym=0.d0,d2tym=0.d0,dty=0.d0,d2ty=0.d0,dtyp=0.d0,   &
  &d2typ=0.d0,dtypp=0.d0,d2typp=0.d0,dtzm=0.d0,d2tzm=0.d0,dtz=0.d0,d2tz=0.d0,  &
  &dtzp=0.d0,d2tzp=0.d0,dtzpp=0.d0,d2tzpp=0.d0
 COMPLEX(dp) :: d1(16)

 ix=int(mod(r(1)+abs(int(r(1)))+1,1.d0)*nrbwf(1))
 iy=int(mod(r(2)+abs(int(r(2)))+1,1.d0)*nrbwf(2))
 iz=int(mod(r(3)+abs(int(r(3)))+1,1.d0)*nrbwf(3))

 if(ix<0.or.iy<0.or.iz<0)call errstop('BLIP3D','Negative index found.')

! The blips are defined as the product of one-dimensional cubic splines, These
! are different from zero only on four adjacent grid points.  It follows that
! for any vector r there are only 64 overlapping three-dimensional blips, which
! are the product of all possible combinations of the three one-dimensional
! splines.

! These are the extra 3 coefficients for each dimension needed
 ix=mod(ix,nrbwf(1))
 iy=mod(iy,nrbwf(2))
 iz=mod(iz,nrbwf(3))
 ixp=mod(ix+1,nrbwf(1))
 ixpp=mod(ix+2,nrbwf(1))
 ixm=mod(ix-1+nrbwf(1),nrbwf(1))
 iyp=mod(iy+1,nrbwf(2))
 iypp=mod(iy+2,nrbwf(2))
 iym=mod(iy-1+nrbwf(2),nrbwf(2))
 izp=mod(iz+1,nrbwf(3))
 izpp=mod(iz+2,nrbwf(3))
 izm=mod(iz-1+nrbwf(3),nrbwf(3))

! Now calculate the 12 monodimensional blip functions
 t=mod(r+abs(int(r))+1,1.d0)*nrbwf
 t=mod(t,dble(nrbwf))

 x=t(1)-ix+1.d0
 txm=2.d0-3*x+1.5d0*x*x-0.25d0*x*x*x
 if(igl==1)then
  dtxm=(-3.d0+3*x-0.75d0*x*x)*nrbwf(1)
  d2txm=(3.d0-1.5d0*x)*nrbwf(1)*nrbwf(1)
 endif
 x=t(1)-ix
 tx=1.d0-1.5d0*x*x+0.75d0*x*x*x
 if(igl==1)then
  dtx=(-3.d0*x+2.25d0*x*x)*nrbwf(1)
  d2tx=(-3.d0+4.5d0*x)*nrbwf(1)*nrbwf(1)
 endif
 x=t(1)-ix-1.d0
 txp=1.d0-1.5d0*x*x-0.75d0*x*x*x
 if(igl==1)then
  dtxp=(-3.d0*x-2.25d0*x*x)*nrbwf(1)
  d2txp=(-3.d0-4.5d0*x)*nrbwf(1)*nrbwf(1)
 endif
 x=t(1)-ix-2.d0
 txpp=2.d0+3*x+1.5d0*x*x+0.25d0*x*x*x
 if(igl==1)then
  dtxpp=(3.d0+3*x+0.75d0*x*x)*nrbwf(1)
  d2txpp=(3.d0+1.5d0*x)*nrbwf(1)*nrbwf(1)
 endif

 y=t(2)-iy+1.d0
 tym=2.d0-3*y+1.5d0*y*y-0.25d0*y*y*y
 if(igl==1)then
  dtym=(-3.d0+3*y-0.75d0*y*y)*nrbwf(2)
  d2tym=(3.d0-1.5d0*y)*nrbwf(2)*nrbwf(2)
 endif
 y=t(2)-iy
 ty=1.d0-1.5d0*y*y+0.75d0*y*y*y
 if(igl==1)then
  dty=(-3.d0*y+2.25d0*y*y)*nrbwf(2)
  d2ty=(-3.d0+4.5d0*y)*nrbwf(2)*nrbwf(2)
 endif
 y=t(2)-iy-1.d0
 typ=1.d0-1.5d0*y*y-0.75d0*y*y*y
 if(igl==1)then
  dtyp=(-3.d0*y-2.25d0*y*y)*nrbwf(2)
  d2typ=(-3.d0-4.5d0*y)*nrbwf(2)*nrbwf(2)
 endif
 y=t(2)-iy-2.d0
 typp=2.d0+3*y+1.5d0*y*y+0.25d0*y*y*y
 if(igl==1) then
  dtypp=(3.d0+3*y+0.75d0*y*y)*nrbwf(2)
  d2typp=(3.d0+1.5d0*y)*nrbwf(2)*nrbwf(2)
 endif

 z=t(3)-iz+1.d0
 tzm=2.d0-3*z+1.5d0*z*z-0.25d0*z*z*z
 if(igl==1)then
  dtzm=(-3.d0+3*z-0.75d0*z*z)*nrbwf(3)
  d2tzm=(3.d0-1.5d0*z)*nrbwf(3)*nrbwf(3)
 endif
 z=t(3)-iz
 tz=1.d0-1.5d0*z*z+0.75d0*z*z*z
 if(igl==1)then
  dtz=(-3.d0*z+2.25d0*z*z)*nrbwf(3)
  d2tz=(-3.d0+4.5d0*z)*nrbwf(3)*nrbwf(3)
 endif
 z=t(3)-iz-1.d0
 tzp=1.d0-1.5d0*z*z-0.75d0*z*z*z
 if(igl==1)then
  dtzp=(-3.d0*z-2.25d0*z*z)*nrbwf(3)
  d2tzp=(-3.d0-4.5d0*z)*nrbwf(3)*nrbwf(3)
 endif
 z=t(3)-iz-2.d0
 tzpp=2.d0+3*z+1.5d0*z*z+0.25d0*z*z*z
 if(igl==1)then
  dtzpp=(3.d0+3*z+0.75d0*z*z)*nrbwf(3)
  d2tzpp=(3.d0+1.5d0*z)*nrbwf(3)*nrbwf(3)
 endif

 b(1)=(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2)
 b(2)=(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2)
 b(3)=(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2)
 b(4)=2*(bg(1,1)*bg(1,2)+bg(2,1)*bg(2,2)+bg(3,1)*bg(3,2))
 b(5)=2*(bg(1,2)*bg(1,3)+bg(2,2)*bg(2,3)+bg(3,2)*bg(3,3))
 b(6)=2*(bg(1,3)*bg(1,1)+bg(2,3)*bg(2,1)+bg(3,3)*bg(3,1))


 d1(1)=avc(ix,iy,iz)*tz+avc(ix,iy,izp)*tzp+                     &
      avc(ix,iy,izpp)*tzpp+avc(ix,iy,izm)*tzm
 d1(2)=avc(ix,iyp,iz)*tz+avc(ix,iyp,izp)*tzp+                   &
      avc(ix,iyp,izpp)*tzpp+avc(ix,iyp,izm)*tzm
 d1(3)=avc(ix,iypp,iz)*tz+avc(ix,iypp,izp)*tzp+                 &
      avc(ix,iypp,izpp)*tzpp+avc(ix,iypp,izm)*tzm
 d1(4)=avc(ix,iym,iz)*tz+avc(ix,iym,izp)*tzp+                   &
      avc(ix,iym,izpp)*tzpp+avc(ix,iym,izm)*tzm

 d1(5)=avc(ixp,iy,iz)*tz+avc(ixp,iy,izp)*tzp                    &
      +avc(ixp,iy,izpp)*tzpp+avc(ixp,iy,izm)*tzm
 d1(6)=avc(ixp,iyp,iz)*tz+avc(ixp,iyp,izp)*tzp                  &
      +avc(ixp,iyp,izpp)*tzpp+avc(ixp,iyp,izm)*tzm
 d1(7)=avc(ixp,iypp,iz)*tz+avc(ixp,iypp,izp)*tzp                &
      +avc(ixp,iypp,izpp)*tzpp+avc(ixp,iypp,izm)*tzm
 d1(8)=avc(ixp,iym,iz)*tz+avc(ixp,iym,izp)*tzp                  &
      +avc(ixp,iym,izpp)*tzpp+avc(ixp,iym,izm)*tzm

 d1(9)=avc(ixpp,iy,iz)*tz+avc(ixpp,iy,izp)*tzp                  &
      +avc(ixpp,iy,izpp)*tzpp+avc(ixpp,iy,izm)*tzm
 d1(10)=avc(ixpp,iyp,iz)*tz+avc(ixpp,iyp,izp)*tzp               &
      +avc(ixpp,iyp,izpp)*tzpp+avc(ixpp,iyp,izm)*tzm
 d1(11)=avc(ixpp,iypp,iz)*tz+avc(ixpp,iypp,izp)*tzp             &
      +avc(ixpp,iypp,izpp)*tzpp+avc(ixpp,iypp,izm)*tzm
 d1(12)=avc(ixpp,iym,iz)*tz+avc(ixpp,iym,izp)*tzp               &
      +avc(ixpp,iym,izpp)*tzpp+avc(ixpp,iym,izm)*tzm

 d1(13)=avc(ixm,iy,iz)*tz+avc(ixm,iy,izp)*tzp                   &
      +avc(ixm,iy,izpp)*tzpp+avc(ixm,iy,izm)*tzm
 d1(14)=avc(ixm,iyp,iz)*tz+avc(ixm,iyp,izp)*tzp                 &
      +avc(ixm,iyp,izpp)*tzpp+avc(ixm,iyp,izm)*tzm
 d1(15)=avc(ixm,iypp,iz)*tz+avc(ixm,iypp,izp)*tzp               &
      +avc(ixm,iypp,izpp)*tzpp+avc(ixm,iypp,izm)*tzm
 d1(16)=avc(ixm,iym,iz)*tz+avc(ixm,iym,izp)*tzp                 &
    +avc(ixm,iym,izpp)*tzpp+avc(ixm,iym,izm)*tzm

! The function
 rpsi=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*typ &
      +d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)* &
      txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm
 if(igl == 1) then
! The gradient, first term involving Theta'(x)
     grad(1)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+(d1(5)*ty+d1(6)*    &
          typ+d1(7)*typp+d1(8)*tym)*dtxp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+     &
          d1(12)*tym)*dtxpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm

! Second term of the gradient involving Theta'(y)
     grad(2)=(d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+(d1(5)*dty+      &
          d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+(d1(9)*dty+d1(10)*dtyp+       &
          d1(11)*dtypp+d1(12)*dtym)*txpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+ &
          d1(16)*dtym)*txm

! Third term of the gradient involving Theta'(z)
     grad(3)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*     &
          typ+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+      &
          d1(12)*tym)*txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

! The Laplacian: first term involving Theta''(x)
     lap=((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*d2tx+(d1(5)*ty+d1(6)*typ+ &
          d1(7)*typp+d1(8)*tym)*d2txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)* &
          tym)*d2txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*d2txm)*b(1)

! The Laplacian: term involving Theta'(x)Theta'(y)
     lap=lap+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*dtx+(d1(5)*dty+  &
          d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*dtxp+(d1(9)*dty+d1(10)*dtyp+      &
          d1(11)*dtypp+d1(12)*dtym)*dtxpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*      &
          dtypp+d1(16)*dtym)*dtxm)*b(4)

! Second term of the laplacian involving Theta''(y)
     lap=lap+((d1(1)*d2ty+d1(2)*d2typ+d1(3)*d2typp+d1(4)*d2tym)*tx+          &
          (d1(5)*d2ty+d1(6)*d2typ+d1(7)*d2typp+d1(8)*d2tym)*txp+(d1(9)*d2ty+   &
          d1(10)*d2typ+d1(11)*d2typp+d1(12)*d2tym)*txpp+(d1(13)*d2ty+d1(14)*   &
          d2typ+d1(15)*d2typp+d1(16)*d2tym)*txm)*b(2)

! And now the third term of the laplacian involving Theta''(z)
     d1(1)=avc(ix,iy,iz)*d2tz+avc(ix,iy,izp)*d2tzp+                &
          avc(ix,iy,izpp)*d2tzpp+avc(ix,iy,izm)*d2tzm
     d1(2)=avc(ix,iyp,iz)*d2tz+avc(ix,iyp,izp)*d2tzp+              &
          avc(ix,iyp,izpp)*d2tzpp+avc(ix,iyp,izm)*d2tzm
     d1(3)=avc(ix,iypp,iz)*d2tz+avc(ix,iypp,izp)*d2tzp+            &
          avc(ix,iypp,izpp)*d2tzpp+avc(ix,iypp,izm)*d2tzm
     d1(4)=avc(ix,iym,iz)*d2tz+avc(ix,iym,izp)*d2tzp+              &
          avc(ix,iym,izpp)*d2tzpp+avc(ix,iym,izm)*d2tzm

     d1(5)=avc(ixp,iy,iz)*d2tz+avc(ixp,iy,izp)*d2tzp+              &
          avc(ixp,iy,izpp)*d2tzpp+avc(ixp,iy,izm)*d2tzm
     d1(6)=avc(ixp,iyp,iz)*d2tz+avc(ixp,iyp,izp)*d2tzp+            &
          avc(ixp,iyp,izpp)*d2tzpp+avc(ixp,iyp,izm)*d2tzm
     d1(7)=avc(ixp,iypp,iz)*d2tz+avc(ixp,iypp,izp)*d2tzp+          &
          avc(ixp,iypp,izpp)*d2tzpp+avc(ixp, iypp, izm)*d2tzm
     d1(8)=avc(ixp,iym,iz)*d2tz+avc(ixp,iym,izp)*d2tzp+            &
          avc(ixp,iym,izpp)*d2tzpp+avc(ixp,iym,izm)*d2tzm

     d1(9)=avc(ixpp,iy,iz)*d2tz+avc(ixpp,iy,izp)*d2tzp+            &
          avc(ixpp,iy,izpp)*d2tzpp+avc(ixpp,iy,izm)*d2tzm
     d1(10)=avc(ixpp,iyp,iz)*d2tz+avc(ixpp,iyp,izp)*d2tzp+         &
          avc(ixpp,iyp,izpp)*d2tzpp+avc(ixpp,iyp,izm)*d2tzm
     d1(11)=avc(ixpp,iypp,iz)*d2tz+avc(ixpp,iypp,izp)*d2tzp+       &
          avc(ixpp,iypp,izpp)*d2tzpp+avc(ixpp,iypp,izm)*d2tzm
     d1(12)=avc(ixpp,iym,iz)*d2tz+avc(ixpp,iym,izp)*d2tzp+         &
          avc(ixpp,iym,izpp)*d2tzpp+avc(ixpp,iym,izm)*d2tzm

     d1(13)=avc(ixm,iy,iz)*d2tz+avc(ixm,iy,izp)*d2tzp+             &
          avc(ixm,iy,izpp)*d2tzpp+avc(ixm,iy,izm)*d2tzm
     d1(14)=avc(ixm,iyp,iz)*d2tz+avc(ixm,iyp,izp)*d2tzp+           &
          avc(ixm,iyp,izpp)*d2tzpp+avc(ixm,iyp,izm)*d2tzm
     d1(15)=avc(ixm,iypp,iz)*d2tz+avc(ixm,iypp,izp)*d2tzp+         &
          avc(ixm,iypp,izpp)*d2tzpp+avc(ixm,iypp,izm)*d2tzm
     d1(16)=avc(ixm,iym,iz)*d2tz+avc(ixm,iym,izp)*d2tzp+           &
          avc(ixm,iym,izpp)*d2tzpp+avc(ixm,iym,izm)*d2tzm

! Theta''(z)
     lap=lap+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+(d1(5)*ty+d1(6)*  &
          typ+d1(7)*typp+d1(8)*tym)*txp+(d1(9)*ty+d1(10)*typ+d1(11)*typp+      &
          d1(12)* tym)*txpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm) &
          *b(3)

! And the third term of the gradient involving Theta'(z)
     d1(1)=avc(ix,iy,iz)*dtz+avc(ix,iy,izp)*dtzp+                  &
          avc(ix,iy,izpp)*dtzpp+avc(ix, iy, izm)*dtzm
     d1(2)=avc(ix,iyp,iz)*dtz+avc(ix,iyp,izp)*dtzp+                &
          avc(ix,iyp,izpp)*dtzpp+avc(ix,iyp,izm)*dtzm
     d1(3)=avc(ix,iypp,iz)*dtz+avc(ix,iypp,izp)*dtzp+              &
          avc(ix,iypp,izpp)*dtzpp+avc(ix,iypp,izm)*dtzm
     d1(4)=avc(ix,iym,iz)*dtz+avc(ix,iym,izp)*dtzp+                &
          avc(ix,iym,izpp)*dtzpp+avc(ix,iym,izm)*dtzm

     d1(5)=avc(ixp,iy,iz)*dtz+avc(ixp,iy,izp)*dtzp+                &
          avc(ixp,iy,izpp)*dtzpp+avc(ixp,iy,izm)*dtzm
     d1(6)=avc(ixp,iyp,iz)*dtz+avc(ixp,iyp,izp)*dtzp+              &
          avc(ixp,iyp,izpp)*dtzpp+avc(ixp,iyp,izm)*dtzm
     d1(7)=avc(ixp,iypp,iz)*dtz+avc(ixp,iypp,izp)*dtzp+            &
          avc(ixp,iypp,izpp)*dtzpp+avc(ixp,iypp,izm)*dtzm
     d1(8)=avc(ixp,iym,iz)*dtz+avc(ixp,iym,izp)*dtzp+              &
          avc(ixp,iym,izpp)*dtzpp+avc(ixp,iym,izm)*dtzm

     d1(9)=avc(ixpp,iy,iz)*dtz+avc(ixpp,iy,izp)*dtzp+              &
          avc(ixpp,iy,izpp)*dtzpp+avc(ixpp,iy,izm)*dtzm
     d1(10)=avc(ixpp,iyp,iz)*dtz+avc(ixpp,iyp,izp)*dtzp+           &
          avc(ixpp,iyp,izpp)*dtzpp+avc(ixpp,iyp,izm)*dtzm
     d1(11)=avc(ixpp,iypp,iz)*dtz+avc(ixpp,iypp,izp)*dtzp+         &
          avc(ixpp,iypp,izpp)*dtzpp+avc(ixpp,iypp,izm)*dtzm
     d1(12)=avc(ixpp,iym,iz)*dtz+avc(ixpp,iym,izp)*dtzp+           &
          avc(ixpp,iym,izpp)*dtzpp+avc(ixpp,iym,izm)*dtzm

     d1(13)=avc(ixm,iy,iz)*dtz+avc(ixm,iy,izp)*dtzp+               &
          avc(ixm,iy,izpp)*dtzpp+avc(ixm,iy,izm)*dtzm
     d1(14)=avc(ixm,iyp,iz)*dtz+avc(ixm,iyp,izp)*dtzp+             &
          avc(ixm,iyp,izpp)*dtzpp+avc(ixm,iyp,izm)*dtzm
     d1(15)=avc(ixm,iypp,iz)*dtz+avc(ixm,iypp,izp)*dtzp+           &
          avc(ixm,iypp,izpp)*dtzpp+avc(ixm,iypp,izm)*dtzm
     d1(16)=avc(ixm,iym,iz)*dtz+avc(ixm,iym,izp)*dtzp+             &
          avc(ixm,iym,izpp)*dtzpp+avc(ixm,iym,izm)*dtzm

! The Laplacian: term involving Theta'(x)Theta'(z)
     lap=lap+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+(d1(5)*ty+       &
          d1(6)*typ+d1(7)*typp+d1(8)*tym)*dtxp+(d1(9)*ty+d1(10)*typ+d1(11)*    &
          typp+d1(12)*tym)*dtxpp+(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym) &
          *dtxm)*b(6)

! the Laplacian: term involving Theta'(y)Theta'(z)
     lap=lap+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+(d1(5)        &
          *dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+(d1(9)*dty+d1(10)*dtyp+  &
          d1(11)*dtypp+d1(12)*dtym)*txpp+(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+ &
          d1(16)*dtym)*txm)*b(5)

     grad=matmul(bg,grad)
  end if !il == 1

 END SUBROUTINE blip_one_band


 SUBROUTINE sum_orbs_over_grid_with_blips(r,sum_orbs,ispin)
! This routine called repeatedly by bwfdet_setup to accumulate the squared
! real and imaginary parts of all orbitals summed over points on a grid.
! in order to check whether BZ edge orbitals are pure real, pure imaginary, or
! neither.
 IMPLICIT NONE
 INTEGER,INTENT(in) :: ispin
 REAL(dp),INTENT(in) :: r(3)
 COMPLEX(dp),INTENT(inout) :: sum_orbs(maxband,nkvec_bwfdet)
 INTEGER band,k
 REAL(dp) kdotrprime,r1(3)
 COMPLEX(dp) tot,totg(3),totl,ekdr
 LOGICAL spin1

 spin1=.true. ; if(ispin==2.and.spin_polarized)spin1=.false.

 r1=matmul(transpose(painv),r)
 do k=1,nkvec_bwfdet
  if(lkedge(k))then
   kdotrprime=kvec_bwfdet(1,k)*r(1)+kvec_bwfdet(2,k)*r(2)+kvec_bwfdet(3,k)*r(3)
   ekdr=exp(kdotrprime*zi)
   do band=1,boccband(k,ispin)
    tot=0.d0
    if(spin1)then ! spin 1 (and spin 2 in spin restricted case)
     if(.not.pwreal)call blip_one_band(tot,totg,totl,r1,cavc(0,0,0,band,k),nrbwf,painv,0)
    else ! spin 2 in spin unrbwfestricted case
     if(.not.pwreal)call blip_one_band(tot,totg,totl,r1,cavc2(0,0,0,band,k),nrbwf,painv,0)
    endif
    tot=ekdr*tot ! = u(r) * exp[ik.r]
! 'Accumulate' the squared components of the orbital
    sum_orbs(band,k)=sum_orbs(band,k)+(dble(tot)**2)+zi*(aimag(tot)**2)
   enddo ! band
  endif
 enddo ! k
 END SUBROUTINE sum_orbs_over_grid_with_blips

 SUBROUTINE get_real_and_imaginary_parts_of_blips(r,rpsi,lap,ispin)
! This routine called repeatedly by bwfdet_setup to accumulate the squared
! real and imaginary parts of all orbitals summed over points on a grid.
! in order to check whether BZ edge orbitals are pure real, pure imaginary, or
! neither.
 IMPLICIT NONE
 INTEGER,INTENT(in) :: ispin
 REAL(dp),INTENT(in) :: r(3)
 COMPLEX(dp),INTENT(out) :: rpsi(nemax,nkvec_bwfdet),lap(nemax,nkvec_bwfdet)
 INTEGER band,k
 REAL(dp) kdotrprime,r1(3)
 COMPLEX(dp) tot,totg(3),totl,ekdr,kekdr(3),k2ekdr
 LOGICAL spin1

 spin1=.true. ; if(ispin==2.and.spin_polarized)spin1=.false.

 r1=matmul(transpose(painv),r)
 do k=1,nkvec_bwfdet
  if(lkedge(k))then
   kdotrprime=kvec_bwfdet(1,k)*r(1)+kvec_bwfdet(2,k)*r(2)+kvec_bwfdet(3,k)*r(3)
   ekdr=exp(kdotrprime*zi)
   kekdr=zi*kvec_bwfdet(:,k)*ekdr
   k2ekdr=-dot_product(kvec_bwfdet(:,k),kvec_bwfdet(:,k))*ekdr
   do band=1,boccband(k,ispin)
    tot=0.d0
    if(spin1)then ! spin 1 (and spin 2 in spin restricted case)
     if(.not.pwreal)call blip_one_band(tot,totg,totl,r1,cavc(0,0,0,band,k),nrbwf,painv,1)
    else ! spin 2 in spin unrbwfestricted case
     if(.not.pwreal)call blip_one_band(tot,totg,totl,r1,cavc2(0,0,0,band,k),nrbwf,painv,1)
    endif
    rpsi(band,k)=ekdr*tot         ! = u(r) * exp[ik.r]
    lap(band,k)=ekdr*totl  &      ! = Lap[u(r)] * exp[ik.r] 
    &+2*(kekdr(1)*totg(1)+ &      !+ 2*Grad[u(r)].k*exp[ik.r] 
    &    kekdr(2)*totg(2)+ &
    &    kekdr(3)*totg(3)) &
    &+k2ekdr*tot                   !+ u(r) * -k.k exp[ik.r]
   enddo ! band
  endif
 enddo ! k
 END SUBROUTINE get_real_and_imaginary_parts_of_blips


 SUBROUTINE bwfdet_wrapper(rvec,iw,igl,ispin,rpsi,grad,lap)
! Not used in CHAMP
! Wrapper to improve behaviour of the Hitachi SR2201  - now sadly defunct :-)
 IMPLICIT NONE
 INTEGER igl,iw,ispin
 REAL(dp) grad(3,nemax,ndet_bwfdet),lap(nemax,ndet_bwfdet),rpsi(nemax,ndet_bwfdet),rvec(3)
 call bwfdet_main(rvec,iw,igl,ispin,rpsi,grad,lap)
 END SUBROUTINE bwfdet_wrapper

 SUBROUTINE indexx(arr,index)
!------------------------------------------------------------!
! Simple sorting routine based on Numerical Recipes version. !
! Creates index table for the array arr.                     !
!------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,DIMENSION(:),INTENT(out) :: index
 REAL(dp),DIMENSION(:),INTENT(in) :: arr
 INTEGER n,k,i,j,indext,jstack,l,r,n_size_index,n_size_arr
 INTEGER,PARAMETER :: nn=15,nstack=50,npar_arth=16,npar2_arth=8
 INTEGER,DIMENSION(nstack) :: istack
 REAL(dp) a
 n_size_index=size(index)
 n_size_arr=size(arr)
 n=assert_eq(n_size_index,n_size_arr)
 index=arth(1,1,n)
 jstack=0
 l=1
 r=n
 do
  if(r-l<nn)then
   do j=l+1,r
    indext=index(j)
    a=arr(indext)
    do i=j-1,l,-1
     if(arr(index(i))<=a)exit
     index(i+1)=index(i)
    enddo
    index(i+1)=indext
   enddo
   if(jstack==0)return
   r=istack(jstack)
   l=istack(jstack-1)
   jstack=jstack-2
  else
   k=(l+r)/2
   call swap(index(k),index(l+1))
   call icomp_xchg(index(l),index(r))
   call icomp_xchg(index(l+1),index(r))
   call icomp_xchg(index(l),index(l+1))
   i=l+1
   j=r
   indext=index(l+1)
   a=arr(indext)
   do
    do
     i=i+1
     if(arr(index(i))>=a)exit
    enddo
    do
     j=j-1
     if(arr(index(j))<=a)exit
    enddo
    if(j<i)exit
    call swap(index(i),index(j))
   enddo
   index(l+1)=index(j)
   index(j)=indext
   jstack=jstack+2
   if(jstack>nstack)call errstop('INDEXX','Stack too small - increase nstack&
    & parameter.')
   if(r-i+1>=j-l)then
    istack(jstack)=r
    istack(jstack-1)=i
    r=j-1
   else
    istack(jstack)=j-1
    istack(jstack-1)=l
    l=i
   endif
  endif
 enddo


 CONTAINS


  SUBROUTINE icomp_xchg(i,j)
  INTEGER,INTENT(inout) :: i,j
  INTEGER swp
  if(arr(j)<arr(i))then
   swp=i
   i=j
   j=swp
  endif
  END SUBROUTINE icomp_xchg


  FUNCTION arth(first,increment,n)
  INTEGER,INTENT(in) :: first,increment,n
  INTEGER,DIMENSION(n) :: arth
  INTEGER k,k2,temp
  if(n>0)arth(1)=first
  if(n<=npar_arth)then
   do k=2,n
    arth(k)=arth(k-1)+increment
   enddo
  else
   do k=2,npar2_arth
    arth(k)=arth(k-1)+increment
   enddo
   temp=increment*npar2_arth
   k=npar2_arth
   do
    if(k>=n)exit
    k2=k+k
    arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
    temp=temp+temp
    k=k2
   enddo
  endif
  END FUNCTION arth


  FUNCTION assert_eq(n1,n2)
  INTEGER,INTENT(in) :: n1,n2
  INTEGER assert_eq
  if(n1==n2)then
   assert_eq=n1
  else
   assert_eq=0
   call errstop('INDEXX','An assert_eq failed.')
  endif
  END FUNCTION assert_eq


  SUBROUTINE swap(a,b)
  INTEGER,INTENT(inout) :: a,b
  INTEGER dum
  dum=a
  a=b
  b=dum
  END SUBROUTINE swap


 END SUBROUTINE indexx

 SUBROUTINE open_units(io_no,ierr)
!-----------------------------------------------------------!
! Find a unit number for i/o that isn't already being used. !
!-----------------------------------------------------------!
 IMPLICIT NONE
 INTEGER io_no,ierr
 ierr=0
 do io_no=10,99
  if(.not.open_unit(io_no))exit
 enddo
 open_unit(io_no)=.true.
 if(io_no==99)ierr=1
 END SUBROUTINE open_units

 SUBROUTINE errstop(subroutine,error)
!-------------------------------------------------------!
! Write out routine name and error message then stop.   !
!-------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: subroutine,error
 write(6,1)subroutine,error
1 format(/1x,'ERROR : ',a,/1x,a/)
 write(6,*)
 write(6,'(1x,78(''-''))')
 write(6,*)
 call qmc_abort
 END SUBROUTINE errstop

 SUBROUTINE errwarn(subroutine,warning)
!--------------------------!
! Print a warning message. !
!--------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: subroutine,warning
 if(idtask == 0)then
  write(6,'(1x,78(''*''))')
  write(6,'(/1x,''WARNING : '',a,/1x,a/)')subroutine,warning
  write(6,'(1x,78(''*''),/)')
 endif
 END SUBROUTINE errwarn

 SUBROUTINE skip(iunit,nskip)
!---------------------------------------!
! Skip records in a free format file.   !
!---------------------------------------!
 IMPLICIT NONE
 INTEGER iunit,nskip,iskip
 do iskip=1,nskip
  read(iunit,fmt=*)
 enddo
 END SUBROUTINE skip

 SUBROUTINE qmc_barrier
 IMPLICIT NONE

 END SUBROUTINE qmc_barrier

 SUBROUTINE qmc_abort
 IMPLICIT NONE

 stop
 END SUBROUTINE qmc_abort


 CHARACTER(20) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! USE utilities                                                         !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
  INTEGER i,j,n
  CHARACTER tmp,sign

  if(n==0)then
   i2s='0' ; return
  endif
  sign=' ' ; if(n<0)sign='-'

  do i=1,len(i2s)
   i2s(i:i)=' '
  enddo

  i=abs(n)
  do j=1,len(i2s)
   if(i==0)exit
   i2s(j:j)=achar(ichar('0')+mod(i,10))
   i=i/10
  enddo

  i=1 ; j=len_trim(i2s)
  do
   if(i>=j)exit
   tmp=i2s(j:j)
   i2s(j:j)=i2s(i:i)
   i2s(i:i)=tmp
   i=i+1
   j=j-1
  enddo

  i2s=trim(sign)//i2s

 END FUNCTION i2s

 CHARACTER(r2s_length) FUNCTION r2s(r,real_format)
!-------------------------------------------------------------------------!
! Converts real variable with arbitrary format to string that can be      !
! trimmed and printed in the middle of a sentence without introducing     !
! large amounts of white space, as you would if you did                   !
! write(6,'(f12.6)')12.0 or similar. Note you need to pass through the    !
! format string e.g. f12.6 .                                              !
!                                                                         !
! Calling routine is intended to include something like:                  !
! USE utilities                                                           !
! REAL(dp) r                                                              !
! r=12.d0                                                                 !
! tmpr=r2s(r,'(f12.6)')                                                   !
! write(6,*)'Real number ',trim(tmpr),' with words at the end.'           !
!                                                                         !
! Note : DON'T USE R2S IN A WRITE STATEMENT SINCE THIS IS ILLEGAL         !
! IN FORTRAN90 (ALTHOUGH NOT IN FORTRAN200X). IF ANYONE HAS TIME, FEEL    !
! FREE TO WRITE A VERSION OF THIS WHICH ISN'T ILLEGAL - SIMILAR TO        !
! I2S ABOVE - SO THAT PEOPLE WHO HAVEN'T READ THIS NOTE DON'T FEEL        !
! TEMPTED TO CALL R2S IN A WRITE STATEMENT.                               !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: r
 CHARACTER(*),INTENT(in) :: real_format

! if(len(r2s)>0)then
  write(r2s,real_format)r
  r2s=adjustl(r2s)
! endif

 END FUNCTION r2s

 SUBROUTINE gvsort(r,ip,alow,ahi,istart,nentry,jp,n,nlim,toll)
!-------------------------------------------------------------------!
! Sort r (indexed by ip) by increasing value.                       !
!-------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER ip(*),istart(*),nentry(*),jp(*),n,nlim
 REAL(dp) r(*),alow(*),ahi(*),toll
 INTEGER ist,j,k,l,m,ne,nen,nep,nseg
 REAL(dp) ah,al,al2,rrr,sss
 al=r(1) ; ip(1)=1 ; ah=al
 do l=2,n
  ip(l)=l
  rrr=r(l)
  al=min(al,rrr)
  ah=max(ah,rrr)
 enddo
 nseg=1 ; istart(1)=0 ; nentry(1)=n
 alow(1)=al ; ahi(1)=ah
1 al=alow(nseg)
 ah=ahi(nseg)
 if((ah-al)<toll)goto 2
 nen=nentry(nseg)
 ist=istart(nseg)
 if(nen>=9)then
! Bisect the list
  sss=(al+ah)*0.5d0
  ahi(nseg+1)=ah ; al2=ah ; ah=al
  ne=0 ; nep=0
  do l=1,nen
   k=ip(ist+l)
   rrr=r(k)
   if(rrr>=sss)then
    nep=nep+1
    jp(nep)=k
    al2=min(al2,rrr)
    cycle
   endif
   ne=ne+1
   ip(ist+ne)=k
   ah=max(ah,rrr)
  enddo
  ahi(nseg)=ah
  nentry(nseg)=ne
  ist=ist+ne
  if(ist>=nlim)goto 1
  nseg=nseg+1
  istart(nseg)=ist
  nentry(nseg)=nep
  alow(nseg)=al2
  do l=1,nep
   ip(l+ist)=jp(l)
  enddo
  goto 1
 endif
! Short vector sort
 ne=nen-1
 do l=1,ne
  k=ip(l+ist)
  rrr=r(k)
  nep=l+1
  do m=nep,nen
   j=ip(m+ist)
   sss=r(j)
   if(sss>=rrr)cycle
   ip(m+ist)=k
   k=j ; rrr=sss
  enddo
  ip(l+ist)=k
 enddo
2 nseg=nseg-1
 if(nseg/=0)goto 1
 END SUBROUTINE gvsort

END MODULE bwfdet_mod

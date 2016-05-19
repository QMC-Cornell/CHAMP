      subroutine grad_hess_jas_mpi
! Written by Claudia Filippi.  Modified by Cyrus Umrigar.

# if defined (MPI)

      use optim_mod
      use contrl_opt2_mod
      use contrl_opt_mod
      use gradhessder_mod
      use gradjerr_mod
      implicit real*8(a-h,o-z)
      integer MPARM2

      common /gradjerrb/ ngrad_jas_blocks,ngrad_jas_bcum,nb_current

      dimension collect(nparm),collect2(nparm,nparm)

      if(igradhess.eq.0) return

      MPARM2=nparm*nparm
! Note, to do: error is not collected

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      call mpi_allreduce(dj,collect,nparm &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 10 i=1,nparm
  10    dj(i)=collect(i)

! Note that if we use mpi_reduce rather than mpi_allreduce we need a '0,' before the MPI_COMM_WORLD
      call mpi_allreduce(dj_e,collect,nparm &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 20 i=1,nparm
  20    dj_e(i)=collect(i)

      call mpi_allreduce(de,collect,nparm &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 25 i=1,nparm
  25    de(i)=collect(i)

      call mpi_allreduce(de_e,collect,nparm &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 26 i=1,nparm
  26    de_e(i)=collect(i)

      call mpi_allreduce(e2,collect,nparm &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 27 i=1,nparm
  27    e2(i)=collect(i)

      call mpi_allreduce(dj_e2,collect,nparm &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 28 i=1,nparm
  28    dj_e2(i)=collect(i)

      call mpi_allreduce(dj_de,collect2,MPARM2 &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 30 i=1,nparm
        do 30 j=1,nparm
  30      dj_de(i,j)=collect2(i,j)

      call mpi_allreduce(dj_dj,collect2,MPARM2 &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 40 i=1,nparm
        do 40 j=1,nparm
  40      dj_dj(i,j)=collect2(i,j)

      call mpi_allreduce(dj_dj_e,collect2,MPARM2 &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 50 i=1,nparm
        do 50 j=1,nparm
  50      dj_dj_e(i,j)=collect2(i,j)

      call mpi_allreduce(d2j,collect2,MPARM2 &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 52 i=1,nparm
        do 52 j=1,nparm
  52      d2j(i,j)=collect2(i,j)

      call mpi_allreduce(d2j_e,collect2,MPARM2 &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 54 i=1,nparm
        do 54 j=1,nparm
  54      d2j_e(i,j)=collect2(i,j)

      call mpi_allreduce(de_de,collect2,MPARM2 &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 56 i=1,nparm
        do 56 j=1,nparm
  56      de_de(i,j)=collect2(i,j)

      if(ngrad_jas_blocks.gt.0) then
        call mpi_allreduce(ngrad_jas_bcum,ngrad_jas_collect,1 &
     &     ,mpi_integer,mpi_sum,MPI_COMM_WORLD,ierr)

        ngrad_jas_bcum=ngrad_jas_collect

        call mpi_allreduce(grad_jas_bcum,collect,nparmj &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

        do 60 i=1,nparmj
  60      grad_jas_bcum(i)=collect(i)

        call mpi_allreduce(grad_jas_bcm2,collect,nparmj &
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

        do 70 i=1,nparmj
  70      grad_jas_bcm2(i)=collect(i)

      endif

! these averages should be set to zero on the slaves but grad_hess_jas_allreduce
! is only called at the end of run (differently than prop_allreduce) and
! only the master writes to output and dumper

      call mpi_barrier(MPI_COMM_WORLD,ierr)

!     write(6,'(''nparm='',10i5)') nparm,MPARM2
!     write(6,'(''dj='',10f14.5)') (dj(i),i=1,nparm)
!     write(6,'(''de='',10f14.5)') (de(i),i=1,nparm)
!     do i=1,nparm
!       write(6,'(''dj_dj='',10f14.5)') (dj_dj(i,j),j=1,nparm)
!     enddo
!     do i=1,nparm
!       write(6,'(''dj_dj_e='',10f14.5)') (dj_dj_e(i,j),j=1,nparm)
!     enddo
!     do i=1,nparm
!       write(6,'(''dj_de='',10f14.5)') (dj_de(i,j),j=1,nparm)
!     enddo

# endif

      return
      end

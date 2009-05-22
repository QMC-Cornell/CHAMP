      subroutine grad_hess_jas_mpi
c Written by Claudia Filippi.  Modified by Cyrus Umrigar.

# if defined (MPI)

      use optim_mod
      use contrl_opt2_mod
      use contrl_opt_mod
      implicit real*8(a-h,o-z)

!JT      include '../vmc.h'
!JT      include '../../fit/fit.h'
!JT      include 'mpif.h'

      parameter(MPARM2=MPARM*MPARM)

!JT      common /contrl_opt/ nparm,nsig,ncalls,iopt,ipr_opt
      common /gradhessder/ dj(MPARM),dj_e(MPARM),dj_de(MPARM,MPARM),dj_dj(MPARM,MPARM),dj_dj_e(MPARM,MPARM)
     &,de(MPARM),d2j(MPARM,MPARM),d2j_e(MPARM,MPARM),de_e(MPARM),e2(MPARM),dj_e2(MPARM),de_de(MPARM,MPARM)

      common /gradjerr/ grad_jas_bcum(MPARMJ),grad_jas_bcm2(MPARMJ),
     &dj_e_bsum(MPARMJ),dj_bsum(MPARMJ),dj_e_save(MPARMJ),dj_save(MPARMJ),e_bsum

      common /gradjerrb/ ngrad_jas_blocks,ngrad_jas_bcum,nb_current

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase

!JT      common /contrl_opt2/ igradhess,iadd_diag_opt

      dimension collect(MPARM),collect2(MPARM,MPARM)

      if(igradhess.eq.0) return

c Note, to do: error is not collected

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      call mpi_allreduce(dj,collect,nparm
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 10 i=1,nparm
  10    dj(i)=collect(i)

c Note that if we use mpi_reduce rather than mpi_allreduce we need a '0,' before the MPI_COMM_WORLD
      call mpi_allreduce(dj_e,collect,nparm
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 20 i=1,nparm
  20    dj_e(i)=collect(i)

      call mpi_allreduce(de,collect,nparm
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 25 i=1,nparm
  25    de(i)=collect(i)

      call mpi_allreduce(de_e,collect,nparm
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 26 i=1,nparm
  26    de_e(i)=collect(i)

      call mpi_allreduce(e2,collect,nparm
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 27 i=1,nparm
  27    e2(i)=collect(i)

      call mpi_allreduce(dj_e2,collect,nparm
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 28 i=1,nparm
  28    dj_e2(i)=collect(i)

      call mpi_allreduce(dj_de,collect2,MPARM2
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 30 i=1,nparm
        do 30 j=1,nparm
  30      dj_de(i,j)=collect2(i,j)

      call mpi_allreduce(dj_dj,collect2,MPARM2
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 40 i=1,nparm
        do 40 j=1,nparm
  40      dj_dj(i,j)=collect2(i,j)

      call mpi_allreduce(dj_dj_e,collect2,MPARM2
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 50 i=1,nparm
        do 50 j=1,nparm
  50      dj_dj_e(i,j)=collect2(i,j)

      call mpi_allreduce(d2j,collect2,MPARM2
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 52 i=1,nparm
        do 52 j=1,nparm
  52      d2j(i,j)=collect2(i,j)

      call mpi_allreduce(d2j_e,collect2,MPARM2
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 54 i=1,nparm
        do 54 j=1,nparm
  54      d2j_e(i,j)=collect2(i,j)

      call mpi_allreduce(de_de,collect2,MPARM2
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 56 i=1,nparm
        do 56 j=1,nparm
  56      de_de(i,j)=collect2(i,j)

      if(ngrad_jas_blocks.gt.0) then
        call mpi_allreduce(ngrad_jas_bcum,ngrad_jas_collect,1
     &     ,mpi_integer,mpi_sum,MPI_COMM_WORLD,ierr)

        ngrad_jas_bcum=ngrad_jas_collect

        call mpi_allreduce(grad_jas_bcum,collect,nparmj
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

        do 60 i=1,nparmj
  60      grad_jas_bcum(i)=collect(i)

        call mpi_allreduce(grad_jas_bcm2,collect,nparmj
     &     ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

        do 70 i=1,nparmj
  70      grad_jas_bcm2(i)=collect(i)

      endif

c these averages should be set to zero on the slaves but grad_hess_jas_allreduce
c is only called at the end of run (differently than prop_allreduce) and
c only the master writes to output and dumper

      call mpi_barrier(MPI_COMM_WORLD,ierr)

c     write(6,'(''MPARMJ='',10i5)') MPARM,MPARM2
c     write(6,'(''dj='',10f14.5)') (dj(i),i=1,nparm)
c     write(6,'(''de='',10f14.5)') (de(i),i=1,nparm)
c     do i=1,nparm
c       write(6,'(''dj_dj='',10f14.5)') (dj_dj(i,j),j=1,nparm)
c     enddo
c     do i=1,nparm
c       write(6,'(''dj_dj_e='',10f14.5)') (dj_dj_e(i,j),j=1,nparm)
c     enddo
c     do i=1,nparm
c       write(6,'(''dj_de='',10f14.5)') (dj_de(i,j),j=1,nparm)
c     enddo

# endif

      return
      end

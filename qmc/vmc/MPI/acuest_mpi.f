      subroutine acuest_mpi
! MPI version created by Claudia Filippi starting from serial version
! routine to accumulate estimators for energy etc.

# if defined (MPI)

      use variables_mod
      use constants_mod
      use mpi_mod
      use atom_mod
      use config_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use forcepar_mod
      use pseudo_mod
      use contr3_mod
      use div_v_mod
      use distance_mod
      use qua_mod
      use forcest_mod
      use denupdn_mod
      use stepv_mod
      use forcewt_mod
      use estsum_mod
      use zigzag_mod
      implicit real*8(a-h,o-z)

      real(dp) :: esum_collect(nforce), wcollect(nforce)
      real(dp) :: pesum_collect, peisum_collect, tpbsum_collect, tjfsum_collect, r1sum_collect, r2sum_collect, accsum_collect
      real(dp) :: r3sum_collect, r4sum_collect
      real(dp) :: d_node_log_collect

      dimension zzsum_collect(nzzvars)
! statement function for error calculation
!     err(x,x2)=dsqrt(dabs(x2/iblk-(x/iblk)**2)/iblk)
!      err(x,x2,i)=dsqrt(abs(x2/wcum(i)-(x/wcum(i))**2)/iblk) ! JT: commented out because not used

! xsum = sum of values of x from metrop
! xnow = average of values of x from metrop
! xcum = accumulated sums of xnow
! xcm2 = accumulated sums of xnow**2
! xave = current average value of x
! xerr = current error of x

! eaverage(1-MFORCE) = esum(1-MFORCE)
! eaverage(MFORCE+1) = pesum
! eaverage(MFORCE+2) = peisum
! eaverage(MFORCE+3) = tpbsum
! eaverage(MFORCE+4) = tjfsum
! eaverage(MFORCE+5) = r2sum
! eaverage(MFORCE+6) = accsum


! Note we do not reduce the 1-electron move quantities here because
! they are not printed out from acuest.  So we just cumulate the
! quantities on individual processors, and reduce the cumulated
! quantities in finwrt_mpi

!JT      call mpi_allreduce(eaverage,ecollect,MFORCE+6,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

!JT: break down eaverage into its pieces because the variables are now in a module
      call mpi_allreduce(esum,esum_collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(pesum,pesum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(peisum,peisum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tpbsum,tpbsum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tjfsum,tjfsum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(r1sum,r1sum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(r2sum,r2sum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(r3sum,r3sum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(r4sum,r4sum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      if(izigzag.gt.0) then
        call mpi_allreduce(zzsum,zzsum_collect,nzzvars,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      endif
      call mpi_allreduce(accsum,accsum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(d_node_log_sum,d_node_log_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wsum,wcollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

! Warning this flush and barrier should not be necessary
      call systemflush(6)
      call mpi_barrier(MPI_COMM_WORLD,ierr)

!JT      do 10 ifr=1,MFORCE+6
!JT        eaverage(ifr)=ecollect(ifr)
!JT   10   if(ifr.le.nforce) wsum(ifr)=wcollect(ifr)

      esum(1:nforce) = esum_collect (1:nforce)
      wsum(1:nforce) = wcollect (1:nforce)
      pesum = pesum_collect
      peisum = peisum_collect
      tpbsum = tpbsum_collect
      tjfsum = tjfsum_collect
      r1sum = r1sum_collect
      r2sum = r2sum_collect
      r3sum = r3sum_collect
      r4sum = r4sum_collect
      if(izigzag.gt.0) then
      zzsum(:) = zzsum_collect(:)
      endif
      accsum = accsum_collect
      d_node_log_sum = d_node_log_collect

!     pesum=ecollect(MFORCE+1)
!     peisum=ecollect(MFORCE+2)
!     tpbsum=ecollect(MFORCE+3)
!     tjfsum=ecollect(MFORCE+4)
!     r2sum=ecollect(MFORCE+5)
!     accsum=ecollect(MFORCE+6)

# endif

      return
      end

      subroutine acuest_mpi
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

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
      real(dp) :: pesum_collect, peisum_collect, tpbsum_collect, tjfsum_collect, r2sum_collect, accsum_collect
      real(dp) :: d_node_log_collect, walker_weights_sum_block_collect

      dimension zzsum_collect(nzzvars)
c statement function for error calculation
c     err(x,x2)=dsqrt(dabs(x2/iblk-(x/iblk)**2)/iblk)
!      err(x,x2,i)=dsqrt(abs(x2/wcum(i)-(x/wcum(i))**2)/iblk) ! JT: commented out because not used

c xsum = sum of values of x from metrop
c xnow = average of values of x from metrop
c xcum = accumulated sums of xnow
c xcm2 = accumulated sums of xnow**2
c xave = current average value of x
c xerr = current error of x

c eaverage(1-MFORCE) = esum(1-MFORCE)
c eaverage(MFORCE+1) = pesum
c eaverage(MFORCE+2) = peisum
c eaverage(MFORCE+3) = tpbsum
c eaverage(MFORCE+4) = tjfsum
c eaverage(MFORCE+5) = r2sum
c eaverage(MFORCE+6) = accsum


c Note we do not reduce the 1-electron move quantities here because
c they are not printed out from acuest.  So we just cumulate the
c quantities on individual processors, and reduce the cumulated
c quantities in finwrt_mpi

!JT      call mpi_allreduce(eaverage,ecollect,MFORCE+6,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

!JT: break down eaverage into its pieces because the variables are now in a module
      call mpi_allreduce(esum,esum_collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(pesum,pesum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(peisum,peisum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tpbsum,tpbsum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tjfsum,tjfsum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(r2sum,r2sum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(zzsum,zzsum_collect,nzzvars,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(accsum,accsum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(d_node_log_sum,d_node_log_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(walker_weights_sum_block,walker_weights_sum_block_collect,1,mpi_double_precision,
     & mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wsum,wcollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

c Warning this flush and barrier should not be necessary
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
      r2sum = r2sum_collect
      zzsum(:) = zzsum_collect(:)
      accsum = accsum_collect
      d_node_log_sum = d_node_log_collect
      walker_weights_sum_block = walker_weights_sum_block_collect

c     pesum=ecollect(MFORCE+1)
c     peisum=ecollect(MFORCE+2)
c     tpbsum=ecollect(MFORCE+3)
c     tjfsum=ecollect(MFORCE+4)
c     r2sum=ecollect(MFORCE+5)
c     accsum=ecollect(MFORCE+6)

# endif

      return
      end

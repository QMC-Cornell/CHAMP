      subroutine walksav_det_mpi
! Written by Claudia Filippi, modified by Cyrus Umrigar
! Only the entries are called, not the subroutine

# if  defined (MPI)
      use all_tools_mod
      use dets_mod
      use slater_mod
      use slaterw_mod
      use const_mod
      use branch_mod
      use gammaw_mod
      use deriv_fast_mod
      implicit real*8(a-h,o-z)

!      dimension istatus(MPI_STATUS_SIZE)

      entry send_det(irecv)

      itag=1
      call mpi_send_double_arr(orbw(1,1,nwalk),size(orbw(:,:,nwalk)),irecv,itag)
      call mpi_send_double_arr(dorbw(1,1,1,nwalk),size(dorbw(:,:,:,nwalk)),irecv,itag)
 
      call mpi_send_double_arr(aiupw(1,1,nwalk),size(aiupw(:,:,nwalk)),irecv,itag)
      call mpi_send_double_arr(deta_upw(nwalk),1,irecv,itag)
      call mpi_send_double_arr(tupw(1,1,nwalk),size(tupw(:,:,nwalk)),irecv,itag)
      call mpi_send_double_arr(detupw(1,nwalk),size(detupw(:,nwalk)),irecv,itag)
      call mpi_send_double_arr(invupw(1,1,nwalk),size(invupw(:,:,nwalk)),irecv,itag)

      call mpi_send_double_arr(aidnw(1,1,nwalk),size(aidnw(:,:,nwalk)),irecv,itag)
      call mpi_send_double_arr(deta_dnw(nwalk),1,irecv,itag)
      call mpi_send_double_arr(tdnw(1,1,nwalk),size(tdnw(:,:,nwalk)),irecv,itag)
      call mpi_send_double_arr(detdnw(1,nwalk),size(detdnw(:,nwalk)),irecv,itag)
      call mpi_send_double_arr(invdnw(1,1,nwalk),size(invdnw(:,:,nwalk)),irecv,itag)

      call mpi_send_double_arr(chiw(nwalk),1,irecv,itag)
      call mpi_send_double_arr(yupw(1,1,nwalk),size(yupw(:,:,nwalk)),irecv,itag)
      call mpi_send_double_arr(ydnw(1,1,nwalk),size(ydnw(:,:,nwalk)),irecv,itag)
      return

      entry recv_det(isend)

      itag=1
      call mpi_recv_double_arr(orbw(1,1,nwalk),size(orbw(:,:,nwalk)),isend,itag)
      call mpi_recv_double_arr(dorbw(1,1,1,nwalk),size(dorbw(:,:,:,nwalk)),isend,itag)
 
      call mpi_recv_double_arr(aiupw(1,1,nwalk),size(aiupw(:,:,nwalk)),isend,itag)
      call mpi_recv_double_arr(deta_upw(nwalk),1,isend,itag)
      call mpi_recv_double_arr(tupw(1,1,nwalk),size(tupw(:,:,nwalk)),isend,itag)
      call mpi_recv_double_arr(detupw(1,nwalk),size(detupw(:,nwalk)),isend,itag)
      call mpi_recv_double_arr(invupw(1,1,nwalk),size(invupw(:,:,nwalk)),isend,itag)

      call mpi_recv_double_arr(aidnw(1,1,nwalk),size(aidnw(:,:,nwalk)),isend,itag)
      call mpi_recv_double_arr(deta_dnw(nwalk),1,isend,itag)
      call mpi_recv_double_arr(tdnw(1,1,nwalk),size(tdnw(:,:,nwalk)),isend,itag)
      call mpi_recv_double_arr(detdnw(1,nwalk),size(detdnw(:,nwalk)),isend,itag)
      call mpi_recv_double_arr(invdnw(1,1,nwalk),size(invdnw(:,:,nwalk)),isend,itag)

      call mpi_recv_double_arr(chiw(nwalk),1,isend,itag)
      call mpi_recv_double_arr(yupw(1,1,nwalk),size(yupw(:,:,nwalk)),isend,itag)
      call mpi_recv_double_arr(ydnw(1,1,nwalk),size(ydnw(:,:,nwalk)),isend,itag)
      return

    contains

        subroutine mpi_send_double_arr(dat, dat_size, irecv, itag)
          use types_mod, only: dp
          use all_tools_mod
          implicit none
        
          real(dp), intent(in) :: dat
          integer, intent(in) :: dat_size, irecv
          integer, intent(inout) :: itag
          integer :: irequest, ierr, istatus(MPI_STATUS_SIZE)
        
          call mpi_isend(dat, dat_size, mpi_double_precision, irecv, itag, MPI_COMM_WORLD, irequest, ierr)
          call MPI_Wait(irequest, istatus, ierr)
          itag=itag+1
        end subroutine mpi_send_double_arr

        subroutine mpi_recv_double_arr(dat, dat_size, isend, itag)
          use types_mod, only: dp
          use all_tools_mod
          implicit none
          real(dp), intent(in) :: dat
          integer, intent(in) :: dat_size, isend
          integer, intent(inout) :: itag
          integer :: istatus(MPI_STATUS_SIZE), ierr
        
          call mpi_recv(dat, dat_size, mpi_double_precision, isend, itag, MPI_COMM_WORLD, istatus, ierr)
          itag=itag+1
        end subroutine mpi_recv_double_arr

# endif
      end subroutine walksav_det_mpi

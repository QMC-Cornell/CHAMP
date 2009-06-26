      subroutine walksav_jas_mpi3
c Written by Claudia Filippi
c Only the entries are called, not the subroutine

# if defined (MPI)
      use all_tools_mod
      use const_mod
      use forcepar_mod
      use force_dmc_mod
      use forcest_dmc_mod
      use jaso_mod
      use jasow_mod
      use branch_mod
      implicit real*8(a-h,o-z)

      common /mpitype/ jas_type1,jas_type2

      dimension istatus(MPI_STATUS_SIZE)

      entry send_jas_mpi3(irecv)

c     write(6,'(''send_jas_mpi3'')')
      itag=0
      call mpi_isend(fsumow(nwalk),1,mpi_double_precision,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(fjow(1,1,nwalk),3*nelec,mpi_double_precision,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(fsow(1,1,nwalk),1,jas_type1,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(fijow(1,1,1,nwalk),1,jas_type2,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Wait(irequest,istatus,ierr)

      return

      entry recv_jas_mpi3(isend)

c     write(6,'(''recv_jas_mpi3'')')
      itag=0
      call mpi_recv(fsumow(nwalk),1,mpi_double_precision,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fjow(1,1,nwalk),3*nelec,mpi_double_precision,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fsow(1,1,nwalk),1,jas_type1,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fijow(1,1,1,nwalk),1,jas_type2,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)

# endif
      return
      end

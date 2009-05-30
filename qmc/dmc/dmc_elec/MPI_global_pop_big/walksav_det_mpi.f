      subroutine walksav_det_mpi
c     subroutine walksav_det_mpi(iw)
c Written by Claudia Filippi, modified by Cyrus Umrigar

# if  defined (MPI)
      use all_tools_mod
      use dets_mod
      use slater_mod
      use slaterw_mod
      use const_mod
      use branch_mod
      implicit real*8(a-h,o-z)

      dimension istatus(MPI_STATUS_SIZE)

c     call walksav_det

c     return

      entry send_det(irecv)

      itag=0
      call mpi_isend(detuw(1,nwalk),ndet,mpi_double_precision,irecv
     &,itag+1,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(detdw(1,nwalk),ndet,mpi_double_precision,irecv
     &,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
      do 150 k=1,ndet
        call mpi_isend(slmuiw(1,k,nwalk),nup*nup,mpi_double_precision
     &  ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(fpuw(1,1,k,nwalk),3*nup*nup,mpi_double_precision
     &  ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(slmdiw(1,k,nwalk),ndn*ndn,mpi_double_precision
     &  ,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(fpdw(1,1,k,nwalk),3*ndn*ndn,mpi_double_precision
     &  ,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(ddeti_detiw(1,1,k,nwalk),3*nelec,mpi_double_precision
     &  ,irecv,itag+5,MPI_COMM_WORLD,irequest,ierr)
  150   itag=itag+5
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
      return

      entry recv_det(isend)

      itag=0
      call mpi_recv(detuw(1,nwalk),ndet,mpi_double_precision,isend
     &,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(detdw(1,nwalk),ndet,mpi_double_precision,isend
     &,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2
      do 160 k=1,ndet
        call mpi_recv(slmuiw(1,k,nwalk),nup*nup,mpi_double_precision
     &  ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(fpuw(1,1,k,nwalk),3*nup*nup,mpi_double_precision
     &  ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(slmdiw(1,k,nwalk),ndn*ndn,mpi_double_precision
     &  ,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(fpdw(1,1,k,nwalk),3*ndn*ndn,mpi_double_precision
     &  ,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(ddeti_detiw(1,1,k,nwalk),3*nelec,mpi_double_precision
     &  ,isend,itag+5,MPI_COMM_WORLD,istatus,ierr)
  160   itag=itag+5

# endif
      return
      end

      subroutine walksav_jas_mpi2
c     subroutine walksav_jas_mpi(iw)
c Written by Claudia Filippi

# if defined (MPI)

      use all_tools_mod
      use const_mod
      use jaso_mod
      use branch_mod
      implicit real*8(a-h,o-z)

!JT!JT     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),

      common /jasow/ fsow(MELEC,MELEC,MWALK),fijow(3,MELEC,MELEC,MWALK),fsumow(MWALK),fjow(3,MELEC,MWALK)

      dimension istatus(MPI_STATUS_SIZE)

c     save fsow,fijow,fsumow,fjow

      entry send_jas_mpi2(irecv)

      itag=0
      call mpi_isend(fsumow(nwalk),1,mpi_double_precision,irecv
     &,itag+1,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(fjow(1,1,nwalk),3*nelec,mpi_double_precision,irecv
     &,itag+2,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(fsow(1,1,nwalk),MELEC*nelec,mpi_double_precision
     &,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(fijow(1,1,1,nwalk),3*MELEC*nelec
     &,mpi_double_precision,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)

      return

      entry recv_jas_mpi2(isend)

      itag=0
      call mpi_recv(fsumow(nwalk),1,mpi_double_precision,isend
     &,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fjow(1,1,nwalk),3*nelec,mpi_double_precision,isend
     &,itag+2,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fsow(1,1,nwalk),MELEC*nelec,mpi_double_precision
     &,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fijow(1,1,1,nwalk),3*MELEC*nelec
     &,mpi_double_precision,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)

# endif
      return
      end

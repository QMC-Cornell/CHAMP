c Named func_qmc_mpi to avoid name conflict with func_mpi in quench_anneal
      subroutine func_qmc_mpi(diff)
# if defined (MPI)
c MPI communication introduced by Claudia Filippi
c MPI communication redone by Richard Hennig to make it work with
c ifort compiler on OSC Pentium4 cluster.  Otherwise it mysteriously
c hangs during the mpi-gatherv on psij.
c Richard Hennig also fixed the load balancing.

      use all_tools_mod
      implicit real*8(a-h,o-z)

      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata

      common /mpioffset/ ircounts(0:MPROC),idispls(0:MPROC)

      dimension diff(*)

      ni=idispls(idtask)+1
      nscounts=ircounts(idtask)

c     RGH
c     call mpi_barrier(MPI_COMM_WORLD, ierr)

      call mpi_gatherv(diff(ni),nscounts,mpi_double_precision
     &     ,diff,ircounts,idispls,mpi_double_precision,0
     &     ,MPI_COMM_WORLD,ierr)
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(diff,ndata,mpi_double_precision,0
     &     ,MPI_COMM_WORLD,ierr)

c     RGH
c     call mpi_barrier(MPI_COMM_WORLD, ierr)

      call mpi_gatherv(psid(ni),nscounts,mpi_double_precision
     &     ,psid,ircounts,idispls,mpi_double_precision,0
     &     ,MPI_COMM_WORLD,ierr)
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(psid,ndata,mpi_double_precision,0
     &     ,MPI_COMM_WORLD,ierr)

c     RGH
c     call mpi_barrier(MPI_COMM_WORLD, ierr)

      call mpi_gatherv(psij(ni),nscounts,mpi_double_precision
     &     ,psij,ircounts,idispls,mpi_double_precision,0
     &     ,MPI_COMM_WORLD,ierr)
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(psij,ndata,mpi_double_precision,0
     &     ,MPI_COMM_WORLD,ierr)

c     RGH
c     call mpi_barrier(MPI_COMM_WORLD, ierr)

# endif
      return
      end

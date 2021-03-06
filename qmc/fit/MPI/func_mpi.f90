! Named func_qmc_mpi to avoid name conflict with func_mpi in quench_anneal
      subroutine func_qmc_mpi(diff)
# if defined (MPI)
! MPI communication introduced by Claudia Filippi

      use all_tools_mod
      use confg_mod
      use mpi_mod
      use mpioffset_mod
      implicit real*8(a-h,o-z)

      dimension diff(*)

      ni=idispls(idtask)+1
      nscounts=ircounts(idtask)

      if(idtask.eq.0) then
        call mpi_gatherv(MPI_IN_PLACE,nscounts,mpi_double_precision, &
            diff,ircounts,idispls,mpi_double_precision, 0,MPI_COMM_WORLD,ierr)
      else
        call mpi_gatherv(diff(ni),nscounts,mpi_double_precision, &
            diff,ircounts,idispls,mpi_double_precision, 0,MPI_COMM_WORLD,ierr)
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(diff,ndata,mpi_double_precision,0, MPI_COMM_WORLD,ierr)

      if(idtask.eq.0) then
        call mpi_gatherv(MPI_IN_PLACE,nscounts,mpi_double_precision, & 
            psid,ircounts,idispls,mpi_double_precision, 0,MPI_COMM_WORLD,ierr)
      else
        call mpi_gatherv(psid(ni),nscounts,mpi_double_precision, & 
            psid,ircounts,idispls,mpi_double_precision, 0,MPI_COMM_WORLD,ierr)
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(psid,ndata,mpi_double_precision,0, MPI_COMM_WORLD,ierr)

      if(idtask.eq.0) then
        call mpi_gatherv(MPI_IN_PLACE,nscounts,mpi_double_precision, &
            psij,ircounts,idispls,mpi_double_precision,0, MPI_COMM_WORLD,ierr)
      else
        call mpi_gatherv(psij(ni),nscounts,mpi_double_precision, &
            psij,ircounts,idispls,mpi_double_precision,0, MPI_COMM_WORLD,ierr)
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(psij,ndata,mpi_double_precision,0, MPI_COMM_WORLD,ierr)

# endif
      return
      end

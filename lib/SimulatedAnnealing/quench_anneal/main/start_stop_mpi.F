#ifdef MPI
      subroutine start_mpi(initialize,jproc,mproc,show_load_chunks1)
c start-up routine
c initialize	= .true. to start MPI and IO up (in)
c jproc	= number of current processor (out)
c mproc	= total number of processors  (out)
c show_load_chunks1	= .true. to show load ballance output (in)
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      include '../include/IO_stop.h'
      logical initialize,show_load_chunks1

      if(initialize) then
        call MPI_INIT(ierr_init)
        call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
        call open_stio(iproc)
      else
        call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
      endif
      call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
      jproc=iproc
      mproc=nproc
      show_load_chunks=show_load_chunks1
      return
      end
c-----------------------------------------------------------------------

      subroutine stop_mpi(finalize)

      implicit real*8(a-h,o-z)
      include '../include/IO_stop.h'
      include 'mpif.h'
      logical finalize

      if(show_load_chunks) write(6,'((a),(10g9.2))')'stop_mpi: times ',
     &  (time_g_tot(i),i=1,nproc)
      if(finalize) call mpi_finalize(ierr_final)
      more_data_qa=0
      iatom=-1
      first_call=.true.
      return
      end
#else
c For those compilers that don't like codeless files
      subroutine never_call()
      write(6,*)'never_call: this routine should never be called'
      stop 'never_call was called'
      end
#endif


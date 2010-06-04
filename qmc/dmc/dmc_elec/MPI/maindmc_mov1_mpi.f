      subroutine maindmc_mov1_mpi
c MPI version created by Claudia Filippi and Cyrus Umrigar

# if defined(MPI)
      use all_tools_mod
      use dmc_mod
      use contr3_mod
      use const_mod
      implicit real*8(a-h,o-z)

      common /mpitype/ jas_type1,jas_type2

      dimension iblocklen(nelec),idispl(nelec)

      call MPI_ATTR_GET(MPI_COMM_WORLD,MPI_TAG_UB,ivalue,flag,ierr)
c     write(6,*) 'In main.f from MPI_ATTR_GET',ivalue
      if(ierr.ne.0) write(6,*) 'Warning:? in main.f from MPI_ATTR_GET',ierr

c dmc_init and opt_wf are the same for serial and parallel, open_files is different.
      call dmc_init

c To minimize the length of what is passed, use jas_typ to pick certain
c elements of matrix
      if(mode.eq.'dmc_mov1_mpi3') then
        do 1 i=1,nelec
          iblocklen(i)=nelec-(i-1)
   1      idispl(i)=nelec*(i-1)+(i-1)
        call mpi_type_indexed(nelec,iblocklen,idispl,mpi_double_precision,jas_type1,ierr)
        call mpi_type_commit(jas_type1,ierr)

        do 2 i=1,nelec
          iblocklen(i)=ndim*nelec
   2      idispl(i)=ndim*nelec*(i-1)
        call mpi_type_indexed(nelec,iblocklen,idispl,mpi_double_precision,jas_type2,ierr)
        call mpi_type_commit(jas_type2,ierr)
      endif

!      call open_files_mpi
      call opt_wf_dmc

# endif
      end

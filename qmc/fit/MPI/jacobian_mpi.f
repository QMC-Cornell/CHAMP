      subroutine jacobian_mpi(ndata2,nanalytic,nparm,ajac)
c Written by Claudia Filippi and Cyrus Umrigar
# if defined (MPI)

      use all_tools_mod
      implicit real*8(a-h,o-z)

      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata

      common /mpioffset/ ircounts(0:nprocx),idispls(0:nprocx)

      dimension ajac(ndata2,nparm)

      ni=idispls(idtask)+1
      nscounts=ircounts(idtask)

      nnumerical=nparm-nanalytic
      do 130 iparm=nnumerical+1,nparm
        call mpi_gatherv(ajac(ni,iparm),nscounts,mpi_double_precision
     &  ,ajac(1,iparm),ircounts,idispls,mpi_double_precision,0
     &  ,MPI_COMM_WORLD,ierr)

  130   call mpi_bcast(ajac(1,iparm),ndata,mpi_double_precision,0
     &  ,MPI_COMM_WORLD,ierr)
# endif
      return
      end

      subroutine send_walker(irecv)
! Written by Claudia Filippi

# if defined (MPI)
      use all_tools_mod
      use const_mod
      use forcepar_mod
      use force_dmc_mod
      use jacobsave_mod
      use forcest_dmc_mod
      use config_dmc_mod
      use branch_mod
      use stats_mod
      use age_mod
      use qua_mod, only: quadr, quadx, nquad !TA
      use contrldmc_mod, only: tmoves !TA
      use atom_mod, only: ncent !TA
      use dim_mod, only: ndim !TA
      use velratio_mod
      implicit real*8(a-h,o-z)

      dimension istatus(MPI_STATUS_SIZE)

      call mpi_isend(wt(nwalk),1,mpi_double_precision,irecv,1 &
     &,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(iage(nwalk),1,mpi_integer,irecv,2 &
     &,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)

! MPI_Request_free is faster than MPI_wait and is OK if there is a barrier
! before variables are reused.
      itag=2
      do 15 ifr=1,nforce
        call mpi_isend(ajacold(nwalk,ifr),1,mpi_double_precision,irecv &
     &  ,itag+1,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        call mpi_isend(eoldw(nwalk,ifr),1,mpi_double_precision,irecv &
     &  ,itag+2,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        call mpi_isend(psidow(nwalk,ifr),1,mpi_double_precision,irecv &
     &  ,itag+3,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        call mpi_isend(psijow(nwalk,ifr),1,mpi_double_precision,irecv &
     &  ,itag+4,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        call mpi_isend(peow(nwalk,ifr),1,mpi_double_precision,irecv &
     &  ,itag+5,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(peiow(nwalk,ifr),1,mpi_double_precision,irecv &
     &  ,itag+5,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        call mpi_isend(d2ow(nwalk,ifr),1,mpi_double_precision,irecv &
     &  ,itag+6,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        call mpi_isend(pwt(nwalk,ifr),1,mpi_double_precision,irecv &
     &  ,itag+7,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        call mpi_isend(fratio(nwalk,ifr),1,mpi_double_precision,irecv &
     &  ,itag+8,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        call mpi_isend(voldw(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
     &  ,irecv,itag+9,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        call mpi_isend(xoldw(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
     &  ,irecv,itag+10,MPI_COMM_WORLD,irequest,ierr)
        itag=itag+10
        call MPI_Request_free(irequest,ierr)
        do 15 ip=0,nwprod-1
        itag=itag+1
        call mpi_isend(wthist(nwalk,ip,ifr),1,mpi_double_precision,irecv &
     &  ,itag,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
   15 continue

      if (tmoves) then
        itag=itag+1
        call mpi_isend(quadr(1,1,nwalk),nquad*ncent*nelec,mpi_double_precision,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
        itag=itag+1
        call mpi_isend(quadx(1,1,1,nwalk),ndim*nquad*ncent*nelec,mpi_double_precision,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Request_free(irequest,ierr)
      endif

!     call send_det(itag,irecv)
!     call send_jas(itag,irecv)

!     nwalk=nwalk-1

      return

      entry recv_walker(isend)

!     nwalk=nwalk+1

      call mpi_recv(wt(nwalk),1,mpi_double_precision,isend,1 &
     &,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(iage(nwalk),1,mpi_integer,isend,2 &
     &,MPI_COMM_WORLD,istatus,ierr)

      itag=2
      do 25 ifr=1,nforce
        call mpi_recv(ajacold(nwalk,ifr),1,mpi_double_precision,isend &
     &  ,itag+1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(eoldw(nwalk,ifr),1,mpi_double_precision,isend &
     &  ,itag+2,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(psidow(nwalk,ifr),1,mpi_double_precision,isend &
     &  ,itag+3,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(psijow(nwalk,ifr),1,mpi_double_precision,isend &
     &  ,itag+4,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(peow(nwalk,ifr),1,mpi_double_precision,isend &
     &  ,itag+5,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(peiow(nwalk,ifr),1,mpi_double_precision,isend &
     &  ,itag+5,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(d2ow(nwalk,ifr),1,mpi_double_precision,isend &
     &  ,itag+6,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(pwt(nwalk,ifr),1,mpi_double_precision,isend &
     &  ,itag+7,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(fratio(nwalk,ifr),1,mpi_double_precision,isend &
     &  ,itag+8,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(voldw(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
     &  ,isend,itag+9,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(xoldw(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
     &  ,isend,itag+10,MPI_COMM_WORLD,istatus,ierr)
        itag=itag+10
        do 25 ip=0,nwprod-1
        itag=itag+1
        call mpi_recv(wthist(nwalk,ip,ifr),1,mpi_double_precision,isend &
     &  ,itag,MPI_COMM_WORLD,istatus,ierr)
   25 continue

      if (tmoves) then
        itag=itag+1
        call mpi_recv(quadr(1,1,nwalk),nquad*ncent*nelec,mpi_double_precision,isend,itag,MPI_COMM_WORLD,istatus,ierr)
        itag=itag+1
        call mpi_recv(quadx(1,1,1,nwalk),ndim*nquad*ncent*nelec,mpi_double_precision,isend,itag,MPI_COMM_WORLD,istatus,ierr)
      endif

!     call recv_det(itag,isend)
!     call recv_jas(itag,isend)

# endif
      return
      end

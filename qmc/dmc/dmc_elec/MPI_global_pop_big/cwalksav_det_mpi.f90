      subroutine cwalksav_det_mpi(iw)
! Same as walksav_det, adapted to complex orbitals
! by A.D.Guclu, Oct2005.

# if defined (MPI)

      use all_tools_mod
      use cslater_mod
      use cslaterw_mod
      use dorb_mod
      use dets_mod
      use const_mod
      use dim_mod
      use forcepar_mod
      use force_dmc_mod
      use forcest_dmc_mod
      use branch_mod
      implicit real*8(a-h,o-z)

      dimension istatus(MPI_STATUS_SIZE)

! allocate memory (will allocate only if it is not already allocated):
      n2=nelec*nelec
      call alloc('cslmui',cslmui,n2,ndetup)
      call alloc('cslmdi',cslmdi,n2,ndetdn)
      call alloc('cfpu',cfpu,ndim,n2,ndetup)
      call alloc('cfpd',cfpd,ndim,n2,ndetdn)
      call alloc('cdetu',cdetu,ndetup)
      call alloc('cdetd',cdetd,ndetdn)
      call alloc('cddeti_deti',cddeti_deti,ndim,nelec,ndet)

! nwalk does not represent the actual size of the arrays.
! it is not clear to me what should be the actual size so
! so for the moment I use MWALK instead of nwalk :
      ndimw=MWALK
      call alloc('cslmuiw',cslmuiw,n2,ndetup,ndimw)
      call alloc('cslmdiw',cslmdiw,n2,ndetdn,ndimw)
      call alloc('cfpuw',cfpuw,ndim,n2,ndetup,ndimw)
      call alloc('cfpdw',cfpdw,ndim,n2,ndetdn,ndimw)
      call alloc('cdetuw',cdetuw,ndetup,ndimw)
      call alloc('cdetdw',cdetdw,ndetdn,ndimw)
      call alloc('cddeti_detiw',cddeti_detiw,ndim,nelec,ndet,ndimw)

!     call cwalksav_det

!     return

      entry csend_det(irecv)

      itag=0
      call mpi_isend(cdetuw(1,nwalk),ndet,mpi_double_complex,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(cdetdw(1,nwalk),ndet,mpi_double_complex,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2
      call MPI_Wait(irequest,istatus,ierr)
      do 150 k=1,ndet
        call mpi_isend(cslmuiw(1,k,nwalk),nup*nup,mpi_double_complex,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(cfpuw(1,1,k,nwalk),3*nup*nup,mpi_double_complex,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(cslmdiw(1,k,nwalk),ndn*ndn,mpi_double_complex,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(cfpdw(1,1,k,nwalk),3*ndn*ndn,mpi_double_complex,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(cddeti_detiw(1,1,k,nwalk),3*nelec,mpi_double_complex,irecv,itag+5,MPI_COMM_WORLD,irequest,ierr)
  150   itag=itag+5
        call MPI_Wait(irequest,istatus,ierr)
      return

      entry crecv_det(isend)

      itag=0
      call mpi_recv(cdetuw(1,nwalk),ndet,mpi_double_complex,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(cdetdw(1,nwalk),ndet,mpi_double_complex,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2
      do 160 k=1,ndet
        call mpi_recv(cslmuiw(1,k,nwalk),nup*nup,mpi_double_complex,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(cfpuw(1,1,k,nwalk),3*nup*nup,mpi_double_complex,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(cslmdiw(1,k,nwalk),ndn*ndn,mpi_double_complex,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(cfpdw(1,1,k,nwalk),3*ndn*ndn,mpi_double_complex,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(cddeti_detiw(1,1,k,nwalk),3*nelec,mpi_double_complex,isend,itag+5,MPI_COMM_WORLD,istatus,ierr)
  160   itag=itag+5

# endif
      return
      end

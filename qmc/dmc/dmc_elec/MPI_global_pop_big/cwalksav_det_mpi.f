      subroutine cwalksav_det_mpi(iw)
c Same as walksav_det, adapted to complex orbitals
c by A.D.Guclu, Oct2005.

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

c complex locals:
c      complex*16 cslmuiw,cslmdiw,cfpuw,cfpdw,cdetuw,cdetdw,cddeti_detiw

c complex commons:
c      complex*16 cslmui,cslmdi,cfpu,cfpd,cfppu,cfppd,cdetu,cdetd,cddeti_deti,cd2edeti_deti
c      complex*16 cdeti_det,cddeti_det,cd2deti_det,cd2det_det


c      common /cslater/ cslmui(MMAT_DIM,MDET),cslmdi(MMAT_DIM,MDET)
c     &,cfpu(3,MMAT_DIM,MDET),cfpd(3,MMAT_DIM,MDET)
c     &,cfppu(MMAT_DIM,MDET),cfppd(MMAT_DIM,MDET)
c     &,cdetu(MDET),cdetd(MDET)
c     &,cddeti_deti(3,MELEC,MDET),cd2edeti_deti(MELEC,MDET),cdeti_det(MCSF),cddeti_det(3,MELEC,MCSF),cd2deti_det(MCSF),cd2det_det

c      dimension cslmuiw(MMAT_DIM,MDET,MWALK)
c     &,cslmdiw(MMAT_DIM,MDET,MWALK)
c     &,cfpuw(3,MMAT_DIM,MDET,MWALK),cfpdw(3,MMAT_DIM,MDET,MWALK)
c     &,cdetuw(MDET,MWALK),cdetdw(MDET,MWALK)
c     &,cddeti_detiw(3,MELEC,MDET,MWALK)

      dimension istatus(MPI_STATUS_SIZE)

c allocate memory (will allocate only if it is not already allocated):
      n2=nelec*nelec
      call alloc('cslmui',cslmui,n2,ndetup)
      call alloc('cslmdi',cslmdi,n2,ndetdn)
      call alloc('cfpu',cfpu,ndim,n2,ndetup)
      call alloc('cfpd',cfpd,ndim,n2,ndetdn)
      call alloc('cdetu',cdetu,ndetup)
      call alloc('cdetd',cdetd,ndetdn)
      call alloc('cddeti_deti',cddeti_deti,ndim,nelec,ndet)

c nwalk does not represent the actual size of the arrays.
c it is not clear to me what should be the actual size so
c so for the moment I use MWALK instead of nwalk :
      ndimw=MWALK
      call alloc('cslmuiw',cslmuiw,n2,ndetup,ndimw)
      call alloc('cslmdiw',cslmdiw,n2,ndetdn,ndimw)
      call alloc('cfpuw',cfpuw,ndim,n2,ndetup,ndimw)
      call alloc('cfpdw',cfpdw,ndim,n2,ndetdn,ndimw)
      call alloc('cdetuw',cdetuw,ndetup,ndimw)
      call alloc('cdetdw',cdetdw,ndetdn,ndimw)
      call alloc('cddeti_detiw',cddeti_detiw,ndim,nelec,ndet,ndimw)

c     call cwalksav_det

c     return

      entry csend_det(irecv)

      itag=0
      call mpi_isend(cdetuw(1,nwalk),ndet,mpi_double_complex,irecv
     &,itag+1,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(cdetdw(1,nwalk),ndet,mpi_double_complex,irecv
     &,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
      do 150 k=1,ndet
        call mpi_isend(cslmuiw(1,k,nwalk),nup*nup,mpi_double_complex
     &  ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(cfpuw(1,1,k,nwalk),3*nup*nup,mpi_double_complex
     &  ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(cslmdiw(1,k,nwalk),ndn*ndn,mpi_double_complex
     &  ,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(cfpdw(1,1,k,nwalk),3*ndn*ndn,mpi_double_complex
     &  ,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(cddeti_detiw(1,1,k,nwalk),3*nelec,mpi_double_complex
     &  ,irecv,itag+5,MPI_COMM_WORLD,irequest,ierr)
  150   itag=itag+5
C     RGH
        call MPI_Wait(irequest,istatus,ierr)
      return

      entry crecv_det(isend)

      itag=0
      call mpi_recv(cdetuw(1,nwalk),ndet,mpi_double_complex,isend
     &,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(cdetdw(1,nwalk),ndet,mpi_double_complex,isend
     &,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2
      do 160 k=1,ndet
        call mpi_recv(cslmuiw(1,k,nwalk),nup*nup,mpi_double_complex
     &  ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(cfpuw(1,1,k,nwalk),3*nup*nup,mpi_double_complex
     &  ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(cslmdiw(1,k,nwalk),ndn*ndn,mpi_double_complex
     &  ,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(cfpdw(1,1,k,nwalk),3*ndn*ndn,mpi_double_complex
     &  ,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(cddeti_detiw(1,1,k,nwalk),3*nelec,mpi_double_complex
     &  ,isend,itag+5,MPI_COMM_WORLD,istatus,ierr)
  160   itag=itag+5

# endif
      return
      end

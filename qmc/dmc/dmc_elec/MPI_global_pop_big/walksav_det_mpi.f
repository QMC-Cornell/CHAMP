      subroutine walksav_det_mpi
c     subroutine walksav_det_mpi(iw)
c Written by Claudia Filippi, modified by Cyrus Umrigar

# if  defined (MPI)
      use all_tools_mod
      use dets_mod
      use slater_mod
      use const_mod
      implicit real*8(a-h,o-z)

c     common /forcepar/ deltot(MFORCE),nforce,istrech
c     common /forcest_dmc/ fgcum(MFORCE),fgcm2(MFORCE)
c     common /force_dmc/ itausec,nwprod

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
!JT      common /slater/ slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD)
!JT     &,fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
!JT     &,fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
!JT     &,detu(MDETUD),detd(MDETUD)
!JT     &,ddeti_deti(3,MELEC,MDETUD),d2edeti_deti(MELEC,MDETUD),deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD),d2deti_det(MPARMD),d2det_det
!JT     &,detij_det(MPARMD,MPARMD)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eoldw(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      common /slaterw/ slmuiw(MMAT_DIM,MDETUD,MWALK),slmdiw(MMAT_DIM,MDETUD,MWALK)
     &,fpuw(3,MMAT_DIM,MDETUD,MWALK),fpdw(3,MMAT_DIM,MDETUD,MWALK)
     &,detuw(MDETUD,MWALK),detdw(MDETUD,MWALK)
     &,ddeti_detiw(3,MELEC,MDETUD,MWALK)

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

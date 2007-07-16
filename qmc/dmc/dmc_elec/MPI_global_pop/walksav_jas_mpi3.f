      subroutine walksav_jas_mpi3
c     subroutine walksav_jas_mpi(iw)
c Written by Claudia Filippi

# if defined (MPI)
      use all_tools_mod
      implicit real*8(a-h,o-z)

      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /jaso/ fso(MELEC,MELEC),fijo(3,MELEC,MELEC)
     &,d2ijo(MELEC,MELEC),d2o,fsumo,fjo(3,MELEC)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eoldw(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      common /mpitype/ jas_type1,jas_type2

      common /jasow/ fsow(MELEC,MELEC,MWALK),fijow(3,MELEC,MELEC,MWALK)
     &,fsumow(MWALK),fjow(3,MELEC,MWALK)
c     dimension istatus(MPI_STATUS_SIZE),idispl(MELEC),iblocklen(MELEC)
      dimension istatus(MPI_STATUS_SIZE)


      entry send_jas_mpi3(irecv)

      itag=0
      call mpi_isend(fsumow(nwalk),1,mpi_double_precision,irecv
     &,itag+1,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
      call mpi_isend(fjow(1,1,nwalk),3*nelec,mpi_double_precision,irecv
     &,itag+2,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
c     do 70 i=1,nelec
c       iblocklen(i)=nelec-(i-1)
c 70    idispl(i)=MELEC*(i-1)+(i-1)
c     call mpi_type_indexed(nelec,iblocklen,idispl,mpi_double_precision
c    &,jas_type1,ierr)
      call mpi_isend(fsow(1,1,nwalk),1,jas_type1,irecv,itag+3
     &,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
c     do 75 i=1,nelec
c       iblocklen(i)=3*nelec
c 75    idispl(i)=3*MELEC*(i-1)
c     call mpi_type_indexed(nelec,iblocklen,idispl,mpi_double_precision
c    &,jas_type2,ierr)
      call mpi_isend(fijow(1,1,1,nwalk),1,jas_type2,irecv,itag+4
     &,MPI_COMM_WORLD,irequest,ierr)
C     RGH
      call MPI_Wait(irequest,istatus,ierr)
c     call mpi_isend(fsow(1,1,nwalk),MELEC*nelec,mpi_double_precision
c    &,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
c     call mpi_isend(fijow(1,1,1,nwalk),3*MELEC*nelec
c    &,mpi_double_precision,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)

      return

      entry recv_jas_mpi3(isend)

      itag=0
      call mpi_recv(fsumow(nwalk),1,mpi_double_precision,isend
     &,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fjow(1,1,nwalk),3*nelec,mpi_double_precision,isend
     &,itag+2,MPI_COMM_WORLD,istatus,ierr)

c     do 80 i=1,nelec
c       iblocklen(i)=nelec-(i-1)
c 80    idispl(i)=MELEC*(i-1)+(i-1)
c     call mpi_type_indexed(nelec,iblocklen,idispl,mpi_double_precision
c    &,jas_type1,ierr)
c     call mpi_type_commit(jas_type1,ierr)
      call mpi_recv(fsow(1,1,nwalk),1,jas_type1,isend,itag+3
     &,MPI_COMM_WORLD,istatus,ierr)

c     do 85 i=1,nelec
c       iblocklen(i)=3*nelec
c 85    idispl(i)=3*MELEC*(i-1)
c     call mpi_type_indexed(nelec,iblocklen,idispl,mpi_double_precision
c    &,jas_type2,ierr)
c     call mpi_type_commit(jas_type2,ierr)
      call mpi_recv(fijow(1,1,1,nwalk),1,jas_type2,isend,itag+4
     &,MPI_COMM_WORLD,istatus,ierr)

c     call mpi_recv(fsow(1,1,nwalk),MELEC*nelec,mpi_double_precision
c    &,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
c     call mpi_recv(fijow(1,1,1,nwalk),3*MELEC*nelec
c    &,mpi_double_precision,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)

# endif
      return
      end

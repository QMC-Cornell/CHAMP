      subroutine dumper_dmc_mov1_mpi2
! MPI version created by Claudia Filippi starting from serial version
! routine to pick up and dump everything needed to restart
! job where it left off

# if defined (MPI)

      use all_tools_mod
      use mpi_mod
      use control_mod
      use atom_mod
      use coefs_mod
      use dets_mod
      use basis1_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use numbas_mod
      use forcepar_mod
      use pseudo_mod
      use contrl_per_mod
      use delocc_mod
      use jaspar_mod
      use force_dmc_mod
      use iterat_mod
      use qua_mod
      use jacobsave_mod
      use forcest_dmc_mod
      use denupdn_mod
      use stepv_mod
      use jaspar1_mod
      use config_dmc_mod
      use branch_mod
      use estsum_dmc_mod
      use estcum_dmc_mod
      use div_v_dmc_mod
      use contrldmc_mod
      use estcm2_mod
      use stats_mod
      use age_mod
      use velratio_mod
      implicit real*8(a-h,o-z)
      parameter (small=1.d-6)

      character*13 filename

      dimension irand_seed_loc(4,0:nproc),istatus(MPI_STATUS_SIZE)
      dimension coefx(nbasis,norb),zexx(nbasis),centx(3,ncent),znucx(ncent) &
     &,n1sx(nctype),n2sx(nctype),n2px(-1:1,nctype) &
     &,n3sx(nctype),n3px(-1:1,nctype),n3dx(-2:2,nctype) &
!    &,n4sx(nctype),n4px(-1:1,nctype),n4dx(-2:2,nctype)
     &,n4sx(nctype),n4px(-1:1,nctype) &
     &,n4fx(-3:3,nctype),n5gx(-4:4,nctype),n6hx(-5:5,nctype) &
     &,nsax(nctype),npax(-1:1,nctype),ndax(-2:2,nctype) &
     &,csf_coefx(ndet)
      integer nprocx

      if(nforce.gt.1) call strech(xoldw,xoldw,ajacob,1,0)

      call savern(irand_seed_loc(1,idtask))

      nscounts=4
      call mpi_gather(irand_seed_loc(1,idtask),nscounts,mpi_integer &
     &,irand_seed_loc,nscounts,mpi_integer,0,MPI_COMM_WORLD,ierr)

      if(idtask.ne.0) then
        call mpi_isend(nwalk,1,mpi_integer,0,1,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(xoldw,3*nelec*nwalk,mpi_double_precision,0,2,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(wt,nwalk,mpi_double_precision,0,3,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(fratio,MWALK*nforce,mpi_double_precision,0,4,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(iage,nwalk,mpi_integer,0,5,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(xq,nquad,mpi_double_precision,0,6,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(yq,nquad,mpi_double_precision,0,7,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(zq,nquad,mpi_double_precision,0,8,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
       else
        open(10,status='unknown',form='unformatted',file='restart_dmc')
        write(10) nproc
        write(10) nfprod,(ff(i),i=0,nfprod),fprod,eigv,eest,wdsumo,ioldest,ioldestmx
        write(10) nwalk
        write(10) (wt(i),i=1,nwalk),(iage(i),i=1,nwalk)
        write(10) (((xoldw(k,i,iw,1),k=1,ndim),i=1,nelec),iw=1,nwalk)
        write(10) nforce,((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
        if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
!       if(nforce.gt.1) write(10) nwprod
!    & ,((pwt(i,j),i=1,nwalk),j=1,nforce)
!    & ,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
        do 450 id=1,nproc-1
          call mpi_recv(nwalk,1,mpi_integer,id,1,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(xoldw,3*nelec*nwalk,mpi_double_precision,id,2,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(wt,nwalk,mpi_double_precision,id,3,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(fratio,MWALK*nforce,mpi_double_precision,id,4,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(iage,nwalk,mpi_integer,id,5,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(xq,nquad,mpi_double_precision,id,6,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(yq,nquad,mpi_double_precision,id,7,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(zq,nquad,mpi_double_precision,id,8,MPI_COMM_WORLD,istatus,ierr)
          write(10) nwalk
          write(10) (wt(i),i=1,nwalk),(iage(i),i=1,nwalk)
          write(10) (((xoldw(k,i,iw,1),k=1,ndim),i=1,nelec),iw=1,nwalk)
          write(10) nforce,((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
          if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
!         if(nforce.gt.1) write(10) nwprod
!    &    ,((pwt(i,j),i=1,nwalk),j=1,nforce)
!    &    ,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
  450   continue
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(idtask.ne.0) return

      write(10) (wgcum(i),egcum(i),pecum(i),tpbcum(i),tjfcum(i) &
     &,wgcm2(i),egcm2(i),pecm2(i),tpbcm2(i),tjfcm2(i),taucum(i) &
     &,i=1,nforce)
      write(10) ((irand_seed_loc(i,j),i=1,4),j=0,nproc-1)
      write(10) hb
      write(10) tau,rttau,idmc
      write(10) nelec,nconf_global
      write(10) (wtgen(i),i=0,nfprod),wgdsumo
      write(10) wcum,wfcum,wdcum,wgdcum &
     &,wcum1,wfcum1,(wgcum1(i),i=1,nforce),wdcum1 &
     &,ecum,efcum,ecum1,efcum1,(egcum1(i),i=1,nforce) &
     &,ei1cum,ei2cum,ei3cum,r2cum,ricum
      write(10) ipass,iblk
      write(10) wcm2,wfcm2,wdcm2,wgdcm2,wcm21 &
     &,wfcm21,(wgcm21(i),i=1,nforce),wdcm21, ecm2,efcm2 &
     &,ecm21,efcm21,(egcm21(i),i=1,nforce) &
     &,ei1cm2,ei2cm2,ei3cm2,r2cm2,ricm2
      write(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
      write(10) (rprob(i)/nproc,rprobup(i),rprobdn(i),i=1,NRAD)
      write(10) dfus2ac,dfus2un,dr2ac,dr2un,acc,acc_int,try_int,nbrnch,nodecr

      write(10) ((coef(ib,i,1),ib=1,nbasis),i=1,norb)
      write(10) nbasis
      write(10) (zex(ib,1),ib=1,nbasis)
      write(10) nctype,ncent,(iwctype(i),i=1,ncent)
      write(10) ((cent(k,i),k=1,ndim),i=1,ncent)
      write(10) pecent
      write(10) (znuc(i),i=1,nctype)
      write(10) (n1s(i),i=1,nctype)
      if(numr.le.0) write(10) (n2s(i),i=1,nctype)
      write(10) ((n2p(m,i),m=-1,1),i=1,nctype)
      if(numr.le.0) then
        write(10) (n3s(i),i=1,nctype)
        write(10) ((n3p(m,i),m=-1,1),i=1,nctype)
      endif
      write(10) ((n3d(m,i),m=-2,2),i=1,nctype)
      if(numr.le.0) then
        write(10) (n4s(i),i=1,nctype)
        write(10) ((n4p(m,i),m=-1,1),i=1,nctype)
!       write(10) ((n4d(m,i),m=-2,2),i=1,nctype)
        write(10) (nsa(i),i=1,nctype)
        write(10) ((npa(m,i),m=-1,1),i=1,nctype)
        write(10) ((nda(m,i),m=-2,2),i=1,nctype)
       else
        write(10) ((n4f(m,i),m=-3,3),i=1,nctype)
        write(10) ((n5g(m,i),m=-4,4),i=1,nctype)
        write(10) ((n6h(m,i),m=-5,5),i=1,nctype)
      endif
      write(10) (csf_coef(i,1),i=1,ncsf)
      write(10) ncsf,ndet,nup,ndn
      write(10) cjas1(1),cjas2(1)
      close (unit=10)
      write(6,'(1x,''successful dump to restart_dmc'')')
      return

      entry startr_dmc_mov1_mpi2
      write(6,'(1x,''attempting restart from restart_dmc'')')
      rewind 10
      read(10) nprocx
      read(10) nfprod,(ff(i),i=0,nfprod),fprod,eigv,eest,wdsumo &
     &,ioldest,ioldestmx
      if(nprocx.ne.nproc) stop 'nproc does not match that in restart file'
      do 4 id=0,idtask
        read(10) nwalk
        read(10) (wt(i),i=1,nwalk),(iage(i),i=1,nwalk)
        read(10) (((xoldw(k,i,iw,1),k=1,ndim),i=1,nelec),iw=1,nwalk)
        read(10) nforce,((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
    4   if(nloc.gt.0) &
     &  read(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
      do 5 id=idtask+1,nproc-1
        read(10) nwalk_id
        read(10) (wt_id,i=1,nwalk_id),(iage_id,i=1,nwalk_id)
        read(10) (xold_id,i=1,ndim*nelec*nwalk_id)
        read(10) nf_id,((fratio_id,iw=1,nwalk_id),ifr=1,nf_id)
    5   if(nloc.gt.0) &
     &  read(10) nq_id,(xq_id,yq_id,zq_id,wq_id,i=1,nquad)
!     if(nforce.gt.1) read(10) nwprod
!    &,((pwt(i,j),i=1,nwalk),j=1,nforce)
!    &,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
      call alloc ('wgcm2', wgcm2, nforce)
      call alloc ('egcm2', egcm2, nforce)
      call alloc ('pecm2', pecm2, nforce)
      call alloc ('tpbcm2', tpbcm2, nforce)
      call alloc ('tjfcm2', tjfcm2, nforce)
      call alloc ('wgcum', wgcum, nforce)
      call alloc ('egcum', egcum, nforce)
      call alloc ('pecum', pecum, nforce)
      call alloc ('peicum', peicum, nforce)
      call alloc ('tpbcum', tpbcum, nforce)
      call alloc ('tjfcum', tjfcum, nforce)
      call alloc ('taucum', taucum, nforce)
      read(10) (wgcum(i),egcum(i),pecum(i),tpbcum(i),tjfcum(i), &
     &wgcm2(i),egcm2(i),pecm2(i),tpbcm2(i),tjfcm2(i),taucum(i), &
     &i=1,nforce)
      read(10) ((irand_seed_loc(i,j),i=1,4),j=0,nproc-1)
      call setrn(irand_seed_loc(1,idtask))
      read(10) hbx
      read(10) taux,rttau,idmc
      read(10) nelecx,nconf_global
      if(dabs(hbx-hb).gt.small) stop 'hb'
      if(dabs(taux-tau).gt.small) stop 'tau'
      if(nelecx.ne.nelec) stop 'nelec'
      read(10) (wtgen(i),i=0,nfprod),wgdsumo
      call alloc ('wgcum1', wgcum1, nforce)
      call alloc ('egcum1', egcum1, nforce)
      read(10) wcum,wfcum,wdcum,wgdcum,wcum1 &
     &,wfcum1,(wgcum1(i),i=1,nforce),wdcum1, ecum,efcum &
     &,ecum1,efcum1,(egcum1(i),i=1,nforce) &
     &,ei1cum,ei2cum,ei3cum, r2cum,ricum
      read(10) ipass,iblk
      call alloc ('wgcm21', wgcm21, nforce)
      call alloc ('egcm21', egcm21, nforce)
      read(10) wcm2,wfcm2,wdcm2,wgdcm2,wcm21 &
     &,wfcm21,(wgcm21(i),i=1,nforce),wdcm21, ecm2,efcm2 &
     &,ecm21,efcm21,(egcm21(i),i=1,nforce) &
     &,ei1cm2,ei2cm2,ei3cm2,r2cm2,ricm2
      call alloc ('fgcum', fgcum, nforce)
      call alloc ('fgcm2', fgcm2, nforce)
      read(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
      call alloc ('rprob', rprob, NRAD)
      call alloc ('rprobup', rprobup, NRAD)
      call alloc ('rprobdn', rprobdn, NRAD)
      read(10) (rprob(i),rprobup(i),rprobdn(i),i=1,NRAD)
      read(10) dfus2ac,dfus2un,dr2ac,dr2un,acc,acc_int,try_int,nbrnch,nodecr
      if(idtask.ne.0) then
        acc=0
        acc_int=0
        try_int=0
        dr2ac=0
        dr2un=0
        nodecr=0
      endif

      read(10) ((coefx(ib,i),ib=1,nbasis),i=1,norb)
      read(10) nbasx
      do 10 j=1,norb
      do 10 i=1,nbasis
      if(dabs(coefx(i,j)-coef(i,j,1)).gt.small) stop 'coef'
   10 continue
      if(nbasx.ne.nbasis) stop 'nbasis'
      read(10) (zexx(ib),ib=1,nbasis)
      read(10) nctypex,ncentx,(iwctype(i),i=1,ncentx)
      read(10) ((centx(k,i),k=1,ndim),i=1,ncentx)
      read(10) pecent
      read(10) (znucx(i),i=1,nctypex)
      read(10) (n1sx(i),i=1,nctype)
      if(numr.le.0) read(10) (n2sx(i),i=1,nctype)
      read(10) ((n2px(m,i),m=-1,1),i=1,nctype)
      if(numr.le.0) then
        read(10) (n3sx(i),i=1,nctype)
        read(10) ((n3px(m,i),m=-1,1),i=1,nctype)
      endif
      read(10) ((n3dx(m,i),m=-2,2),i=1,nctype)
      if(numr.le.0) then
        read(10) (n4sx(i),i=1,nctype)
        read(10) ((n4px(m,i),m=-1,1),i=1,nctype)
!       read(10) ((n4dx(m,i),m=-2,2),i=1,nctype)
        read(10) (nsax(i),i=1,nctype)
        read(10) ((npax(m,i),m=-1,1),i=1,nctype)
        read(10) ((ndax(m,i),m=-2,2),i=1,nctype)
       else
        read(10) ((n4fx(m,i),m=-3,3),i=1,nctype)
        read(10) ((n5gx(m,i),m=-4,4),i=1,nctype)
        read(10) ((n6hx(m,i),m=-5,5),i=1,nctype)
      endif

      if(ncentx.ne.ncent) stop 'ncent'
      if(nctypex.ne.nctype) stop 'nctype'
      do 20 i=1,nbasis
      if(dabs(zexx(i)-zex(i,1)).gt.small) stop 'zex'
   20 continue
      do 30 i=1,ncent
      do 30 k=1,ndim
      if(dabs(cent(k,i)-centx(k,i)).gt.small) stop 'cent'
   30 continue
      do 49 i=1,nctype
        if(dabs(znucx(i)-znuc(i)).gt.small) stop 'znuc'
        if(n1s(i).ne.n1sx(i)) stop 'n1s'
        if(numr.le.0) then
          if(n2s(i).ne.n2sx(i)) stop 'n2s'
          if(n3s(i).ne.n3sx(i)) stop 'n3s'
          if(n4s(i).ne.n4sx(i)) stop 'n4s'
          if(nsa(i).ne.nsax(i)) stop 'nsa'
        endif
        do 40 m=-1,1
          if(n2p(m,i).ne.n2px(m,i)) stop 'n2p'
          if(numr.le.0) then
            if(n3p(m,i).ne.n3px(m,i)) stop 'n3p'
            if(n4p(m,i).ne.n4px(m,i)) stop 'n4p'
            if(npa(m,i).ne.npax(m,i)) stop 'npa'
          endif
   40   continue
        do 45 m=-2,2
          if(n3d(m,i).ne.n3dx(m,i)) stop 'n3d'
          if(numr.le.0) then
!           if(n4d(m,i).ne.n4dx(m,i)) stop 'n4d'
            if(nda(m,i).ne.ndax(m,i)) stop 'nda'
          endif
   45   continue
        if(numr.gt.0) then
          do 46 m=-3,3
   46       if(n4f(m,i).ne.n4fx(m,i)) stop 'n4f'
          do 47 m=-4,4
   47       if(n5g(m,i).ne.n5gx(m,i)) stop 'n5g'
          do 48 m=-5,5
   48       if(n6h(m,i).ne.n6hx(m,i)) stop 'n6h'
        endif
   49 continue
      read(10) (csf_coefx(i),i=1,ncsf)
      read(10) ncsfx,ndetx,nupx,ndnx
      do 50 i=1,ncsf
      if(dabs(csf_coefx(i)-csf_coef(i,1)).gt.small) stop 'csf_coef'
   50 continue
      if(ncsfx.ne.ncsf) stop 'ncsf'
      if(ndetx.ne.ndet) stop 'ndet'
      if(nupx.ne.nup) stop 'nup'
      if(ndnx.ne.ndn) stop 'ndn'
      read(10) cjas1x,cjas2x
      if(dabs(cjas1x-cjas1(1)).gt.small) stop 'cjas1'
      if(dabs(cjas2x-cjas2(1)).gt.small) stop 'cjas2'
      write(6,'(1x,''succesful read from restart_dmc'')')
      write(6,'(t5,''egnow'',t15,''egave'',t21 &
     &,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55 &
     &,''(tpberr)'' ,t66,''tjfave'',t72,''(tjferr)'',t83,''npass'',t93 &
     &,''wgsum'',t103 ,''ioldest'')')

      do 70 iw=1,nwalk
        if(istrech.eq.0) then
          do 60 ifr=2,nforce
            do 60 ie=1,nelec
              do 60 k=1,ndim
   60           xoldw(k,ie,iw,ifr)=xoldw(k,ie,iw,1)
        endif
        do 70 ifr=1,nforce
          if(nforce.gt.1) then
            if(ifr.eq.1.or.istrech.eq.0) then
              call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,ifr),ajacob,ifr,0)
               else
              call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,ifr),ajacob,ifr,1)
            endif
           else
            ajacob=one
          endif
          ajacold(iw,ifr)=ajacob
          call hpsi(xoldw(1,1,iw,ifr),psidow(iw,ifr),psijow(iw,ifr),voldw(1,1,iw,ifr),div_vow(1,iw),d2ow(iw,ifr), &
     &    peow(iw,ifr),peiow(iw,ifr),eoldw(iw,ifr),denergy,ifr)
          if(ifr.eq.1) then
            if(ibasis.eq.3) then                ! complex calculation
              call cwalksav_det(iw)
            else
              call walksav_det(iw)
            endif
            call walksav_jas(iw)
          endif
   70 continue

! zero out xsum variables

      wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r2sum=zero
      risum=zero

      call alloc ('wgsum', wgsum, nforce)
      call alloc ('egsum', egsum, nforce)
      call alloc ('pesum', pesum, nforce)
      call alloc ('peisum', peisum, nforce)
      call alloc ('tpbsum', tpbsum, nforce)
      call alloc ('tjfsum', tjfsum, nforce)
      call alloc ('tausum', tausum, nforce)
      do 80 ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum(ifr)=zero
        peisum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
   80   tausum(ifr)=zero

      if(ipr.gt.-2) then
        if(idtask.le.9) then
          write(filename,'(''walkalize.'',i1)') idtask
         elseif(idtask.le.99) then
          write(filename,'(''walkalize.'',i2)') idtask
         elseif(idtask.le.999) then
          write(filename,'(''walkalize.'',i3)') idtask
         else
          stop 'idtask > 999'
        endif
        open(11,file=filename,status='old')
        do 90 i=1,2000000000
   90     read(11,fmt=*,end=100)
      endif
  100 backspace 11
      backspace 11

# endif
      return
      end

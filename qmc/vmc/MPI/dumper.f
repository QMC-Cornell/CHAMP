      subroutine dumper_mpi
c MPI version created by Claudia Filippi starting from Cyrus' serial version.
c Routine to pick up and dump everything needed to restart job where it left off.

# if defined (MPI)

      use all_tools_mod
      use atom_mod
      use config_mod
      use coefs_mod
      use dets_mod
      use basis1_mod
      use const_mod
      use dim_mod
      use numbas_mod
      use forcepar_mod
      use wfsec_mod
      use doefp_mod
      use pseudo_mod
      use delocc_mod
      use jaspar_mod
      use div_v_mod
      use qua_mod
      use forcest_mod
      use denupdn_mod
      use stepv_mod
      use jaspar1_mod
      implicit real*8(a-h,o-z)

      parameter(small=1.d-6)

!JT      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
c     common /contrl/ nstep,nblk,nblkeq,nconf,nconf_global,nconf_new,isite,idump,irstar
!JT      common /config/ xold(3,MELEC),xnew(3,MELEC),vold(3,MELEC)
!JT     &,vnew(3,MELEC),psi2o(MFORCE),psi2n(MFORCE),eold(MFORCE),enew(MFORCE)
!JT     &,peo,pen,peio,pein,tjfn,tjfo,psido,psijo
!JT     &,rmino(MELEC),rminn(MELEC),rvmino(3,MELEC),rvminn(3,MELEC)
!JT     &,rminon(MELEC),rminno(MELEC),rvminon(3,MELEC),rvminno(3,MELEC)
!JT     &,nearesto(MELEC),nearestn(MELEC),delttn(MELEC)
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere
!JT      common /basis/ zex(MBASIS,MWF),betaq
!JT     &,n1s(MCTYPE)
!JT     &,n2s(MCTYPE),n2p(-1:1,MCTYPE)
!JT     &,n3s(MCTYPE),n3p(-1:1,MCTYPE),n3d(-2:2,MCTYPE)
!JT     &,n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE),n4f(-3:3,MCTYPE)
!JT     &,n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE),n5f(-3:3,MCTYPE)
!JT     &,n5g(-4:4,MCTYPE)
!JT     &,n6d(-2:2,MCTYPE),n6f(-3:3,MCTYPE),n6g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
!JT     &,n7g(-4:4,MCTYPE),n7h(-5:5,MCTYPE),n7i(-6:6,MCTYPE)
!JT     &,n8i(-6:6,MCTYPE),n8j(-7:7,MCTYPE)
!JT     &,n9k(-8:8,MCTYPE)
!JT     &,n10l(-9:9,MCTYPE)
!JT     &,n11m(-10:10,MCTYPE)
!JT     &,n12n(-11:11,MCTYPE)
!JT     &,n13o(-12:12,MCTYPE)
!JT     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
!JT      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
!JT     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
!JT     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)
!JT      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
!JT     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
!JT      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
!JT     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
!JT      common /jaspar1/ cjas1(MWF),cjas2(MWF)
!JT      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
!JT      common /delocc/ denergy(MPARM)
      common /estsum/ esum1,esum(MFORCE),pesum,peisum,tpbsum,tjfsum,r2sum,accsum
      common /estcum/ ecum1,ecum(MFORCE),pecum,peicum,tpbcum,tjfcum,r2cum,acccum,iblk
      common /est2cm/ ecm21,ecm2(MFORCE),pecm2,peicm2,tpbcm2,tjfcm2,r2cm2
      common /estsig/ wsum1s(MFORCE),esum1s(MFORCE),ecum1s(MFORCE),ecm21s(MFORCE)
      common /stats_vmc/ rejmax
!JT      common /stepv/ try(NRAD),suc(NRAD),trunfb(NRAD),rprob(NRAD),
!JT     &ekin(NRAD),ekin2(NRAD)
!JT      common /denupdn/ rprobup(NRAD),rprobdn(NRAD)
!JT      common /wfsec/ iwftype(MFORCE),iwf,nwftype
!JT      common /forcepar/ deltot(MFORCE),nforce,istrech
!JT      common /forcest/ fcum(MFORCE),fcm2(MFORCE)
      common /forcewt/ wsum(MFORCE),wcum(MFORCE)
      common /forcjac/ ajacob
!JT      common /doefp/ nefp
!JT      common /div_v/ div_vo(MELEC)

      dimension irn(4,0:MPROC),istatus(MPI_STATUS_SIZE)
      dimension ircounts(0:MPROC),idispls(0:MPROC)

      dimension coefx(MBASIS,MORB),zexx(MBASIS),centx(3,MCENT),znucx(MCTYPE)
     &,n1sx(MCTYPE),n2sx(MCTYPE),n2px(-1:1,MCTYPE)
     &,n3sx(MCTYPE),n3px(-1:1,MCTYPE),n3dx(-2:2,MCTYPE)
c    &,n4sx(MCTYPE),n4px(-1:1,MCTYPE),n4dx(-2:2,MCTYPE)
     &,n4sx(MCTYPE),n4px(-1:1,MCTYPE)
     &,n4fx(-3:3,MCTYPE),n5gx(-4:4,MCTYPE),n6hx(-5:5,MCTYPE)
     &,nsax(MCTYPE),npax(-1:1,MCTYPE),ndax(-2:2,MCTYPE)
     &,csf_coefx(MDET)
      dimension xstrech(3,MELEC)

      rewind 10
      if(idtask.eq.0) write(10) nproc
      call savern(irn(1,idtask))

      do 124 i=0,nproc-1
        ircounts(i)=4
  124   idispls(i)=i*4
      idispls(nproc)=4*nproc

      nscounts=ircounts(idtask)
      call mpi_gatherv(irn(1,idtask),nscounts,mpi_integer
     &,irn,ircounts,idispls,mpi_integer,0,MPI_COMM_WORLD,ierr)

      if(idtask.ne.0) then
        call mpi_isend(xold,3*nelec,mpi_double_precision,0
     &  ,1,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(xq,nquad,mpi_double_precision,0
     &  ,2,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(yq,nquad,mpi_double_precision,0
     &  ,3,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
        call mpi_isend(zq,nquad,mpi_double_precision,0
     &  ,4,MPI_COMM_WORLD,irequest,ierr)
         call MPI_Wait(irequest,istatus,ierr)
       else
        write(10) ((xold(k,i),k=1,ndim),i=1,nelec)
        if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
        do 450 id=1,nproc-1
          call mpi_recv(xold,3*nelec,mpi_double_precision,id
     &    ,1,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(xq,nquad,mpi_double_precision,id
     &    ,2,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(yq,nquad,mpi_double_precision,id
     &    ,3,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(zq,nquad,mpi_double_precision,id
     &    ,4,MPI_COMM_WORLD,istatus,ierr)
          write(10) ((xold(k,i),k=1,ndim),i=1,nelec)
  450     if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(idtask.ne.0) return
      write(10) ((irn(i,j),i=1,4),j=0,nproc-1)
      write(10) hb,delta
      write(10) nelec,nforce
      write(10) ((rvmino(k,i),k=1,ndim),i=1,nelec)
      write(10) (rmino(i),i=1,nelec)
      write(10) (nearesto(i),i=1,nelec)
      write(10) (psi2o(i),eold(i),i=1,nforce),peo,tjfo
c     write(10) ecum1/nproc,(ecum(i),i=1,nforce),pecum,tpbcum,
      write(10) ecum1,(ecum(i),i=1,nforce),pecum,tpbcum,
     &tjfcum,r2cum,acccum
      write(10) iblk
c     write(10) ecm21/nproc,(ecm2(i),i=1,nforce),pecm2,tpbcm2,
      write(10) ecm21,(ecm2(i),i=1,nforce),pecm2,tpbcm2,
     &tjfcm2,r2cm2
      if(nforce.gt.0) write(10) (wcum(i),fcum(i),fcm2(i),i=1,nforce)
      do 3 i=1,nforce
c   3   write(10) wsum1s(i),esum1s(i),ecum1s(i),ecm21s(i)
    3   write(10) ecum1s(i),ecm21s(i)

      write(10) rejmax
      write(10) (try(i)/nproc,suc(i)/nproc,trunfb(i)/nproc,
     &rprob(i)/nproc,rprobup(i)/nproc,rprobdn(i)/nproc,
     &ekin(i)/nproc,ekin2(i)/nproc,i=1,NRAD)
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
c       write(10) ((n4d(m,i),m=-2,2),i=1,nctype)
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
      rewind 10
      write(6,'(1x,''successful dump to restart_vmc'')')

      return
c-----------------------------------------------------------------------

      entry startr_mpi
      write(6,'(1x,''attempting restart from restart_vmc'')')

      call pot_nn(cent,znuc,iwctype,ncent,pecent)
c     pecent=0
c     if(ncent.lt.2) goto 3
c     do 2 i=2,ncent
c       j1=i-1
c       do 2 j=1,j1
c         r2=0
c         do 1 k=1,ndim
c   1       r2=r2+(cent(k,i)-cent(k,j))**2
c       r=dsqrt(r2)
c   2   pecent=pecent+znuc(iwctype(i))*znuc(iwctype(j))/r
c   3 continue

      rewind 10
      read(10) nprocx
      if(nprocx.ne.nproc) stop 'nproc does not match that in restart file'
      do 4 id=0,idtask
        read(10) ((xold(k,i),k=1,ndim),i=1,nelec)
    4   if(nloc.gt.0) read(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
      do 5 id=idtask+1,nproc-1
        read(10) (x_id,i=1,ndim*nelec)
    5   if(nloc.gt.0) read(10) nq_id,(xq_id,yq_id,zq_id,wq_id,i=1,nq_id)
      read(10) ((irn(i,j),i=1,4),j=0,nproc-1)
      call setrn(irn(1,idtask))
      read(10) hbx,deltax
      read(10) nelecx,nforce
      if(dabs(hbx-hb).gt.small) stop 'hb'
      if(dabs(deltax-delta).gt.small) stop 'delta'
      if(nelecx.ne.nelec) stop 'nelec'
      read(10) ((rvmino(k,i),k=1,ndim),i=1,nelec)
      read(10) (rmino(i),i=1,nelec)
      read(10) (nearesto(i),i=1,nelec)
      read(10) (psi2o(i),eold(i),i=1,nforce),peo,tjfo
      read(10) ecum1,(ecum(i),i=1,nforce),pecum,tpbcum,
     &tjfcum,r2cum,acccum
      read(10) iblk
      read(10) ecm21,(ecm2(i),i=1,nforce),pecm2,tpbcm2,
     &tjfcm2,r2cm2
      if(nforce.gt.0) read(10) (wcum(i),fcum(i),fcm2(i),i=1,nforce)
      do 6 i=1,nforce
c   6   read(10) wsum1s(i),esum1s(i),ecum1s(i),ecm21s(i)
    6   read(10) ecum1s(i),ecm21s(i)

      read(10) rejmax
      read(10) (try(i),suc(i),trunfb(i),rprob(i),
     &rprobup(i),rprobdn(i),
     &ekin(i),ekin2(i),i=1,NRAD)
      read(10) ((coefx(ib,i),ib=1,nbasis),i=1,norb)
      read(10) nbasx
      if(nbasx.ne.nbasis) stop 'nbasis'
      do 10 j=1,norb
        do 10 i=1,nbasis
          if(dabs(coefx(i,j)-coef(i,j,1)).gt.small) stop 'coef'
   10 continue
      read(10) (zexx(ib),ib=1,nbasis)
      read(10) nctypex,ncentx,(iwctype(i),i=1,ncentx)
      read(10) ((centx(k,i),k=1,ndim),i=1,ncentx)
      read(10) pecx
      if(ncentx.ne.ncent) stop 'ncent'
      if(nctypex.ne.nctype) stop 'nctype'
      read(10) (znucx(i),i=1,nctype)
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
c       read(10) ((n4dx(m,i),m=-2,2),i=1,nctype)
        read(10) (nsax(i),i=1,nctype)
        read(10) ((npax(m,i),m=-1,1),i=1,nctype)
        read(10) ((ndax(m,i),m=-2,2),i=1,nctype)
       else
        read(10) ((n4fx(m,i),m=-3,3),i=1,nctype)
        read(10) ((n5gx(m,i),m=-4,4),i=1,nctype)
        read(10) ((n6hx(m,i),m=-5,5),i=1,nctype)
      endif
      do 20 i=1,nbasis
        if(dabs(zexx(i)-zex(i,1)).gt.small) stop 'zex'
   20 continue
      do 30 i=1,ncent
        do 30 k=1,ndim
          if(dabs(cent(k,i)-centx(k,i)).gt.small) stop 'cent'
   30 continue
      if(pecx.ne.pecent) stop 'pec'
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
c           if(n4d(m,i).ne.n4dx(m,i)) stop 'n4d'
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
      write(6,'(1x,''succesful read from restart_vmc'')')
      write(6,'(t5,''enow'',t15,''eave'',t25,''eerr'',t35,''peave'',
     &t45,''peerr'',t55,''tpbave'',t65,''tpberr'',t75,''tjfave'',
     &t85,''tjferr'',t95,''accept'',t105,''iter'')')

      if(nforce.gt.1) then
        call readforce
        call wf_secondary
       else
        nwftype=1
        iwftype(1)=1
      endif

c loop over secondary config
      do 80 ifr=2,nforce
c set n- and e-coord and n-n potential
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido,psijo,vold,div_vo,d2,peo,peio,eold(ifr),denergy,ifr)
   80   psi2o(ifr)=2*(dlog(dabs(psido))+psijo)+dlog(ajacob)

c primary config
c set n-coord and n-n potential
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido,psijo,vold,div_vo,d2,peo,peio,eold(1),denergy,1)
      psi2o(1)=2*(dlog(dabs(psido))+psijo)
      tjfo=d2
      tjfo=-tjfo*half*hb

      if(nefp.gt.0) then
        call restartefp
        call sample_efp(1,xold,eold(1),0.d0)
        call efpsav
      endif

      do 400 i=1,nelec
        do 400 k=1,ndim
  400     xnew(k,i)=xold(k,i)

      do 86 i=1,nelec
        rmino(i)=99.d9
        do 85 j=1,ncent
          dist=0
          do 84 k=1,ndim
   84       dist=dist+(xold(k,i)-cent(k,j))**2
          if(dist.lt.rmino(i)) then
            rmino(i)=dist
            nearesto(i)=j
          endif
   85     continue
        rmino(i)=dsqrt(rmino(i))
        do 86 k=1,ndim
   86     rvmino(k,i)=xold(k,i)-cent(k,nearesto(i))

      do 87 ifr=1,nforce
        esum1s(ifr)=0
        wsum1s(ifr)=0
        esum(ifr)=0
   87   wsum(ifr)=0
      pesum=0
      peisum=0
      tpbsum=0
      tjfsum=0
      r2sum=0

c Zero out all except master processor
      if(idtask.eq.0) return

      acccum=0
      pecum=0
      peicum=0
      tpbcum=0
      tjfcum=0
      r2cum=0

      pecm2=0
      peicm2=0
      tpbcm2=0
      tjfcm2=0
      r2cm2=0

      ecum1=0
      ecm21=0

      do 180 ifr=1,nforce
        ecum1s(ifr)=0
        ecm21s(ifr)=0
        ecum(ifr)=0
        ecm2(ifr)=0
        wcum(ifr)=0
        fcum(ifr)=0
  180   fcm2(ifr)=0

      do 190 i=1,NRAD
        try(i)=0
        suc(i)=0
        trunfb(i)=0
        ekin(i)=0
        ekin2(i)=0
  190   rprob(i)=0

c     call prop_init(0)
c     call efpci_init(0)
c     call eh_ex_init(0)
      call grad_hess_jas_init

# endif
      return
      end

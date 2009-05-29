      subroutine dumper
#     if defined (MPI)
        call dumper_mpi
#     else
        call dumper_serial
#     endif
      end
      subroutine startr
#     if defined (MPI)
        call startr_mpi
#     else
        call startr_serial
#     endif
      end

      subroutine dumper_serial
c Written by Cyrus Umrigar, modified by Claudia Filippi
c Routine to pick up and dump everything needed to restart job where it left off.
      use constants_mod
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
!      use doefp_mod
      use pseudo_mod
      use delocc_mod
      use div_v_mod
      use qua_mod
      use forcest_mod
      use denupdn_mod
      use stepv_mod
      use jaspar1_mod
      implicit real*8(a-h,o-z)

      parameter(small=1.d-6)

      common /estsum/ esum1,esum(MFORCE),pesum,peisum,tpbsum,tjfsum,r2sum,accsum
      common /estcum/ ecum1,ecum(MFORCE),pecum,peicum,tpbcum,tjfcum,r2cum,acccum,iblk
      common /est2cm/ ecm21,ecm2(MFORCE),pecm2,peicm2,tpbcm2,tjfcm2,r2cm2
      common /estsig/ wsum1s(MFORCE),esum1s(MFORCE),ecum1s(MFORCE),ecm21s(MFORCE)
      common /stats_vmc/ rejmax
      common /forcewt/ wsum(MFORCE),wcum(MFORCE)

      dimension irn(4)
      dimension coefx(nbasis,norb),zexx(nbasis),centx(3,ncent),znucx(nctype)
     &,n1sx(nctype),n2sx(nctype),n2px(-1:1,nctype)
     &,n3sx(nctype),n3px(-1:1,nctype),n3dx(-2:2,nctype)
c    &,n4sx(nctype),n4px(-1:1,nctype),n4dx(-2:2,nctype)
     &,n4sx(nctype),n4px(-1:1,nctype)
     &,n4fx(-3:3,nctype),n5gx(-4:4,nctype),n6hx(-5:5,nctype)
     &,nsax(nctype),npax(-1:1,nctype),ndax(-2:2,nctype)
     &,csf_coefx(MDET)
      dimension xstrech(3,MELEC)

      call savern(irn)
      rewind 10
      write(10) irn
      write(10) hb,delta
      write(10) nelec,nforce
      write(10) ((xold(k,i),k=1,ndim),i=1,nelec)
      write(10) ((vold(k,i),k=1,ndim),i=1,nelec)
      write(10) ((rvmino(k,i),k=1,ndim),i=1,nelec)
      write(10) (rmino(i),i=1,nelec)
      write(10) (nearesto(i),i=1,nelec)
      write(10) (psi2o(i),eold(i),i=1,nforce),peo,tjfo
      write(10) ecum1,(ecum(i),i=1,nforce),pecum,tpbcum,
     &tjfcum,r2cum,acccum
      write(10) iblk
      write(10) ecm21,(ecm2(i),i=1,nforce),pecm2,tpbcm2,
     &tjfcm2,r2cm2
      if(nforce.gt.0) write(10) (wcum(i),fcum(i),fcm2(i),i=1,nforce)
      do 3 i=1,nforce
c   3   write(10) wsum1s(i),esum1s(i),ecum1s(i),ecm21s(i)
    3   write(10) ecum1s(i),ecm21s(i)
      write(10) rejmax
      write(10) (try(i),suc(i),trunfb(i),rprob(i),
     &rprobup(i),rprobdn(i),
     &ekin(i),ekin2(i),i=1,NRAD)
      if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
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

      call grad_hess_jas_dump(10)

      rewind 10
      write(6,'(1x,''successful dump to restart_vmc'')')

      call grad_hess_jas_dump(10)

      return

      entry startr_serial
      write(6,'(1x,''attempting restart from restart_vmc'')')

c     call pot_nn(cent,znuc,iwctype,ncent,pecent)

      rewind 10
      read(10) irn
      call setrn(irn)
      read(10) hbx,deltax
      read(10) nelecx,nforce
      if(dabs(hbx-hb).gt.small) stop 'hb'
      if(dabs(deltax-delta).gt.small) stop 'delta'
      if(nelecx.ne.nelec) stop 'nelec'
      read(10) ((xold(k,i),k=1,ndim),i=1,nelec)
      read(10) ((vold(k,i),k=1,ndim),i=1,nelec)
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
      if(nloc.gt.0) read(10) nqx,(xq(i),yq(i),zq(i),wq(i),i=1,nqx)
      if(nqx.ne.nquad) stop 'nquad'
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

      call grad_hess_jas_rstrt(10)

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

!      if(nefp.gt.0) then
!        call restartefp
!        call sample_efp(1,xold,eold(1),0.d0)
!        call efpsav
!      endif

      call grad_hess_jas_save

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

      return
      end

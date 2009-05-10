      subroutine cwalksav_det(iw)
c Same as walksav_det, adapted to complex orbitals
c by A.D.Guclu, Feb2004.

      use basic_tools_mod
      use cslater_mod
      use cslaterw_mod
      use dorb_mod

      implicit real*8(a-h,o-z)

c complex locals:
c      complex*16 cslmuiw,cslmdiw,cfpuw,cfpdw,cdetuw,cdetdw,cddeti_detiw

c complex commons:
c      complex*16 cslmui,cslmdi,cfpu,cfpd,cfppu,cfppd,cdetu,cdetd,cddeti_deti,cd2edeti_deti
c      complex*16 cdeti_det,cddeti_det,cd2deti_det,cd2det_det

c following commons don't seem to be used?
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eoldw(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
!JT      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
!JT     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
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

      do 45 idet=1,ndetup
        cdetuw(idet,iw)=cdetu(idet)
        do 40 j=1,nup*nup
          cslmuiw(j,idet,iw)=cslmui(j,idet)
          do 40 idim=1,ndim
   40       cfpuw(idim,j,idet,iw)=cfpu(idim,j,idet)
          do 45 i=1,nup
            do 45 idim=1,ndim
   45         cddeti_detiw(idim,i,idet,iw)=cddeti_deti(idim,i,idet)

      do 55 idet=1,ndetdn
        cdetdw(idet,iw)=cdetd(idet)
        do 50 j=1,ndn*ndn
          cslmdiw(j,idet,iw)=cslmdi(j,idet)
          do 50 idim=1,ndim
   50       cfpdw(idim,j,idet,iw)=cfpd(idim,j,idet)
          do 55 i=nup+1,nelec
            do 55 idim=1,ndim
   55         cddeti_detiw(idim,i,idet,iw)=cddeti_deti(idim,i,idet)

      return
      entry cwalkstrdet(iw)

       do 85 idet=1,ndetup
         cdetu(idet)=cdetuw(idet,iw)
         do 80 j=1,nup*nup
           cslmui(j,idet)=cslmuiw(j,idet,iw)
           do 80 idim=1,ndim
   80        cfpu(idim,j,idet)=cfpuw(idim,j,idet,iw)
         do 85 i=1,nup
           do 85 idim=1,ndim
   85        cddeti_deti(idim,i,idet)=cddeti_detiw(idim,i,idet,iw)

       do 95 idet=1,ndetdn
         cdetd(idet)=cdetdw(idet,iw)
         do 90 j=1,ndn*ndn
           cslmdi(j,idet)=cslmdiw(j,idet,iw)
           do 90 idim=1,ndim
   90        cfpd(idim,j,idet)=cfpdw(idim,j,idet,iw)
         do 95 i=nup+1,nelec
           do 95 idim=1,ndim
   95        cddeti_deti(idim,i,idet)=cddeti_detiw(idim,i,idet,iw)

      return
      entry csplitjdet(iw,iw2)

      do 125 idet=1,ndetup
        cdetuw(idet,iw2)=cdetuw(idet,iw)
        do 120 j=1,nup*nup
          cslmuiw(j,idet,iw2)=cslmuiw(j,idet,iw)
          do 120 idim=1,ndim
  120       cfpuw(idim,j,idet,iw2)=cfpuw(idim,j,idet,iw)
        do 125 i=1,nup
          do 125 idim=1,ndim
  125       cddeti_detiw(idim,i,idet,iw2)=cddeti_detiw(idim,i,idet,iw)

      do 135 idet=1,ndetdn
        cdetdw(idet,iw2)=cdetdw(idet,iw)
        do 130 j=1,ndn*ndn
          cslmdiw(j,idet,iw2)=cslmdiw(j,idet,iw)
          do 130 idim=1,ndim
  130       cfpdw(idim,j,idet,iw2)=cfpdw(idim,j,idet,iw)
        do 135 i=nup+1,nelec
          do 135 idim=1,ndim
  135       cddeti_detiw(idim,i,idet,iw2)=cddeti_detiw(idim,i,idet,iw)

      return
      end

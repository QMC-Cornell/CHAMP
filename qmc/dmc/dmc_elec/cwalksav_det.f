      subroutine cwalksav_det(iw)
! Same as walksav_det, adapted to complex orbitals
! by A.D.Guclu, Feb2004.

      use basic_tools_mod
      use cslater_mod
      use cslaterw_mod
      use dorb_mod
      use dets_mod
      use const_mod
      use dim_mod
      use forcepar_mod
      use force_dmc_mod
      use branch_mod
      implicit real*8(a-h,o-z)

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

      subroutine cdetsav(iel)
c same as detsav
c adapted to complex calculations by A.D.Guclu Feb2004

      use basic_tools_mod
      use cslater_mod
      use dorb_mod
      use cslatn_mod
      use dets_mod
      use const_mod
      use dim_mod
      implicit real*8(a-h,o-z)

c complex commons:
c      complex*16 cslmui,cslmdi,cfpu,cfpd,cfppu,cfppd,cdetu,cdetd,cddeti_deti,cd2edeti_deti
c      complex*16 cdeti_det,cddeti_det,cd2deti_det,cd2det_det
!JT      complex*16 cslmin,cdetn,cddeti_detin,cd2edeti_detin,cdorb,cddorb

c      common /cslater/ cslmui(MMAT_DIM,MDET),cslmdi(MMAT_DIM,MDET)
c     &,cfpu(3,MMAT_DIM,MDET),cfpd(3,MMAT_DIM,MDET)
c     &,cfppu(MMAT_DIM,MDET),cfppd(MMAT_DIM,MDET)
c     &,cdetu(MDET),cdetd(MDET)
c     &,cddeti_deti(3,MELEC,MDET),cd2edeti_deti(MELEC,MDET),cdeti_det(MCSF),cddeti_det(3,MELEC,MCSF),cd2deti_det(MCSF),cd2det_det
c cd2edeti_detin and cddorb dont need to be in common?

c allocate memory (will allocate only if it is not already allocated):
      n2=nelec*nelec
      call alloc('cslmui',cslmui,n2,ndetup)
      call alloc('cslmdi',cslmdi,n2,ndetdn)
      call alloc('cfpu',cfpu,ndim,n2,ndetup)
      call alloc('cfpd',cfpd,ndim,n2,ndetdn)
      call alloc('cdetu',cdetu,ndetup)
      call alloc('cdetd',cdetd,ndetdn)
      call alloc('cddeti_deti',cddeti_deti,ndim,nelec,ndet)



      if(iel.le.nup) then
        ikel=nup*(iel-1)
        do 60 idet=1,ndetup
          cdetu(idet)=cdetn(idet)
          do 45 l=1,nup*nup
   45       cslmui(l,idet)=cslmin(l,idet)
          do 50 j=1,nup
            do 50 idim=1,ndim
   50         cfpu(idim,j+ikel,idet)=cdorb(idim,iworbd(j,idet))
          do 60 i=1,nup
            do 60 idim=1,ndim
   60         cddeti_deti(idim,i,idet)=cddeti_detin(idim,i,idet)
       else
        ikel=ndn*(iel-nup-1)
        do 80 idet=1,ndetdn
          cdetd(idet)=cdetn(idet)
          do 65 j=1,ndn*ndn
   65       cslmdi(j,idet)=cslmin(j,idet)
          do 70 j=1,ndn
            do 70 idim=1,ndim
   70         cfpd(idim,j+ikel,idet)=cdorb(idim,iworbd(j+nup,idet))
          do 80 i=1+nup,ndn+nup
            do 80 idim=1,ndim
   80         cddeti_deti(idim,i,idet)=cddeti_detin(idim,i,idet)
      endif

      return
      end

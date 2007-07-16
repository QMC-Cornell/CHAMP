      subroutine detsav(iel)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
      common /slater/ slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD)
     &,fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
     &,fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
     &,detu(MDETUD),detd(MDETUD)
     &,ddeti_deti(3,MELEC,MDET),d2edeti_deti(MELEC,MDET),deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD),d2deti_det(MPARMD),d2det_det
     &,detij_det(MPARMD,MPARMD)
      common /slatn/ slmin(MMAT_DIM,MDETUD),detn(MDETUD)
     &,ddeti_detin(3,MELEC,MDETUD),d2edeti_detin(MELEC,MDETUD)
     &,dorb(3,MORB),ddorb(MORB)

      if(iel.le.nup) then
        ikel=nup*(iel-1)
        do 60 idet=1,ndetup
          detu(idet)=detn(idet)
          do 45 l=1,nup*nup
   45       slmui(l,idet)=slmin(l,idet)
          do 50 j=1,nup
            fpu(1,j+ikel,idet)=dorb(1,iworbd(j,idet))
            fpu(2,j+ikel,idet)=dorb(2,iworbd(j,idet))
   50       fpu(3,j+ikel,idet)=dorb(3,iworbd(j,idet))
          do 60 i=1,nup
            ddeti_deti(1,i,idet)=ddeti_detin(1,i,idet)
            ddeti_deti(2,i,idet)=ddeti_detin(2,i,idet)
   60       ddeti_deti(3,i,idet)=ddeti_detin(3,i,idet)

       else

        ikel=ndn*(iel-nup-1)
        do 80 idet=1,ndetdn
          detd(idet)=detn(idet)
          do 65 j=1,ndn*ndn
   65       slmdi(j,idet)=slmin(j,idet)
          do 70 j=1,ndn
            fpd(1,j+ikel,idet)=dorb(1,iworbd(j+nup,idet))
            fpd(2,j+ikel,idet)=dorb(2,iworbd(j+nup,idet))
   70       fpd(3,j+ikel,idet)=dorb(3,iworbd(j+nup,idet))
          do 80 i=1+nup,ndn+nup
            ddeti_deti(1,i,idet)=ddeti_detin(1,i,idet)
            ddeti_deti(2,i,idet)=ddeti_detin(2,i,idet)
   80       ddeti_deti(3,i,idet)=ddeti_detin(3,i,idet)
      endif

      return
      end

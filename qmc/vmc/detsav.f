      subroutine detsav(iel)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use dorb_mod
      use slatn_mod
      use orbe_mod
      use dets_mod
      use slater_mod
      implicit real*8(a-h,o-z)


      if(iel.le.nup) then
        ikel=nup*(iel-1)
        do 60 idet=1,ndetup
          detu(idet)=detn(idet)
          do 45 l=1,nup*nup
   45       slmui(l,idet)=slmin(l,idet)
          do 50 j=1,nup
            fpu(1,j+ikel,idet)=dorbe(1,iworbd(j,idet))
            fpu(2,j+ikel,idet)=dorbe(2,iworbd(j,idet))
   50       fpu(3,j+ikel,idet)=dorbe(3,iworbd(j,idet))
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
            fpd(1,j+ikel,idet)=dorbe(1,iworbd(j+nup,idet))
            fpd(2,j+ikel,idet)=dorbe(2,iworbd(j+nup,idet))
   70       fpd(3,j+ikel,idet)=dorbe(3,iworbd(j+nup,idet))
          do 80 i=1+nup,ndn+nup
            ddeti_deti(1,i,idet)=ddeti_detin(1,i,idet)
            ddeti_deti(2,i,idet)=ddeti_detin(2,i,idet)
   80       ddeti_deti(3,i,idet)=ddeti_detin(3,i,idet)
      endif

      return
      end

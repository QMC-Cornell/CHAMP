      subroutine cuspcheck3(diff,iprin)
c Written by Claudia Filippi

      use constants_mod
      use atom_mod
      use contr2_mod
      use jaspar3_mod
      implicit real*8(a-h,o-z)

      character*1 ipw,ipw1,ipw2

      dimension diff(*),o(2*nord),op(0:2*nord)

      isht=0
      do 200 it=1,nctype

      do 20 jp=1,nord
        o(jp)=zero
        op(jp-1)=zero
        o(jp+nord)=zero
  20    op(jp+nord-1)=zero

      ll=0
      do 30 jp=1,nord
        jpm=jp-1
        jpmsh=jpm+nord
        do 30 ju=jp,0,-1
          jsx=jp-ju
          do 30 js=jsx,0,-1
            ll=ll+1
            jt=jsx-js
            jsg=1
            if(mod(jt,2).ne.0) jsg=-1
            o(jp)=o(jp)+jsg*c(ll,it,1)
            if(js.ge.1) op(jpm)=op(jpm)+jsg*js*c(ll,it,1)
            if(jt.ge.1) op(jpm)=op(jpm)-jsg*jt*c(ll,it,1)

            if(jt.eq.0) then
              if(ju.eq.0) o(jp+nord)=o(jp+nord)+c(ll,it,1)
              if(ju.eq.1) op(jpmsh)=op(jpmsh)+c(ll,it,1)
            endif
  30        continue

      nfsh=0
      if(ifock.gt.0) call fockcheck3(o,op,diff,isht,nfsh,it,iprin)

      ishe=nord+nfsh
      do 40 jp=1,nord
        diff(isht+jp)=op(jp-1)
  40    diff(isht+jp+ishe)=op(jp+nord-1)

      if(iprin.ge.1) then
        write(6,'(''atom type'',i3)') it
        do 60 jp=1,nord
          write(ipw,'(i1)') jp-1
  60      write(6,'(''coefs of r**'',a1,9f12.6)') ipw,diff(isht+jp)
        do 70 jp=1,nord
          write(ipw,'(i1)') jp-1
  70      write(6,'(''coefs of s**'',a1,9f12.6)') ipw,diff(isht+jp+ishe)
        do 80 jp=1,nord
          write(ipw1,'(i1)') jp-1
          write(ipw2,'(i1)') jp
          write(6,'(''e-n op('',a1,'')'',f12.6,4x,''o('',a1,'')'',f12.6)') ipw1,op(jp-1),ipw2,o(jp)
  80      write(6,'(''e-e op('',a1,'')'',f12.6,4x,''o('',a1,'')'',f12.6)') ipw1,op(jp+nord-1),ipw2,o(jp+nord)
      endif

 200  isht=isht+2*ishe
      return
      end

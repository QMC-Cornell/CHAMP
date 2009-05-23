      subroutine hpsiedmc(iel,iw,coord,psid,psij,force)
c Written by Claudia Filippi

      use const_mod
      use dim_mod
      use forcepar_mod
      use force_dmc_mod
      use forcest_dmc_mod
      implicit real*8(a-h,o-z)
!JT      include '../../vmc/vmc.h'
!JT      include '../dmc.h'
!JT      include '../../vmc/force.h'
!JT      common /dim/ ndim
!JT      common /forcepar/ deltot(MFORCE),nforce,istrech
!JT      common /forcest_dmc/ fgcum(MFORCE),fgcm2(MFORCE)
!JT      common /force_dmc/ itausec,nwprod

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /config_dmc/ xoldw(3,MELEC,MWALK,MFORCE),voldw(3,MELEC,MWALK,MFORCE),
     &psidow(MWALK,MFORCE),psijow(MWALK,MFORCE),peow(MWALK,MFORCE),peiow(MWALK,MFORCE),d2ow(MWALK,MFORCE)

      dimension coord(3),force(3,*),x(3,MELEC)

      do 10 k=1,ndim
      do 10 i=1,iel-1
  10    x(k,i)=xoldw(k,i,iw,1)

      do 20 k=1,ndim
  20    x(k,iel)=coord(k)

      do 30 k=1,ndim
      do 30 i=iel+1,nelec
  30    x(k,i)=xoldw(k,i,iw,1)

      call hpsie(iel,x,psid,psij,force)

      return
      end

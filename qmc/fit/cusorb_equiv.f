      subroutine cusorb_equiv(icent,orb2)
c Written by Cyrus Umrigar
c Calculate that part, orb2, of orbitals, orb2, at position of nucleus
c icent coming from basis, ibas, and symmetry-related bases.

      use coefs_mod
      implicit real*8(a-h,o-z)
!JT      include '../vmc/vmc.h'
!JT      include 'fit.h'
!JT      include '../vmc/force.h'

      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC)
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /coefs2/ coef2(MBASIS,MORB,MCENT)
      dimension orb2(*)

c Note that there is no need to recalculate basis fns., phin, since
c they were already calculated in cusorb.

c Set coef2=coef for basis fns. related by symmetry to basis ibas, and
c to zero otherwise.
      call equiv_bas

c Calculate that part, orb2, of orbitals, orb2, at position of nucleus
c icent coming from basis, ibas, and symmetry-related bases.
      do 10 iorb=1,norb
        orb2(iorb)=0
        do 10 m=1,nbasis
   10     orb2(iorb)=orb2(iorb)+coef2(m,iorb,icent)*phin(m,1)

      return
      end

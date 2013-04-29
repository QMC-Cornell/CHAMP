      subroutine cusorb_equiv(icent,orb2)
! Written by Cyrus Umrigar
! Calculate that part, orb2, of orbitals, orb2, at position of nucleus
! icent coming from basis, ibas, and symmetry-related bases.
      use coefs_mod
      use coefs2_mod
      use phifun_mod
      implicit real*8(a-h,o-z)

      dimension orb2(*)

! Note that there is no need to recalculate basis fns., phin, since
! they were already calculated in cusorb.

! Set coef2=coef for basis fns. related by symmetry to basis ibas, and
! to zero otherwise.
      call equiv_bas

! Calculate that part, orb2, of orbitals, orb2, at position of nucleus
! icent coming from basis, ibas, and symmetry-related bases.
      do 10 iorb=1,norb
        orb2(iorb)=0
        do 10 m=1,nbasis
   10     orb2(iorb)=orb2(iorb)+coef2(m,iorb,icent)*phin(m,1)

      return
      end

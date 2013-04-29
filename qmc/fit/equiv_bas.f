      subroutine equiv_bas
! Written by Cyrus Umrigar
! Set coef2=coef for basis fns. related by symmetry to basis ibas, and
! to zero otherwise.

! J. Toulouse - 08 Jan 05: change coef(i,j,1) -> coef(i,j,iwf)
      use all_tools_mod
      use atom_mod
      use coefs_mod
      use coefs2_mod
      use optim_mod
      use wfsec_mod
      implicit real*8(a-h,o-z)

      call alloc ('coef2', coef2, nbasis, norb, ncent)

      do 10 icent=1,ncent
        do 10 iorb=1,norb
          do 10 ib=1,nbasis
   10       coef2(ib,iorb,icent)=0

      do 20 icent=1,ncent
        ibas=imnbas(icent)
          do 20 i=1,necn
             if(iebasi(1,i).eq.ibas .or. iebasi(2,i).eq.ibas) then
               if(coef(iebasi(1,i),iabs(ieorb(1,i)),iwf).eq.0)
     &            coef(iebasi(1,i),iabs(ieorb(1,i)),iwf)=1.d-100
               if(coef(iebasi(2,i),iabs(ieorb(2,i)),iwf).eq.0)
     &            coef(iebasi(2,i),iabs(ieorb(2,i)),iwf)=1.d-100
               coef2(iebasi(1,i),iabs(ieorb(1,i)),icent)=
     &         coef(iebasi(1,i),iabs(ieorb(1,i)),iwf)
               coef2(iebasi(2,i),iabs(ieorb(2,i)),icent)=
     &         coef(iebasi(2,i),iabs(ieorb(2,i)),iwf)
             endif
   20     continue

!     do 30 icent=1,ncent
!       do 30 iorb=1,norb
!  30     write(6,'(7f9.3)') (coef2(ibas,iorb,icent),ibas=1,7)

      return
      end

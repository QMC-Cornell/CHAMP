      subroutine equiv_bas
c Written by Cyrus Umrigar
c Set coef2=coef for basis fns. related by symmetry to basis ibas, and
c to zero otherwise.

! J. Toulouse - 08 Jan 05: change coef(i,j,1) -> coef(i,j,iwf)
!
      use atom_mod
      use coefs_mod
      use optim_mod
      implicit real*8(a-h,o-z)
!JT      include '../vmc/vmc.h'
!JT      include 'fit.h'
!JT      include '../vmc/force.h'

!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /coefs2/ coef2(MBASIS,MORB,MCENT)
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase

c     common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /wfsec/ iwftype(MFORCE),iwf,nwftype !JT

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

c     do 30 icent=1,ncent
c       do 30 iorb=1,norb
c  30     write(6,'(7f9.3)') (coef2(ibas,iorb,icent),ibas=1,7)

      return
      end

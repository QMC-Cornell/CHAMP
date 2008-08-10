      integer lo(MORB),npoint(MORB) 
      integer iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE)
      integer iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM)
      integer iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM)
      integer imnbas(MCENT)
      integer nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj
      integer nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE)
      integer necn,nebase

      common /optim/ lo,npoint,iwjasa,iwjasb,iwjasc,iwjasf,iwbase,iwbasi,iworb,iwcsf,iebase,iebasi,ieorb,imnbas,nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,nparma,nparmb,nparmc,nparmf,necn,nebase

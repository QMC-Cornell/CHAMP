      subroutine walksav_det(iw)
! Written by Claudia Filippi, modified by Cyrus Umrigar

      use dorb_mod
      use dets_mod
      use slater_mod
      use slaterw_mod
      use const_mod
      implicit real*8(a-h,o-z)

      do 45 idet=1,ndetup
        detuw(idet,iw)=detu(idet)
        do 40 j=1,nup
          ddeti_detiw(1,j,idet,iw)=ddeti_deti(1,j,idet)
          ddeti_detiw(2,j,idet,iw)=ddeti_deti(2,j,idet)
   40     ddeti_detiw(3,j,idet,iw)=ddeti_deti(3,j,idet)
        do 45 j=1,nup*nup
          slmuiw(j,idet,iw)=slmui(j,idet)
          fpuw(1,j,idet,iw)=fpu(1,j,idet)
          fpuw(2,j,idet,iw)=fpu(2,j,idet)
   45     fpuw(3,j,idet,iw)=fpu(3,j,idet)

      do 55 idet=1,ndetdn
        detdw(idet,iw)=detd(idet)
        do 50 j=nup+1,nelec
          ddeti_detiw(1,j,idet,iw)=ddeti_deti(1,j,idet)
          ddeti_detiw(2,j,idet,iw)=ddeti_deti(2,j,idet)
   50     ddeti_detiw(3,j,idet,iw)=ddeti_deti(3,j,idet)
        do 55 j=1,ndn*ndn
          slmdiw(j,idet,iw)=slmdi(j,idet)
          fpdw(1,j,idet,iw)=fpd(1,j,idet)
          fpdw(2,j,idet,iw)=fpd(2,j,idet)
   55     fpdw(3,j,idet,iw)=fpd(3,j,idet)

      return
      entry walkstrdet(iw)

      do 85 idet=1,ndetup
        detu(idet)=detuw(idet,iw)
        do 80 j=1,nup
          ddeti_deti(1,j,idet)=ddeti_detiw(1,j,idet,iw)
          ddeti_deti(2,j,idet)=ddeti_detiw(2,j,idet,iw)
   80     ddeti_deti(3,j,idet)=ddeti_detiw(3,j,idet,iw)
        do 85 j=1,nup*nup
          slmui(j,idet)=slmuiw(j,idet,iw)
          fpu(1,j,idet)=fpuw(1,j,idet,iw)
          fpu(2,j,idet)=fpuw(2,j,idet,iw)
   85     fpu(3,j,idet)=fpuw(3,j,idet,iw)

      do 95 idet=1,ndetdn
        detd(idet)=detdw(idet,iw)
        do 90 j=nup+1,nelec
          ddeti_deti(1,j,idet)=ddeti_detiw(1,j,idet,iw)
          ddeti_deti(2,j,idet)=ddeti_detiw(2,j,idet,iw)
   90     ddeti_deti(3,j,idet)=ddeti_detiw(3,j,idet,iw)
        do 95 j=1,ndn*ndn
          slmdi(j,idet)=slmdiw(j,idet,iw)
          fpd(1,j,idet)=fpdw(1,j,idet,iw)
          fpd(2,j,idet)=fpdw(2,j,idet,iw)
   95     fpd(3,j,idet)=fpdw(3,j,idet,iw)

      return
      entry splitjdet(iw,iw2)

      do 125 idet=1,ndetup
        detuw(idet,iw2)=detuw(idet,iw)
        do 120 j=1,nup
          ddeti_detiw(1,j,idet,iw2)=ddeti_detiw(1,j,idet,iw)
          ddeti_detiw(2,j,idet,iw2)=ddeti_detiw(2,j,idet,iw)
  120     ddeti_detiw(3,j,idet,iw2)=ddeti_detiw(3,j,idet,iw)
        do 125 j=1,nup*nup
          slmuiw(j,idet,iw2)=slmuiw(j,idet,iw)
          fpuw(1,j,idet,iw2)=fpuw(1,j,idet,iw)
          fpuw(2,j,idet,iw2)=fpuw(2,j,idet,iw)
  125     fpuw(3,j,idet,iw2)=fpuw(3,j,idet,iw)

      do 135 idet=1,ndetdn
        detdw(idet,iw2)=detdw(idet,iw)
        do 130 j=nup+1,nelec
          ddeti_detiw(1,j,idet,iw2)=ddeti_detiw(1,j,idet,iw)
          ddeti_detiw(2,j,idet,iw2)=ddeti_detiw(2,j,idet,iw)
  130     ddeti_detiw(3,j,idet,iw2)=ddeti_detiw(3,j,idet,iw)
        do 135 j=1,ndn*ndn
          slmdiw(j,idet,iw2)=slmdiw(j,idet,iw)
          fpdw(1,j,idet,iw2)=fpdw(1,j,idet,iw)
          fpdw(2,j,idet,iw2)=fpdw(2,j,idet,iw)
  135     fpdw(3,j,idet,iw2)=fpdw(3,j,idet,iw)

      return
      end

      subroutine walksav_det(iw)
! Written by Claudia Filippi, modified by Cyrus Umrigar

      use dorb_mod
      use dets_mod
      use slater_mod
      use slaterw_mod
      use const_mod
      use deriv_fast_mod
      use gammaw_mod
      use qua_mod, only: quadx, quadr
      implicit real*8(a-h,o-z)

      !TA
      orbw(:,:,iw)    = orb
      dorbw(:,:,:,iw) = dorb
 
      aiupw(:,:,iw) = aiup
      deta_upw(iw)  = deta_up
      tupw(:,:,iw)  = tup
      detupw(:,iw)  = detup
      do idet=1,ndetup
        iord = ordup(idet)
        do i=1,iord
          do j=1,iord
            invupw(j+iord*(i-1),idet,iw) = invup(j+iord*(i-1),idet)
          enddo
        enddo
      enddo
 
      aidnw(:,:,iw) = aidn
      deta_dnw(iw)  = deta_dn
      tdnw(:,:,iw)  = tdn
      detdnw(:,iw)  = detdn
      do idet=1,ndetdn
        iord = orddn(idet)
        do i=1,iord
          do j=1,iord
            invdnw(j+iord*(i-1),idet,iw) = invdn(j+iord*(i-1),idet)
          enddo
        enddo
      enddo

      chiw(iw)     = chi
      yupw(:,:,iw) = yup
      ydnw(:,:,iw) = ydn

!      do 45 idet=1,ndetup
!        detuw(idet,iw)=detu(idet)
!        do 40 j=1,nup
!          ddeti_detiw(1,j,idet,iw)=ddeti_deti(1,j,idet)
!          ddeti_detiw(2,j,idet,iw)=ddeti_deti(2,j,idet)
!   40     ddeti_detiw(3,j,idet,iw)=ddeti_deti(3,j,idet)
!        do 45 j=1,nup*nup
!          slmuiw(j,idet,iw)=slmui(j,idet)
!          fpuw(1,j,idet,iw)=fpu(1,j,idet)
!          fpuw(2,j,idet,iw)=fpu(2,j,idet)
!   45     fpuw(3,j,idet,iw)=fpu(3,j,idet)
!
!      do 55 idet=1,ndetdn
!        detdw(idet,iw)=detd(idet)
!        do 50 j=nup+1,nelec
!          ddeti_detiw(1,j,idet,iw)=ddeti_deti(1,j,idet)
!          ddeti_detiw(2,j,idet,iw)=ddeti_deti(2,j,idet)
!   50     ddeti_detiw(3,j,idet,iw)=ddeti_deti(3,j,idet)
!        do 55 j=1,ndn*ndn
!          slmdiw(j,idet,iw)=slmdi(j,idet)
!          fpdw(1,j,idet,iw)=fpd(1,j,idet)
!          fpdw(2,j,idet,iw)=fpd(2,j,idet)
!   55     fpdw(3,j,idet,iw)=fpd(3,j,idet)

      return
      entry walkstrdet(iw)

      !TA
      orb  = orbw(:,:,iw)
      dorb = dorbw(:,:,:,iw)
 
      aiup     = aiupw(:,:,iw)
      deta_up  = deta_upw(iw)
      tup      = tupw(:,:,iw)
      detup    = detupw(:,iw)
      do idet=1,ndetup
        iord = ordup(idet)
        do i=1,iord
          do j=1,iord
            invup(j+iord*(i-1),idet) = invupw(j+iord*(i-1),idet,iw)
          enddo
        enddo
      enddo
 
      aidn     = aidnw(:,:,iw)
      deta_dn  = deta_dnw(iw)
      tdn      = tdnw(:,:,iw)
      detdn    = detdnw(:,iw)
      do idet=1,ndetdn
        iord = orddn(idet)
        do i=1,iord
          do j=1,iord
            invdn(j+iord*(i-1),idet) = invdnw(j+iord*(i-1),idet,iw)
          enddo
        enddo
      enddo
 
      chi = chiw(iw)
      yup = yupw(:,:,iw)
      ydn = ydnw(:,:,iw)

!      do 85 idet=1,ndetup
!        detu(idet)=detuw(idet,iw)
!        do 80 j=1,nup
!          ddeti_deti(1,j,idet)=ddeti_detiw(1,j,idet,iw)
!          ddeti_deti(2,j,idet)=ddeti_detiw(2,j,idet,iw)
!   80     ddeti_deti(3,j,idet)=ddeti_detiw(3,j,idet,iw)
!        do 85 j=1,nup*nup
!          slmui(j,idet)=slmuiw(j,idet,iw)
!          fpu(1,j,idet)=fpuw(1,j,idet,iw)
!          fpu(2,j,idet)=fpuw(2,j,idet,iw)
!   85     fpu(3,j,idet)=fpuw(3,j,idet,iw)
!
!      do 95 idet=1,ndetdn
!        detd(idet)=detdw(idet,iw)
!        do 90 j=nup+1,nelec
!          ddeti_deti(1,j,idet)=ddeti_detiw(1,j,idet,iw)
!          ddeti_deti(2,j,idet)=ddeti_detiw(2,j,idet,iw)
!   90     ddeti_deti(3,j,idet)=ddeti_detiw(3,j,idet,iw)
!        do 95 j=1,ndn*ndn
!          slmdi(j,idet)=slmdiw(j,idet,iw)
!          fpd(1,j,idet)=fpdw(1,j,idet,iw)
!          fpd(2,j,idet)=fpdw(2,j,idet,iw)
!   95     fpd(3,j,idet)=fpdw(3,j,idet,iw)

      return
      entry splitjdet(iw,iw2)

      !TA
      orbw(:,:,iw2)    = orbw(:,:,iw)
      dorbw(:,:,:,iw2) = dorbw(:,:,:,iw)
 
      aiupw(:,:,iw2) = aiupw(:,:,iw)
      deta_upw(iw2)  = deta_upw(iw)
      tupw(:,:,iw2)  = tupw(:,:,iw)
      detupw(:,iw2)  = detupw(:,iw)
      do idet=1,ndetup
        iord = ordup(idet)
        do i=1,iord
          do j=1,iord
            invupw(j+iord*(i-1),idet,iw2) = invupw(j+iord*(i-1),idet,iw)
          enddo
        enddo
      enddo
 
      aidnw(:,:,iw2) = aidnw(:,:,iw)
      deta_dnw(iw2)  = deta_dnw(iw)
      tdnw(:,:,iw2)  = tdnw(:,:,iw)
      detdnw(:,iw2)  = detdnw(:,iw)
      do idet=1,ndetdn
        iord = orddn(idet)
        do i=1,iord
          do j=1,iord
            invdnw(j+iord*(i-1),idet,iw2) = invdnw(j+iord*(i-1),idet,iw)
          enddo
        enddo
      enddo
 
      chiw(iw2)     = chiw(iw)
      yupw(:,:,iw2) = yupw(:,:,iw)
      ydnw(:,:,iw2) = ydnw(:,:,iw)
      quadr(:,:,iw2)=quadr(:,:,iw)
      quadx(:,:,:,iw2)=quadx(:,:,:,iw)

!      do 125 idet=1,ndetup
!        detuw(idet,iw2)=detuw(idet,iw)
!        do 120 j=1,nup
!          ddeti_detiw(1,j,idet,iw2)=ddeti_detiw(1,j,idet,iw)
!          ddeti_detiw(2,j,idet,iw2)=ddeti_detiw(2,j,idet,iw)
!  120     ddeti_detiw(3,j,idet,iw2)=ddeti_detiw(3,j,idet,iw)
!        do 125 j=1,nup*nup
!          slmuiw(j,idet,iw2)=slmuiw(j,idet,iw)
!          fpuw(1,j,idet,iw2)=fpuw(1,j,idet,iw)
!          fpuw(2,j,idet,iw2)=fpuw(2,j,idet,iw)
!  125     fpuw(3,j,idet,iw2)=fpuw(3,j,idet,iw)
!
!      do 135 idet=1,ndetdn
!        detdw(idet,iw2)=detdw(idet,iw)
!        do 130 j=nup+1,nelec
!          ddeti_detiw(1,j,idet,iw2)=ddeti_detiw(1,j,idet,iw)
!          ddeti_detiw(2,j,idet,iw2)=ddeti_detiw(2,j,idet,iw)
!  130     ddeti_detiw(3,j,idet,iw2)=ddeti_detiw(3,j,idet,iw)
!        do 135 j=1,ndn*ndn
!          slmdiw(j,idet,iw2)=slmdiw(j,idet,iw)
!          fpdw(1,j,idet,iw2)=fpdw(1,j,idet,iw)
!          fpdw(2,j,idet,iw2)=fpdw(2,j,idet,iw)
!  135     fpdw(3,j,idet,iw2)=fpdw(3,j,idet,iw)

      return
      end

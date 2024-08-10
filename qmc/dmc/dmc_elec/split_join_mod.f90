module split_join_mod
contains
      subroutine split_join(walkers, nwalk)
! Originally written by Cyrus Umrigar, minor changes by Tyler Anderson

      use all_tools_mod
      use constants_mod
      use control_mod
      use const_mod
      use dim_mod
      use forcepar_mod
      use contrl_per_mod
      use force_dmc_mod
      use jacobsave_mod
      use config_dmc_mod
      use branch_mod, only: MWALK, pwt, wthist
      use stats_mod
      use age_mod
      use velratio_mod
!      use walker_type_mod, only: walker_t
      use walker_type_mod
      implicit real*8(a-h,o-z)

      type(walker_t), target, intent(inout) :: walkers(:)
      integer, intent(inout)        :: nwalk
      integer :: ipr_sav=0
      save ipr_sav

      dimension iwundr(MWALK),wt_sav(MWALK)

!     write(6,'(''nwalk at beg of splitjoin'',9i5)') nwalk

      if(ipr.ge.2) then
        write(6,'(''nwalk,,ioldest,ioldestmx='',9i5)') nwalk,ioldest,ioldestmx
        write(6,'(''wt='',100f6.2)') (walkers(iww)%weight,iww=1,nwalk)
!        write(6,'(''eold='',100f6.2)') (eoldw(iww,1),iww=1,nwalk)
      endif

! Identify the walkers that are to be killed because their wt. is zero, or that are
! to be merged because their wts are less than 1/2.  For each pair whose wts are < 1/2
! give the one that is to be kept the additional wt of the other and put the other one on a stack.
      iunder=0
      ipair=0
      wtsm=zero
      do 10 iw=1,nwalk
        wtsm=wtsm+walkers(iw)%weight
        if(walkers(iw)%weight.lt.half) then
          if(walkers(iw)%weight.eq.zero) then
            nbrnch=nbrnch+1
            iunder=iunder+1
            iwundr(iunder)=iw
!           write(6,'(''Eliminating walker'',i5)') iw
           else
            if(ipair.eq.0) then
              ipair=1
              iw2=iw
             else
              nbrnch=nbrnch+1
              ipair=0
              iunder=iunder+1
              wttot=walkers(iw)%weight+walkers(iw2)%weight
              if(rannyu(0).gt.(walkers(iw)%weight/wttot)) then
                walkers(iw2)%weight=wttot
                iwundr(iunder)=iw
               else
                walkers(iw)%weight=wttot
                iwundr(iunder)=iw2
              endif
!             write(6,'(''Joined walkers'',9i5)') iwundr(iunder),iw2,iw
            endif
          endif
        endif
   10 continue

! Figure out what weight walkers can be split without exceeding MWALK
! Warning: MWALK can still be exceeded on very rare occasions if a walker gets split more than once.
      nwalk_after_join=nwalk-iunder
      nwalk_max_dupl=MWALK-nwalk_after_join
      do iw=1,nwalk
        wt_sav(iw)=walkers(iw)%weight
      enddo
      call shell(wt_sav,nwalk)
      wt_split=max(2.d0,wt_sav(max(1,nwalk-nwalk_max_dupl)))
      if(iunder.gt.0) then
!       write(6,'(''Number of walkers joined,nwalk,nwalk_after_join'',9i5)') iunder,nwalk,nwalk_after_join
!       write(6,'(''Holes are at walkers'',999i5)') (iwundr(i),i=1,iunder)
      endif
!     write(6,'(''wt_split1='',f8.2)') wt_split

! Figure out the number of walkers that could get split more than once.
! Whether it actually gets split again or not depends on where the duplicated walker gets placed in the list of walkers.
! When there is a split, one walker stays where it was and the other fills an empty slot.  It is only the one in the
! empty slot that will get split again (it will split again only if it is put in a higher indexed position than the current walker),
! so the number to add is 1,2,3,... rather than 1,3,7,...
      n_split_more_than_once=0
      do iw=max(1,nwalk-nwalk_max_dupl),nwalk
        if(wt_sav(iw)/wt_split.ge.8.d0) then
          n_split_more_than_once=n_split_more_than_once+3
         elseif(wt_sav(iw)/wt_split.ge.4.d0) then
          n_split_more_than_once=n_split_more_than_once+2
         elseif(wt_sav(iw)/wt_split.ge.2.d0) then
          n_split_more_than_once=n_split_more_than_once+1
        endif
      enddo

      if(n_split_more_than_once.gt.0) then
        ipr_sav=ipr_sav+1
        if(ipr_sav.le.50) then
          write(6,'(''Warning: Number of walkers on master process split more than once (not of concern)='',i4)') n_split_more_than_once
        elseif(ipr_sav.eq.51) then
          write(6,'(''Warning: Additional warning msgs. of Number of walkers on master process split more than once suppressed'')')
        endif
      endif

! Adjust what weight walkers can be split without exceeding MWALK
! If nwalk_max_dupl<n_split_more_than_once then set wt_split to huge number.
! This is not the best choice but it is a conservative choice.
      if(nwalk-nwalk_max_dupl+n_split_more_than_once.lt.nwalk) then
        wt_split=max(2.d0,wt_sav(max(1,nwalk-nwalk_max_dupl+n_split_more_than_once)))
       else
        wt_split=1.d99
      endif
!     write(6,'(''wt_split2='',f8.2)') wt_split

! Split the walkers whose wt is >max(two,wt_split).  If there are walkers that were eliminated, so that iunder>0
! then put the new walker in that location.  Otherwise put it at the end.
      nwalk2=nwalk
      do 20 iw=1,nwalk
        if(walkers(iw)%weight.gt.wt_split) then
          nbrnch=nbrnch+1
          if(iunder.gt.0) then
            iw2=iwundr(iunder)
            iunder=iunder-1
           else
            nwalk2=nwalk2+1
            iw2=nwalk2
            if(nwalk2.gt.MWALK) then
              write(6,'(''iw,nwalk,nwalk2,MWALK='',9i5)') iw,nwalk,nwalk2,MWALK
              write(6,'(''wt='',100f6.2)') (walkers(iww)%weight,iww=1,nwalk)
              stop 'MWALK exceeded in splitj_mov1'
            endif
          endif
          walkers(iw)%weight=walkers(iw)%weight*half
          walkers(iw2)%weight=walkers(iw)%weight
          walkers(iw2)%x=walkers(iw)%x
          walkers(iw2)%psi=walkers(iw)%psi
          iage(iw2)=iage(iw)
          do ifr=1,nforce
            fratio(iw2,ifr)=fratio(iw,ifr)
!            ajacold(iw2,ifr)=ajacold(iw,ifr)
!            eoldw(iw2,ifr)=eoldw(iw,ifr)
!            psidow(iw2,ifr)=psidow(iw,ifr)
!            psijow(iw2,ifr)=psijow(iw,ifr)
!            peow(iw2,ifr)=peow(iw,ifr)
!            peiow(iw2,ifr)=peiow(iw,ifr)
!            d2ow(iw2,ifr)=d2ow(iw,ifr)
            pwt(iw2,ifr)=pwt(iw,ifr)
            do ip=0,nwprod-1
              wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
            enddo
          enddo
!            v2=0
!            do 15 i=1,nelec
!              do 15 k=1,ndim
!                voldw(k,i,iw2,ifr)=voldw(k,i,iw,ifr)
!                v2=v2+voldw(k,i,iw2,ifr)**2
!   15           xoldw(k,i,iw2,ifr)=xoldw(k,i,iw,ifr)
!         if(iage(iw2).ge.1) write(6,'(''iage,wt='',i5,f8.2)') iage(iw2),walkers(iw2)%weight
!         write(6,'(''iage,iw2,wt,e,v2='',2i5,9f9.2)') iage(iw2),iw2,walkers(iw2)%weight,eoldw(iw2,1),v2
        endif
   20 continue

! If more walkers were eliminated than the number duplicated then consolidate
! the remaining walkers so that they are they occupy the first positions.
      do j=iunder,1,-1
        iw2=iwundr(j)
        iw=nwalk2
!       write(6,'(''Moving walker'',i5,'' to walker'',i5)') iw,iw2
        nwalk2=nwalk2-1
!        walkers(iw2)%weight=walkers(iw)%weight
!        walkers(iw2)%x=walkers(iw)%x
!        walkers(iw2)%psi=walkers(iw)%psi
        walkers(iw2)=walkers(iw)
        iage(iw2)=iage(iw)
        do ifr=1,nforce
          fratio(iw2,ifr)=fratio(iw,ifr)
!          ajacold(iw2,ifr)=ajacold(iw,ifr)
!          eoldw(iw2,ifr)=eoldw(iw,ifr)
!          psidow(iw2,ifr)=psidow(iw,ifr)
!          psijow(iw2,ifr)=psijow(iw,ifr)
!          peow(iw2,ifr)=peow(iw,ifr)
!          peiow(iw2,ifr)=peiow(iw,ifr)
!          d2ow(iw2,ifr)=d2ow(iw,ifr)
          pwt(iw2,ifr)=pwt(iw,ifr)
          do ip=0,nwprod-1
            wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
          enddo
        enddo
!          v2=0
!          do 15 i=1,nelec
!            do 15 k=1,ndim
!              voldw(k,i,iw2,ifr)=voldw(k,i,iw,ifr)
!              v2=v2+voldw(k,i,iw2,ifr)**2
!   15         xoldw(k,i,iw2,ifr)=xoldw(k,i,iw,ifr)
      enddo
      nwalk=nwalk2
      call object_modified_by_index (nwalk_index)
!     write(6,'(''nwalk at end of splitjoin'',9i5)') nwalk

      wtsm2=zero
      do 40 iw=1,nwalk
        wtsm2=wtsm2+walkers(iw)%weight
!       if(walkers(iw)%weight.lt.half) write(11,'(i4,9d12.5)') iw,walkers(iw)%weight,eoldw(iw)
!       if(walkers(iw)%weight.gt.two) write(11,'(i4,9d12.5)') iw,walkers(iw)%weight,eoldw(iw)
   40 continue
!     if(dabs(wtsm-wtsm2).gt.1.d-10) write(11,'(2f12.6)') wtsm,wtsm2

      return
      end subroutine split_join
end module split_join_mod

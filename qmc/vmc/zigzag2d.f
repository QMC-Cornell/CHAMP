      subroutine zigzag2d(p,q,xold,xnew,ielec)

c Written by Abhijit Mehta, December 2011
c  Calculates quantities useful for studying zigzag quantum phase
c  transition in rings and wires.  
c  -reduced pair density
c  -"staggered amplitude"
c  xold and xnew are the old and new configurations
c  p and q are the probabilities of accept and reject
c  ielec labels the electron that was moved for 1-electron moves
c  if ielec=0, we are doing an all-electron move.

      use dets_mod
      use const_mod
      use dim_mod
      use contrl_per_mod, only: iperiodic
      use pairden_mod
      use periodic_1d_mod, only: alattice
      use zigzag_mod
      implicit real*8(a-h,o-z)
      logical l_oldneoldsav, l_oldnenewsav
c     common /circularmesh/ delti
      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys
      dimension xold(3,nelec),xnew(3,nelec)
      dimension temppos(2),zzterm(nzzvars)
      dimension zzmaglocal_new(nelec),zzmaglocal_old(nelec)
      dimension zzcorrmat_old(nelec,nelec),zzcorrmat_new(nelec,nelec)

      if(ielec.lt.0 .or. ielec.gt.nelec) then
        write (6,*) 'Bad value of ielec in zigzag2d, ', ielec
        stop "Bad value of ielec in zigzag2d"
      endif

c   zzpos will ultimately hold the sorted electron positions.  The first index labels longitudinal
c      or transverse coordinate (i.e., zzpos(1,nelec) is sorted in ascending order)

c   First, save some sorting work, and see if one of the configurations xnew or xold was already sorted
c    xold_sav and xnew_sav contain the last values of xold and xnew used by this subroutine
c    iold_indices and inew_indices contain the mapping of how the electrons in xold, xnew were sorted
     
      l_oldneoldsav = .false.
      l_oldnenewsav = .false.
      outerloop: do i1 = 1,nelec
        do i2 = 1,2
          if(xold(i2,i1).ne.xold_sav(i2,i1)) then
            l_oldneoldsav = .true.
          endif
          if(xold(i2,i1).ne.xnew_sav(i2,i1)) then
            l_oldnenewsav = .true.
          endif
          if(l_oldneoldsav.and.l_oldnenewsav) exit outerloop 
        enddo
      end do outerloop
      ! if(xold.ne.xold_sav)
      if(l_oldneoldsav) then ! xold has changed
        ! if(xold.eq.xnew_sav)
        if(.not.l_oldnenewsav) then ! The last MC move was accepted, so 'xold' is the old 'xnew'
          zzposold = zzposnew
          iold_indices = inew_indices
        else   ! This branch should probably only happen at the beginning of a block
          if(iperiodic.eq.0) then !rings -> polar coords
            do iering=1,nelec
              zzposold(1,iering) = datan2(xold(2,iering),xold(1,iering)) !theta -> x coordinate
              zzposold(2,iering) = dsqrt(xold(1,iering)**2 + xold(2,iering)**2) ! r -> y coord
            enddo
          elseif(iperiodic.eq.1) then ! wires, keep coords
            zzposold = xold(1:2,:)
          endif
          ! Now, sort zzposold if we can't reuse an already sorted array
          ! iold_indices will contain information on the ordering of electrons in xold
          do i=1,nelec
            iold_indices(i) = i
          enddo
          lognb2 = int(dlog(dfloat(nelec))/dlog(2.D0)+1.D-14)
c   First, we need to sort the electron positions with respect to theta (for rings) or x (for wires)
c   I am using Shell sort for now; it might be more efficient to use something like quicksort or merge sort
c    This implementation of Shell sort is from Cyrus
          M = nelec
          do NN=1,lognb2
            M = M/2
            K = nelec - M
            do J = 1,K
              do I =J,1,-M
                L = I + M
                if(zzposold(1,L).gt.zzposold(1,I)) exit
                temppos = zzposold(:,I)
                itemp = iold_indices(I)
                zzposold(:,I) = zzposold(:,L)
                iold_indices(I) = iold_indices(L)
                zzposold(:,L) = temppos
                iold_indices(L) = itemp
              enddo
            enddo
          enddo
        endif
      endif

c Now, construct zzposnew so that it has about the same order as zzposold

      if(ielec.eq.0) then ! all-electron move; we need to find r and theta 
        if(iperiodic.eq.0) then !rings -> polar coords
          do iering=1,nelec
            ioi = iold_indices(iering)
            zzposnew(1,iering) = datan2(xnew(2,ioi),xnew(1,ioi)) 
            zzposnew(2,iering) = dsqrt(xnew(1,ioi)**2 + xnew(2,ioi)**2)
          enddo
        elseif(iperiodic.eq.1) then ! wires, keep coords
          zzposnew = xnew(1:2,iold_indices)
        endif
      else ! one-electron move
        zzposnew = zzposold ! start by setting new to old, since most of them will be the same
        if(iperiodic.eq.0) then
          temppos(1) = datan2(xnew(2,ielec),xnew(1,ielec))
          temppos(2) = dsqrt(xnew(1,ielec)**2 + xnew(2,ielec)**2)
        elseif(iperiodic.eq.1) then
          temppos = xnew(1:2,ielec)
        endif
        do i=1,nelec !insert the new position into the list at the right place
          if(iold_indices(i).eq.ielec) then
            zzposnew(:,i) = temppos
          endif
        enddo
      endif

c   Now, sort zzposnew.  Hopefully, this will go very quickly for 1-electron moves, since zzposnew 
c      will almost be sorted to begin with.
      
      inew_indices = iold_indices ! initally, zznewpos has same sorting as zzoldpos
      lognb2 = int(dlog(dfloat(nelec))/dlog(2.D0)+1.D-14)
      M = nelec
      do NN=1,lognb2
        M = M/2
        K = nelec - M
        do J = 1,K
          do I =J,1,-M
            L = I + M
            if(zzposnew(1,L).gt.zzposnew(1,I)) exit
            temppos = zzposnew(:,I)
            itemp = inew_indices(I)
            zzposnew(:,I) = zzposnew(:,L)
            inew_indices(I) = inew_indices(L)
            zzposnew(:,L) = temppos
            inew_indices(L) = itemp
          enddo
        enddo
      enddo

c  Now all of the electrons are sorted, and we can calculate observables
      
      zzsumold = 0.d0
      zzsumnew = 0.d0
      stagsignold = 1.0d0/dble(nelec)
      stagsignnew = 1.0d0/dble(nelec)
c     Set the sign of the staggered order such that the largest r (or y) has sign +1
c       i.e., in the zigzag phase, sum_i (-1)^i y_i should always be positive
      imaxold = maxloc(zzposold(2,:),1)
      imaxnew = maxloc(zzposnew(2,:),1)
      if (imaxold.eq.nelec) then
        imaxoldn = 1
      else
        imaxoldn = imaxold+1
      endif
      if (imaxnew.eq.nelec) then
        imaxnewn = 1
      else
        imaxnewn = imaxnew+1
      endif
      if(mod(imaxold,2).eq.0) stagsignold = -stagsignold
      if(mod(imaxnew,2).eq.0) stagsignnew = -stagsignnew
      rave = (q*sum(zzposold(2,:)) + p*sum(zzposnew(2,:)))/dble(nelec)
      do i =1,nelec
        if (iperiodic.eq.0) then
          zzmaglocal_old(i) = stagsignold*(zzposold(2,i)-rave)
          zzmaglocal_new(i) = stagsignnew*(zzposnew(2,i)-rave)
        else
          zzmaglocal_old(i) = stagsignold*zzposold(2,i)
          zzmaglocal_new(i) = stagsignnew*zzposnew(2,i)
        endif
        zzsumold = zzsumold + zzmaglocal_old(i)
        zzsumnew = zzsumnew + zzmaglocal_new(i)
        stagsignold = -stagsignold
        stagsignnew = -stagsignnew
      enddo
c  For debugging:
c      write(6,*) 'in zigzag2d:'
c      write(6,*) (zzposold(1,i),i=1,nelec)
c      write(6,*) (zzposold(2,i),i=1,nelec)
c      write(6,*) zzsumold, zzsumnew, q*dabs(zzsumold)+p*dabs(zzsumnew)
      zzterm(3) = q*zzsumold + p*zzsumnew
      zzterm(1) = q*dabs(zzsumold) + p*dabs(zzsumnew)
      zzterm(2) = q*zzsumold*zzsumold + p*zzsumnew*zzsumnew
c     Calculate the values if we throw out max value of y or r and its neighbor
      zzsumoldred = zzsumold-zzmaglocal_old(imaxold)-zzmaglocal_old(imaxoldn) 
      zzsumnewred = zzsumnew-zzmaglocal_new(imaxnew)-zzmaglocal_new(imaxnewn) 
      zzsumoldred = zzsumoldred*dble(nelec)/dble(nelec-2)
      zzsumnewred = zzsumnewred*dble(nelec)/dble(nelec-2)
      zzterm(6) = q*zzsumoldred + p*zzsumnewred
      zzterm(4) = q*dabs(zzsumoldred) + p*dabs(zzsumnewred)
      zzterm(5) = q*zzsumoldred*zzsumoldred + p*zzsumnewred*zzsumnewred
      
c     This is a kludge to make sure that the averages come out correctly 
c        for single-electron moves.  Since this routine gets called
c        once per electron in the mov1 update, we need to divide by
c        nelec. This is not needed for all-electron updates, though
c        since we just call this routine once after the update.
      corrnorm = dble(nelec) !makes sure corr is counted properly
      if(ielec.gt.0) then
        zzterm(:) = zzterm(:)/dble(nelec)
        corrnorm = 1.0
      endif
      zzsum(:) = zzsum(:) + zzterm(:)
c     'spread(v,dim,ncopies)' copies an array v, ncopies times along dim
      zzcorrmat_old = spread(zzmaglocal_old,dim=2,ncopies=nelec)*spread(zzmaglocal_old,dim=1,ncopies=nelec)
      zzcorrmat_new = spread(zzmaglocal_new,dim=2,ncopies=nelec)*spread(zzmaglocal_new,dim=1,ncopies=nelec)
      
      if(iperiodic.eq.0) then
        delxt = delti
      else
        delxt = delxi(1)
      endif

      do j = 0,nelec-1
        do i = 1,nelec
          i2 = mod(i+j-1,nelec) + 1  !mod returns a number in [0,n-1], array index is [1,n]
          ! compute difference in x or theta
          xtdiffo = zzposold(1,i2) - zzposold(1,i)
          xtdiffn = zzposnew(1,i2) - zzposnew(1,i)
          if(iperiodic.eq.1) then
            xtdiffo = modulo(xtdiffo,alattice)
            xtdiffn = modulo(xtdiffn,alattice)
            if (xtdiffo.ge.(alattice/2.)) xtdiffo = alattice - xtdiffo
            if (xtdiffn.ge.(alattice/2.)) xtdiffn = alattice - xtdiffn
          elseif(iperiodic.eq.0) then
            xtdiffo = modulo(xtdiffo,2.*3.1415926)
            xtdiffn = modulo(xtdiffn,2.*3.1415926)
            if (xtdiffo.ge.3.1415926) xtdiffo = 2.*3.1415926 - xtdiffo
            if (xtdiffn.ge.3.1415926) xtdiffn = 2.*3.1415926 - xtdiffn
          endif
          ixto = nint(delxt*xtdiffo)
          ixtn = nint(delxt*xtdiffn)
          zzcorrtermo = q*corrnorm*zzcorrmat_old(i,i2)
          zzcorrtermn = p*corrnorm*zzcorrmat_new(i,i2)
          zzcorr(ixto) = zzcorr(ixto) + zzcorrtermo
          zzcorr(ixtn) = zzcorr(ixtn) + zzcorrtermn
          zzcorrij(j) = zzcorr(j) + zzcorrtermo + zzcorrtermn
        enddo
      enddo
      
      xold_sav = xold
      xnew_sav = xnew
c      if(izigzag.gt.1) then ! do all of the pair density stuff
c      
c      endif
    
      return
      end

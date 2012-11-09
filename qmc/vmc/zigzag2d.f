      subroutine zigzag2d(p,q,xold,xnew,ielec)

c Written by Abhijit Mehta, December 2011
c  Calculates quantities useful for studying zigzag quantum phase
c  transition in rings and wires.  
c  -reduced pair density
c  -"staggered amplitude"
c  - zigzag correlation functions
c  xold and xnew are the old and new configurations
c  p and q are the probabilities of accept and reject
c  ielec labels the electron that was moved for 1-electron moves
c  if ielec=0, we are doing an all-electron move.

c  If you want to add a new observable, make sure you change
c    nzzvars in zigzag_mod.f90 (and do a make clean!)
c  You will also need to add a print out line to print_zigzag_vars()
c       below
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
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      dimension xold(3,nelec),xnew(3,nelec)
      dimension zztemp2(nelec), ransign(nelec)
      dimension iwouter_old(nelec), iwouter_new(nelec)
      dimension evenfluct_old((nelec+1)/2), outerfluct_old((nelec+1)/2)
      dimension evenfluct_new((nelec+1)/2), outerfluct_new((nelec+1)/2)
      dimension temppos(2),zzterm(nzzvars)
      dimension zzmaglocal_new(nelec),zzmaglocal_old(nelec)
      dimension zzcorrmat_old(nelec,nelec),zzcorrmat_new(nelec,nelec)
      dimension cnndiffo(nelec), cn2ndiffo(nelec), cnndiffn(nelec), cn2ndiffn(nelec)
      dimension cnndiffo2(nelec), cn2ndiffo2(nelec), cnndiffn2(nelec), cn2ndiffn2(nelec)
      dimension cnndiffo3(nelec), cn2ndiffo3(nelec), cnndiffn3(nelec), cn2ndiffn3(nelec)
      dimension cnndiffo4(nelec), cn2ndiffo4(nelec), cnndiffn4(nelec), cn2ndiffn4(nelec)
      dimension zzxthetaold(nelec), zzxthetanew(nelec)

      if(ielec.lt.0 .or. ielec.gt.nelec) then
        write (6,*) 'Bad value of ielec in zigzag2d, ', ielec
        stop "Bad value of ielec in zigzag2d"
      endif

      if(iperiodic.eq.0) then
        delxti = delti
      else
        delxti = delxi(1)
      endif
      delyri = 1.0/zzdelyr

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
c     Set the sign of the staggered order such that the n/3rd largest r (or y) has sign +1
c       i.e., in the zigzag phase, sum_i (-1)^i y_i should always be positive
c       Choosing the n/3rd (and not the largest) should hopefully cause zigzag amp = 0 in linear phase
c     Also, see which electrons are in outer half, which are in inner half
c   iwouter(i) is 1 if ith electron has one of the n/2 largest radii, -1 otherwise


      iwouter_old(:) = -1
      iwouter_new(:) = -1
      zztemp2(:) = zzposold(2,:)
      do ipos = 1,nelec/2
        iztemp = maxloc(zztemp2,1)
        zztemp2(iztemp) = -1.0
        iwouter_old(iztemp) = 1
        if (ipos.le.nelec/3) izagold = iztemp
      enddo
      zztemp2(:) = zzposnew(2,:)
      do ipos = 1,nelec/2
        iztemp = maxloc(zztemp2,1)
        zztemp2(iztemp) = -1.0
        iwouter_new(iztemp) = 1
        if (ipos.le.nelec/3) izagnew = iztemp
      enddo
      if(mod(izagold,2).eq.0) stagsignold = -stagsignold
      if(mod(izagnew,2).eq.0) stagsignnew = -stagsignnew
c      rave = (q*sum(zzposold(2,:)) + p*sum(zzposnew(2,:)))/dble(nelec)

c     Calculate <theta_i+1 - theta_i> and <theta_i+2 - theta_i>
      zzxthetaold = zzposold(1,:)
      zzxthetanew = zzposnew(1,:)
      cnndiffo = cshift(zzxthetaold, SHIFT = 1) - zzxthetaold
      cn2ndiffo = cshift(zzxthetaold, SHIFT = 2) - zzxthetaold
      cnndiffn = cshift(zzxthetanew, SHIFT = 1) - zzxthetanew
      cn2ndiffn = cshift(zzxthetanew, SHIFT = 2) - zzxthetanew
      if (iperiodic.eq.0) then
        cellsize = 2.*pi
      else if (iperiodic.eq.1) then
        cellsize = alattice
      endif
  
      ioutero = 0   
      ieveno = 0
      ioutern = 0   
      ievenn = 0
c     xtspacing = cellsize/dble(nelec)
c     xtshift = cellsize/2. - xtspacing/2. !since theta/x go from -pi..pi or -alattice/2..a/2
      xtspacing = cellsize/dble(nelec/2)
      xtshift = cellsize/2. - xtspacing/2. !since theta/x go from -pi..pi or -alattice/2..a/2
      do i =1,nelec
        if (iperiodic.eq.0) then
          zzmaglocal_old(i) = stagsignold*(zzposold(2,i)-rring)
          zzmaglocal_new(i) = stagsignnew*(zzposnew(2,i)-rring)
c          zzmaglocal_old(i) = stagsignold*(zzposold(2,i)-rave)
c          zzmaglocal_new(i) = stagsignnew*(zzposnew(2,i)-rave)
        else
          zzmaglocal_old(i) = stagsignold*zzposold(2,i)
          zzmaglocal_new(i) = stagsignnew*zzposnew(2,i)
        endif
        zzsumold = zzsumold + zzmaglocal_old(i)
        zzsumnew = zzsumnew + zzmaglocal_new(i)
        if(stagsignold.gt.0) then
          ieveno = ieveno + 1
c         evenfluct_old(ieveno) = (zzposold(1,i)- xtspacing*i + xtshift)/dble(nelec)
          evenfluct_old(ieveno) = (zzposold(1,i)- xtspacing*ieveno + xtshift)/dble(nelec)
        endif 
        if(stagsignnew.gt.0) then
          ievenn = ievenn + 1
c         evenfluct_new(ievenn) = (zzposnew(1,i)- xtspacing*i + xtshift)/dble(nelec)
          evenfluct_new(ievenn) = (zzposnew(1,i)- xtspacing*ievenn + xtshift)/dble(nelec)
        endif 
        if(iwouter_old(i).gt.0) then
          ioutero = ioutero + 1
c         outerfluct_old(ioutero) = (zzposold(1,i)- xtspacing*i + xtshift)/dble(nelec)
          outerfluct_old(ioutero) = (zzposold(1,i)- xtspacing*ioutero + xtshift)/dble(nelec)
        endif 
        if(iwouter_new(i).gt.0) then
          ioutern = ioutern + 1
c         outerfluct_new(ioutern) = (zzposnew(1,i)- xtspacing*i + xtshift)/dble(nelec)
          outerfluct_new(ioutern) = (zzposnew(1,i)- xtspacing*ioutern + xtshift)/dble(nelec)
        endif 
        stagsignold = -stagsignold
        stagsignnew = -stagsignnew
        cnndiffo(i) = dabs(cnndiffo(i))
        cn2ndiffo(i) = dabs(cn2ndiffo(i))
        cnndiffn(i) = dabs(cnndiffn(i))
        cn2ndiffn(i) = dabs(cn2ndiffn(i))
        if (cnndiffo(i).ge.(cellsize/2.)) cnndiffo(i) = cellsize - cnndiffo(i)
        if (cn2ndiffo(i).ge.(cellsize/2.)) cn2ndiffo(i) = cellsize - cn2ndiffo(i)
        if (cnndiffn(i).ge.(cellsize/2.)) cnndiffn(i) = cellsize - cnndiffn(i)
        if (cn2ndiffn(i).ge.(cellsize/2.)) cn2ndiffn(i) = cellsize - cn2ndiffn(i)
        ixto_nn = nint(delxti*cnndiffo(i)*10./dble(nelec))
        ixto_n2n = nint(delxti*cn2ndiffo(i)*10./dble(nelec))
        ixtn_nn = nint(delxti*cnndiffn(i)*10./dble(nelec))
        ixtn_n2n = nint(delxti*cn2ndiffn(i)*10./dble(nelec))
        znncorr(ixto_nn) = znncorr(ixto_nn) + q/dble(nelec)
        znncorr(ixtn_nn) = znncorr(ixtn_nn) + p/dble(nelec)
        zn2ncorr(ixto_n2n) = zn2ncorr(ixto_n2n) + q/dble(nelec)
        zn2ncorr(ixtn_n2n) = zn2ncorr(ixtn_n2n) + p/dble(nelec)
      enddo
c  For debugging:
c      write(6,*) 'in zigzag2d:'
c      write(6,*) (zzposold(1,i),i=1,nelec)
c      write(6,*) (zzposold(2,i),i=1,nelec)
c      write(6,*) (zzposnew(1,i),i=1,nelec)
c      write(6,*) (zzposnew(2,i),i=1,nelec)
c      write(6,*) zzsumold, zzsumnew, q*dabs(zzsumold)+p*dabs(zzsumnew)
c      write(6,*) zzsumold, zzsumnew, q*dabs(zzsumold)+p*dabs(zzsumnew)
c  For odd N, there's an extra rave.  We should correct this.
c   The following 2 lines are to test n=29, rs=3.7, w=0.1 rings
c      zzsumold = zzsumold + (35.486-rring)/29.
c      zzsumnew = zzsumnew + (35.486-rring)/29. 
      zzterm(3) = q*zzsumold + p*zzsumnew
      zzterm(1) = q*dabs(zzsumold) + p*dabs(zzsumnew)
      zzterm(2) = q*zzsumold*zzsumold + p*zzsumnew*zzsumnew
      zzterm(10) = q*(zzsumold**4) + p*(zzsumnew**4)
c     Calculate the values if we throw out max value of y or r and its neighbor
      imaxold = maxloc(zzposold(2,:),1)
      imaxnew = maxloc(zzposnew(2,:),1)
      iminold = minloc(zzposold(2,:),1)
      iminnew = minloc(zzposnew(2,:),1)
c     if (imaxold.eq.nelec) then
c       imaxoldn = 1
c     else
c       imaxoldn = imaxold+1
c     endif
c     if (imaxnew.eq.nelec) then
c       imaxnewn = 1
c     else
c       imaxnewn = imaxnew+1
c     endif
c      zzsumoldred = zzsumold-zzmaglocal_old(imaxold)-zzmaglocal_old(imaxoldn) 
c      zzsumnewred = zzsumnew-zzmaglocal_new(imaxnew)-zzmaglocal_new(imaxnewn) 
      zzsumoldred = zzsumold-zzmaglocal_old(imaxold)-zzmaglocal_old(iminold) 
      zzsumnewred = zzsumnew-zzmaglocal_new(imaxnew)-zzmaglocal_new(iminnew) 
      zzsumoldred = zzsumoldred*dble(nelec)/dble(nelec-2)
      zzsumnewred = zzsumnewred*dble(nelec)/dble(nelec-2)
      zzterm(6) = q*zzsumoldred + p*zzsumnewred
      zzterm(4) = q*dabs(zzsumoldred) + p*dabs(zzsumnewred)
      zzterm(5) = q*zzsumoldred*zzsumoldred + p*zzsumnewred*zzsumnewred
      zzterm(11) = q*(zzsumoldred**4) + p*(zzsumnewred**4)

c     Averages of <x_i+1 - x_i>, <x_i+2 - x_i>, and higher moments
      cnndiffo2(:) = cnndiffo(:)*cnndiffo(:)
      cnndiffo3(:) = cnndiffo2(:)*cnndiffo(:)
      cnndiffo4(:) = cnndiffo3(:)*cnndiffo(:)
      cnndiffn2(:) = cnndiffn(:)*cnndiffn(:)
      cnndiffn3(:) = cnndiffn2(:)*cnndiffn(:)
      cnndiffn4(:) = cnndiffn3(:)*cnndiffn(:)
      cn2ndiffo2(:) = cn2ndiffo(:)*cn2ndiffo(:)
      cn2ndiffo3(:) = cn2ndiffo2(:)*cn2ndiffo(:)
      cn2ndiffo4(:) = cn2ndiffo3(:)*cn2ndiffo(:)
      cn2ndiffn2(:) = cn2ndiffn(:)*cn2ndiffn(:)
      cn2ndiffn3(:) = cn2ndiffn2(:)*cn2ndiffn(:)
      cn2ndiffn4(:) = cn2ndiffn3(:)*cn2ndiffn(:)
      
      zzterm(16) = q*sum(cnndiffo(:)) + p*sum(cnndiffn(:))
      zzterm(17) = q*sum(cnndiffo2(:)) + p*sum(cnndiffn2(:))
      zzterm(18) = q*sum(cnndiffo3(:)) + p*sum(cnndiffn3(:))
      zzterm(19) = q*sum(cnndiffo4(:)) + p*sum(cnndiffn4(:))
      zzterm(20) = q*sum(cn2ndiffo(:)) + p*sum(cn2ndiffn(:))
      zzterm(21) = q*sum(cn2ndiffo2(:)) + p*sum(cn2ndiffn2(:))
      zzterm(22) = q*sum(cn2ndiffo3(:)) + p*sum(cn2ndiffn3(:))
      zzterm(23) = q*sum(cn2ndiffo4(:)) + p*sum(cn2ndiffn4(:))
      zzterm(16:23) = zzterm(16:23)/dble(nelec) 


c     Pick sign randomly, so that N/2 have "-" sign
      ransign(:) = 1.0d0/dble(nelec)
      do itry = 1,nelec/2
        do 
          irand = int(nelec*rannyu(0)) + 1
          if (ransign(irand).gt.0) then
            ransign(irand) = -ransign(irand)
            exit
          endif
        enddo
      enddo
      zzrandsumold = sum(zzposold(2,:)*ransign(:))
      zzrandsumnew = sum(zzposnew(2,:)*ransign(:))
      zzterm(9) = q*zzrandsumold + p*zzrandsumnew
      zzterm(7) = q*dabs(zzrandsumold) + p*dabs(zzrandsumnew)
      zzterm(8) = q*zzrandsumold*zzrandsumold + p*zzrandsumnew*zzrandsumnew
      zzterm(12) = q*(zzrandsumold**4) + p*(zzrandsumnew**4)
     
      zzevenout_numerator_new = sum(evenfluct_new(:)*outerfluct_new(:)) ! don't need to divide by n/2
      zzevenout_numerator_old = sum(evenfluct_old(:)*outerfluct_old(:)) !  bc numerator and denominator
      zzevenout_evensigma_old = sum(evenfluct_old(:)*evenfluct_old(:))  !  cancel out.
      zzevenout_evensigma_new = sum(evenfluct_new(:)*evenfluct_new(:))
      zzevenout_outersigma_old = sum(outerfluct_old(:)*outerfluct_old(:))
      zzevenout_outersigma_new = sum(outerfluct_new(:)*outerfluct_new(:))
!      zzevenout_denominator_new = dsqrt(zzevenout_evensigma_new*zzevenout_outersigma_new)
!      zzevenout_denominator_old = dsqrt(zzevenout_evensigma_old*zzevenout_outersigma_old)

!      zzterm(13) = q*zzevenout_numerator_old/zzevenout_denominator_old + p*zzevenout_numerator_new/zzevenout_denominator_new

      zzterm(13) = q*zzevenout_numerator_old + p*zzevenout_numerator_new
      zzterm(14) = q*zzevenout_evensigma_old + p*zzevenout_evensigma_new
      zzterm(15) = q*zzevenout_outersigma_old + p*zzevenout_outersigma_new

c     This is a kludge to make sure that the averages come out correctly 
c        for single-electron moves.  Since this routine gets called
c        once per electron in the mov1 update, we need to divide by
c        nelec. This is not needed for all-electron updates, though
c        since we just call this routine once after the update.
      corrnorm = dble(nelec) !makes sure corr is counted properly
      pairdennorm = 1.0d0
      if(ielec.gt.0) then
        zzterm(:) = zzterm(:)/dble(nelec)
        corrnorm = 1.0 ! remember that zzmaglocal^2 has a factor of 1/nelec^2 in it!!
        pairdennorm = 1.0d0/dble(nelec)
      endif
      zzsum(:) = zzsum(:) + zzterm(:)
c     'spread(v,dim,ncopies)' copies an array v, ncopies times along dim
c     zzcorrmat_old = spread(zzmaglocal_old,dim=2,ncopies=nelec)*spread(zzmaglocal_old,dim=1,ncopies=nelec)
c     zzcorrmat_new = spread(zzmaglocal_new,dim=2,ncopies=nelec)*spread(zzmaglocal_new,dim=1,ncopies=nelec)
      corrfactor = (corrnorm/dble(nelec))/dble(nelec)
      corrsgn = corrfactor
      do j = 0,nelec-1
        do i = 1,nelec
c          i2 = mod(i+j-1,nelec) + 1  !mod returns a number in [0,n-1], array index is [1,n]
          i2 = i + j
          if (i2.gt.nelec) i2 = i2 - nelec
          ! compute difference in x or theta
          xtdiffo = zzposold(1,i2) - zzposold(1,i)
          xtdiffn = zzposnew(1,i2) - zzposnew(1,i)
          if(iperiodic.eq.1) then
            xtdiffo = modulo(xtdiffo,alattice)
            xtdiffn = modulo(xtdiffn,alattice)
            if (xtdiffo.ge.(alattice/2.d0)) xtdiffo = alattice - xtdiffo
            if (xtdiffn.ge.(alattice/2.d0)) xtdiffn = alattice - xtdiffn
          elseif(iperiodic.eq.0) then
            xtdiffo = modulo(xtdiffo,2.d0*3.1415926d0)
            xtdiffn = modulo(xtdiffn,2.d0*3.1415926d0)
            if (xtdiffo.ge.3.1415926d0) xtdiffo = 2.d0*3.1415926d0 - xtdiffo
            if (xtdiffn.ge.3.1415926d0) xtdiffn = 2.d0*3.1415926d0 - xtdiffn
          endif
          ixto = nint(delxti*xtdiffo)
          ixtn = nint(delxti*xtdiffn)
c         zzcorrtermo = q*corrnorm*zzcorrmat_old(i,i2)
c         zzcorrtermn = p*corrnorm*zzcorrmat_new(i,i2)
c          zzcorrtermo = q*corrnorm*zzmaglocal_old(i)*zzmaglocal_old(i2)
c          zzcorrtermn = p*corrnorm*zzmaglocal_new(i)*zzmaglocal_new(i2)
c         <-1^(i+j) (r_i - rave)*(r_j - rave)> = 
c                <-1^(i-j)(r_i*r_j)> - rave*<-1^(i-j)(r_i + r_j)> + rave^2*<-1^(i-j)>
c         label terms:   zzcorr      - rave*zzcorr1               + rave^2*zzcorr2
c         For yycorr, we do the same thing (only without the -1^(i-j)) 
          zzcorrtermo = q*corrsgn*zzposold(2,i)*zzposold(2,i2)
          zzcorrtermn = p*corrsgn*zzposnew(2,i)*zzposnew(2,i2)
          zzcorr(ixto) = zzcorr(ixto) + zzcorrtermo
          zzcorr(ixtn) = zzcorr(ixtn) + zzcorrtermn
          zzcorrterm1o = q*corrsgn*(zzposold(2,i) + zzposold(2,i2))
          zzcorrterm1n = p*corrsgn*(zzposnew(2,i) + zzposnew(2,i2))
          zzcorr1(ixto) = zzcorr1(ixto) + zzcorrterm1o
          zzcorr1(ixtn) = zzcorr1(ixtn) + zzcorrterm1n
          zzcorrterm2o = q*corrsgn
          zzcorrterm2n = p*corrsgn
          zzcorr2(ixto) = zzcorr2(ixto) + zzcorrterm2o
          zzcorr2(ixtn) = zzcorr2(ixtn) + zzcorrterm2n
          zzcorrij(j) = zzcorrij(j) + zzcorrtermo + zzcorrtermn
          yycorrtermo = q*corrfactor*zzposold(2,i)*zzposold(2,i2)
          yycorrtermn = p*corrfactor*zzposnew(2,i)*zzposnew(2,i2)
          yycorr(ixto) = yycorr(ixto) + yycorrtermo
          yycorr(ixtn) = yycorr(ixtn) + yycorrtermn
          yycorrterm1o = q*corrfactor*(zzposold(2,i) + zzposold(2,i2))
          yycorrterm1n = p*corrfactor*(zzposnew(2,i) + zzposnew(2,i2))
          yycorr1(ixto) = yycorr1(ixto) + yycorrterm1o
          yycorr1(ixtn) = yycorr1(ixtn) + yycorrterm1n
          yycorrterm2o = q*corrfactor
          yycorrterm2n = p*corrfactor
          yycorr2(ixto) = yycorr2(ixto) + yycorrterm2o
          yycorr2(ixtn) = yycorr2(ixtn) + yycorrterm2n
          yycorrij(j) = yycorrij(j) + yycorrtermo + yycorrtermn
          if(izigzag.gt.1 .and. j.ne.0) then ! do all of the pair density stuff
            yrdiffo = zzposold(2,i2) - zzposold(2,i)
            yrdiffn = zzposnew(2,i2) - zzposnew(2,i)
            iyro = min(max(nint(delyri*yrdiffo),-NAX),NAX)
            iyrn = min(max(nint(delyri*yrdiffn),-NAX),NAX)
            zzpairden_t(iyro,ixto) = zzpairden_t(iyro,ixto) + q*pairdennorm*0.5
            zzpairden_t(iyro,-ixto) = zzpairden_t(iyro,-ixto) + q*pairdennorm*0.5
            zzpairdenij_t(iyro,j) = zzpairdenij_t(iyro,j) + q*pairdennorm
            zzpairden_t(iyrn,ixtn) = zzpairden_t(iyrn,ixtn) + p*pairdennorm*0.5
            zzpairden_t(iyrn,-ixtn) = zzpairden_t(iyrn,-ixtn) + p*pairdennorm*0.5
            zzpairdenij_t(iyrn,j) = zzpairdenij_t(iyrn,j) + p*pairdennorm
          endif
        enddo
        corrsgn = -corrsgn
      enddo
      
      xold_sav = xold
      xnew_sav = xnew

    
      return
      end

c-------------------------------------------------------------------

      subroutine print_zigzag_vars(zzave,zzerr,rtpass)

c     Written by Abhijit Mehta, May 2012
c      Routine to print out all the zigzag variables in the
c       various finwrt routines.
c      Inputs: 
c       zzave - array of size nzzvars with averages of all zigzag vars
c       zzerr - array of size nzzvars with errors in all of the averages
c       rtpass - sqrt of the number of passes
c
c        We put all the print out statements here since this code was
c         basically repeated several times throughout CHAMP
c        Now, if we add a new variable, we only need to change the above
c        subroutine (zigzag2d), this subroutine, and the parameter
c        'nzzvars' in zigzag_mod.f90
      
      use zigzag_mod, only: nzzvars

      implicit real*8(a-h,o-z)
      dimension zzave(nzzvars), zzerr(nzzvars)

c     Fourth cumulants:
      zzu4 = 1.0d0 - zzave(10)/(3.0d0*zzave(2)*zzave(2))
      zzu4err = 1.0d0/(3.0d0*zzave(2)*zzave(2)) * ((2.0d0*zzave(10)/zzave(2))*zzerr(2) - zzerr(10))

      zzu4red = 1.0d0 - zzave(11)/(3.0d0*zzave(5)*zzave(5))
      zzu4rederr = 1.0d0/(3.0d0*zzave(5)*zzave(5)) * ((2.0d0*zzave(11)/zzave(5))*zzerr(5) - zzerr(11))

      zzu4rand = 1.0d0 - zzave(12)/(3.0d0*zzave(8)*zzave(8))
      zzu4randerr = 1.0d0/(3.0d0*zzave(8)*zzave(8)) * ((2.0d0*zzave(12)/zzave(8))*zzerr(8) - zzerr(12))

c  Even-outer correlation

      zzeocorr = zzave(13)/dsqrt(zzave(14)*zzave(15))

c  This line is in the finwrt routines:      
c     write(6,'(''physical variable'',t20,''average'',t34,''rms error''
c    &,t47,''rms er*rt(pass)'',t65,''sigma'',t86,''Tcor'')')  !JT

      write(6,'(''<ZigZag Amp> ='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(3),zzerr(3),zzerr(3)*rtpass
      write(6,'(''<|ZigZag Amp|> ='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(1),zzerr(1),zzerr(1)*rtpass
      write(6,'(''<ZigZag Amp^2> ='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(2),zzerr(2),zzerr(2)*rtpass
      write(6,'(''<ZigZag Amp (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(6),zzerr(6),zzerr(6)*rtpass
      write(6,'(''<|ZigZag Amp| (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(4),zzerr(4),zzerr(4)*rtpass
      write(6,'(''<ZigZag Amp^2 (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(5),zzerr(5),zzerr(5)*rtpass
      write(6,'(''<ZigZag rand Amp>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(9),zzerr(9),zzerr(9)*rtpass
      write(6,'(''<|ZigZag rand Amp|>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(7),zzerr(7),zzerr(7)*rtpass
      write(6,'(''<ZigZag rand Amp^2>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(8),zzerr(8),zzerr(8)*rtpass
      write(6,'(''<ZigZag Amp^4> ='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(10),zzerr(10),zzerr(10)*rtpass
      write(6,'(''<ZigZag Amp^4 (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(11),zzerr(11),zzerr(11)*rtpass
      write(6,'(''<ZigZag rand Amp^4>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(12),zzerr(12),zzerr(12)*rtpass
      write(6,'(''<U4> = '',t22,f12.7,'' +-'',f11.7,f9.5)') zzu4, zzu4err, zzu4err*rtpass
      write(6,'(''<U4 (red)> = '',t22,f12.7,'' +-'',f11.7,f9.5)') zzu4red, zzu4rederr, zzu4rederr*rtpass
      write(6,'(''<U4 rand> = '',t22,f12.7,'' +-'',f11.7,f9.5)') zzu4rand, zzu4randerr, zzu4randerr*rtpass
      write(6,'(''<Even-Outer Covariance>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(13),zzerr(13),zzerr(13)*rtpass
      write(6,'(''<Even-elec fluct^2>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(14),zzerr(14),zzerr(14)*rtpass
      write(6,'(''<Outer-elec fluct^2>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(15),zzerr(15),zzerr(15)*rtpass
      write(6,'(''<Even-Outer Corr>='',t22,f12.7)') zzeocorr
      write(6,'(''<th_(i+1) - th_i>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(16),zzerr(16),zzerr(16)*rtpass
      write(6,'(''<(th_(i+1) - th_i)^2>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(17),zzerr(17),zzerr(17)*rtpass
      write(6,'(''<(th_(i+1) - th_i)^3>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(18),zzerr(18),zzerr(18)*rtpass
      write(6,'(''<(th_(i+1) - th_i)^4>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(19),zzerr(19),zzerr(19)*rtpass
      write(6,'(''<th_(i+2) - th_i>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(20),zzerr(20),zzerr(20)*rtpass
      write(6,'(''<(th_(i+2) - th_i)^2>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(21),zzerr(21),zzerr(21)*rtpass
      write(6,'(''<(th_(i+2) - th_i)^3>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(22),zzerr(22),zzerr(22)*rtpass
      write(6,'(''<(th_(i+2) - th_i)^4>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(23),zzerr(23),zzerr(23)*rtpass
      
      return
      end

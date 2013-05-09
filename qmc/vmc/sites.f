      subroutine sites(x,nelec,nsite)
! Written by Cyrus Umrigar
      use constants_mod
      use atom_mod
      use dim_mod
      use pseudo_mod
      use jel_sph2_mod
      use contrl_per_mod
      use orbpar_mod
      use wfsec_mod
      use dorb_mod
      implicit real*8(a-h,o-z)

! Routine to put electrons down around centers for a VERY crude initial
! configuration if nothing else is available.  It is better to put them
! too close than to put them too far away because they equilibrate faster
! when they are too close.

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
       common /wire/ wire_w,wire_length,wire_length2,wire_radius2, wire_potential_cutoff,wire_prefactor,wire_root1

      dimension x(3,*),nsite(*)

! Loop over spins and centers. If odd number of electrons on all
! atoms then the up-spins have an additional electron.
! So assumption is that system is not strongly polarized.

!     gauss()=dcos(two*pi*rannyu(0))*dsqrt(-two*dlog(rannyu(0)))
      pi=4*datan(1.d0)

      if((nloc.eq.-1).or.(nloc.eq.-5)) then ! parabolic quantum dot
        if(we.eq.0.d0) stop 'we should not be 0 in sites for quantum dots (nloc=-1)'
       elseif(nloc.eq.-3) then ! jellium RM
        if(zconst.eq.0.d0) stop 'zconst should not be 0 in sites for atoms in jellium (nloc=-3)'
      endif

      ielec=0
      do 10 ispin=1,2
        do 10 i=1,ncent
          if((nloc.eq.-1).or.(nloc.eq.-5)) then ! parabolic quantum dot
            znucc=dsqrt(we)
           elseif(nloc.eq.-3) then ! jellium RM
            znucc=zconst
           elseif(nloc.eq.-4) then ! quantum wire
            znucc=dsqrt(wire_w)
           else ! atoms and molecules
            if(znuc(iwctype(i)).eq.0.d0)
     &stop 'znuc should not be 0 in sites for atoms and molecules'
            znucc=znuc(iwctype(i))
          endif
          ju=(nsite(i)+2-ispin)/2
          do 10 j=1,ju
            ielec=ielec+1
            if(ielec.gt.nelec) return
            if(nloc.eq.-1 .or. nloc.eq.-5 .or. nloc.eq.-4) then
              sitsca=1/znucc
             elseif(j.eq.1) then
              sitsca=1/max(znucc,1.d0)
             elseif(j.le.5) then
              sitsca=2/max(znucc-2,1.d0)
             elseif(j.le.9) then
              sitsca=3/max(znucc-10,1.d0)
             elseif(j.le.18) then
              sitsca=4/max(znucc-18,1.d0)
             else
              sitsca=5/max(znucc-36,1.d0)
            endif


! sample position from exponentials or gaussian around center
! A.D.Guclu 5/2008: need circular coo. for ring shaped quantum dots
            if((nloc.eq.-1 .or. nloc.eq.-5) .and. rring.gt.0.d0) then
              if(ibasis.eq.5) then
                site = (0.5d0 - rannyu(0))/dsqrt(we*oparm(3, iworbd(ielec,1), iwf))
                angle = (0.5d0 - rannyu(0))/dsqrt(oparm(4, iworbd(ielec,1), iwf))
                site = site + oparm(1, iworbd(ielec,1), iwf)
                angle = angle + oparm(2, iworbd(ielec,1), iwf)
!  Make sure electron is near the center of some gaussian - might not work
!     if there's more than 1 slater determinant
                x(1,ielec)=site*dcos(angle)
                x(2,ielec)=site*dsin(angle)
              else
!               This code sampled from a gaussian:
!                site=-dlog(rannyu(0))
!                site=dsqrt(site)
!                site=sign(site,(rannyu(0)-half))
!               This code samples from a smaller, uniform region:
!                site = 2.0d0*(0.5d0 - rannyu(0))
                site = (0.5d0 - rannyu(0))/dsqrt(we)
                angle=2.0d0*pi*(dble(ielec) - rannyu(0))/dble(nelec)
                x(1,ielec)=(site+rring)*dcos(angle)
                x(2,ielec)=(site+rring)*dsin(angle)
!               x(1,ielec)=(sitsca*site+rring)*dcos(angle)
!               x(2,ielec)=(sitsca*site+rring)*dsin(angle)
              endif
             else
               do 5 k=1,ndim
! sample position from exponentials or gaussian around center
! a.d.guclu: for wires distribute electrons linearly in y direction
! a.c.mehta: unless floating gaussians, then make sure electrons
!             are close to centers of gaussians
!  Warning:  this might not work if we have multiple slater determinants
                 site=-dlog(rannyu(0))
                 if(nloc.eq.-1 .or. nloc.eq.-4 .or. nloc.eq.-5) site=dsqrt(site)
                 site=sign(site,(rannyu(0)-half))

                 if(nloc.eq.-4) then
                   if (ibasis.eq.6 .or. ibasis.eq.7) then
                     site = (0.5d0 - rannyu(0))/dsqrt(we*oparm(k+2, iworbd(ielec,1), iwf))
!  Make sure electron is near the center of some gaussian - might not work
!     if there's more than 1 slater determinant
                     x(k,ielec) = site + oparm(k, iworbd(ielec,1), iwf)
                   else
                     if(k.eq.2) then
                       x(k,ielec)=sitsca*site+cent(k,i)
                     elseif(iperiodic.eq.0) then
                       x(k,ielec)=wire_length*(0.5d0-rannyu(0))
                     else
                       x(k,ielec)=wire_length*rannyu(0)
                     endif
                   endif
                 else
                   x(k,ielec)=sitsca*site+cent(k,i)
                 endif
   5           enddo
            endif

   10     continue


!      write(6,*)
      write(6,'(a,i3,a)') '1 configuration for',ielec,' electrons has been generated by routine sites.'
      write(6,'(''sites:'',1000d12.4)') ((x(k,i),k=1,ndim),i=1,nelec)

      if(ielec.lt.nelec) stop 'bad input to sites'
      return
      end


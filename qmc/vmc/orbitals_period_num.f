      subroutine orbitals_period_num(x,orb,dorb,ddorb)
! Written by Cyrus Umrigar.  Modified by William Parker to interface to Einspline.
! Calculate pw orbitals, gradient and laplacian
! by interpolating on a grid.
      use bwfdet_mod
      use bsplines_mod
      use orbital_grid_mod
      use coefs_mod
      use dets_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use periodic_mod
      use periodic2_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*),r_basis(3)
      dimension orb(nelec,*),dorb(3,nelec,*),ddorb(nelec,*)
      dimension orb_tmp(norb),dorb_tmp(3,norb),ddorb_tmp(norb)

      dimension orb_blip_tmp(norb,ndet),dorb_blip_tmp(3,norb,ndet),ddorb_blip_tmp(norb,ndet)


!wparker Set up the ict array to tell PSPLINE to return the value
!        of the function, its first derivative in every direction
!        and its second derivative in every direction and the mixed
!        derivatives

!     do i=1,10
!       ict(i)=1
!     enddo

      do iel=1,nelec

! Determine position in lattice coordinates in line 10, same modulo 1 in line 20, and whether phase factor, isgn, is +1 or -1.
! Note that modulo and mod are different functions in fortran.
! Note we add 1 to r_basis because interpol_orb expects a positive input
! The 'if' after 10 does the foll.  Since the k-vectors can be half sim-cell reciprocal lattice vectors
! and we reduce r_basis(k) to be between 0 and 1 in line 20, we need to multiply by -1 for
! each rbasis(k) in line 10 for which floor(r_basis) is odd.
! r_basis in the line just before 20 is between  0 and 1 if r_basis in line 10 is > 0
!                                   but between -1 and 0 if r_basis in line 10 is < 0
! r_basis in line 20 is always between 0 and 1.
        isgn=1
        do 20 k=1,ndim
          r_basis(k)=0
          do 10 i=1,ndim
   10       r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*x(i,iel)
          if(rkvec_shift_latt(k).ne.0.d0) then
            if(r_basis(k).ge.0.d0) then
              isgn=isgn*(-1)**int(r_basis(k))
             else
              isgn=isgn*(-1)**(int(r_basis(k))+1)
            endif
          endif
          r_basis(k)=r_basis(k)-int(r_basis(k))
   20     if(r_basis(k).lt.0.d0) r_basis(k)=r_basis(k)+1

! Warning: I could replace the above with the foll., but since I do not have the time to test it I am
! for the moment just putting in the foll. commented lines:
! Determine position in lattice coordinates in line 10, same modulo 1 in line 20, and whether phase factor, isgn, is +1 or -1.
! Note that modulo and mod are different functions in fortran.
! Since the k-vectors can be half sim-cell reciprocal lattice vectors
! and we reduce r_basis(k) to be between 0 and 1 in line 20, we need to multiply by -1 for
! each rbasis(k) in line 10 for which floor(r_basis) is odd.
!       isgn=1
!       do 20 k=1,ndim
!         r_basis(k)=0
!         do 10 i=1,ndim
!  10       r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*x(i,iel)
!         isgn=isgn*(-1)**floor(r_basis(k))
!  20     r_basis(k)=r_basis(k)-floor(r_basis(k))


        if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
          xi=r_basis(1)*ngrid_orbx
          yi=r_basis(2)*ngrid_orby
          zi=r_basis(3)*ngrid_orbz

!       write(6,'(''r_basis'',9f9.4)') r_basis,xi,yi,zi

! Warning:  If we changed the order of iel and iorb then we would not need _tmp arrays
!       call interpol_orb(ngrid_orbx,ngrid_orby,ngrid_orbz,xi,yi,zi,orb(1,iel),dorb(1,1,iel),ddorb(1,iel))
        call interpol_orb(ngrid_orbx,ngrid_orby,ngrid_orbz,xi,yi,zi,orb_tmp,dorb_tmp,ddorb_tmp)

          do iorb=1,norb
            orb(iel,iorb)=orb_tmp(iorb)*isgn
            do k=1,ndim
              dorb(k,iel,iorb)=dorb_tmp(k,iorb)*isgn
            enddo
            ddorb(iel,iorb)=ddorb_tmp(iorb)*isgn
!            if((iorb.eq.1).and.(iel.eq.1)) then
!              write(6,'(''r_basis='',3f9.4)')r_basis(1),r_basis(2),r_basis(3)
!             write(6,*)'ddorb_lagrange=',ddorb(iel,iorb)
!            endif
          enddo
        endif

!wparker Have to think about whether it works for more than one det.
        if(inum_orb.eq.6 .or. inum_orb.eq.-6) then

          if(ndet.gt.1) stop 'Smoothing splines implementation can handle only one determinant'

!wparker Use x since bwfdet_main converts to crystal lattice units for us
! Pass 1 for iw because we need the wavefunction here
! Pass 1 for igl because we need the gradient and Laplacian here
! Pass 1 for spin because we don't differentiate between spins yet

!wparker Loop over orbitals and set temporary arrays equal to permanent ones

          call bwfdet_main(x(:,iel),1,1,1,orb_blip_tmp,dorb_blip_tmp,ddorb_blip_tmp)

          do iorb=1,norb

             orb(iel,iorb)=orb_blip_tmp(iorb,1)
             dorb(:,iel,iorb)=dorb_blip_tmp(:,iorb,1)
             ddorb(iel,iorb)=ddorb_blip_tmp(iorb,1)

!wparker End of loop over orbitals
          enddo

!wparker End of if statement for blips
        endif

        if(inum_orb.eq.8 .or. inum_orb.eq.-8) then
          if(ndet.gt.1) stop 'Interpolating b-splines implementation can handle only one determinant currently'

#ifdef NOEINSPLINE
          write(6,*) 'the code must be linked to einspline library'
          stop 'the code must be linked to einspline library'
#else
          call evaluate_bsplines_with_derivatives(x(:,iel),orb_tmp,dorb_tmp,ddorb_tmp)
#endif

           do iorb=1,norb
             orb(iel,iorb)=orb_tmp(iorb)
             do k=1,ndim
                dorb(k,iel,iorb)=dorb_tmp(k,iorb)
             enddo
             ddorb(iel,iorb)=ddorb_tmp(iorb)

!wparker End of do loop over the orbitals
           enddo

!wparker End of if statement for interpolating B-splines
        endif

!wparker End of loop over the electrons
      enddo

!wp   write(6,*)'End of orbitals_period_num'

      return
      end
!-----------------------------------------------------------------------
      subroutine orbitals_period_nume(x,orb)
! Written by Cyrus Umrigar
! Calculate pw orbitals for electron iel by interpolating on a grid.
      use bwfdet_mod
      use bsplines_mod
      use orbital_grid_mod
      use coefs_mod
      use dets_mod
      use dim_mod
      use contr2_mod
      use periodic_mod
      use periodic2_mod
      implicit real*8(a-h,o-z)

      dimension x(3),orb(*),r_basis(3)

      dimension orb_blip_tmp(norb,ndet),dorb_blip_tmp(3,norb,ndet),ddorb_blip_tmp(norb,ndet)

! Determine position in lattice coordinates
! Find vector in basis coordinates
! Note we add 1 to r_basis because interpol_orb expects a positive input
! r_basis in the line after 10 is between  0 and 1 if r_basis in line 10 is > 0
!                             but between -1 and 0 if r_basis in line 10 is < 0
! r_basis in line 20 is always between 0 and 1
        isgn=1
        do 20 k=1,ndim
          r_basis(k)=0
          do 10 i=1,ndim
   10       r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*x(i)
          if(rkvec_shift_latt(k).ne.0.d0) then
            if(r_basis(k).ge.0.d0) then
              isgn=isgn*(-1)**int(r_basis(k))
             else
              isgn=isgn*(-1)**(int(r_basis(k))+1)
            endif
          endif
          r_basis(k)=r_basis(k)-int(r_basis(k))
   20     if(r_basis(k).lt.0.d0) r_basis(k)=r_basis(k)+1

      xi=r_basis(1)*ngrid_orbx
      yi=r_basis(2)*ngrid_orby
      zi=r_basis(3)*ngrid_orbz

!     write(6,'(''r_basis'',9f9.4)') r_basis,xi,yi,zi

      if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
         call interpol_orbe(ngrid_orbx,ngrid_orby,ngrid_orbz,xi,yi,zi,orb)

         do 30 iorb=1,norb
   30       orb(iorb)=orb(iorb)*isgn

      endif

      if(inum_orb.eq.6 .or. inum_orb.eq.-6) then

         if(ndet.gt.1) stop 'can only handle one determinant'

         call bwfdet_main(x,1,0,1,orb_blip_tmp,dorb_blip_tmp,
     &                    ddorb_blip_tmp)

         do iorb=1,norb
            orb(iorb)=orb_blip_tmp(iorb,1)
         enddo

      endif

! Interpolating B-splines
      if(inum_orb.eq.8 .or. inum_orb.eq.-8) then
#ifdef NOEINSPLINE
         write(6,*) 'the code must be linked to einspline library'
         stop 'the code must be linked to einspline library'
#else
         call evaluate_bsplines_function_only(x,orb)
#endif
      endif

      return
      end

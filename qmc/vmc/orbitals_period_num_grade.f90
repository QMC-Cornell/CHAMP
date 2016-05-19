      subroutine orbitals_period_num_grade(x,orb,dorb,ddorb)
! Written by Cyrus Umrigar. Modified by William Parker to interface to Einspline.
! Calculate pw orbitals, gradient and laplacian for electron at x
! by interpolating on a grid.
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

      dimension x(3),r_basis(3),orb(*),dorb(3,*),ddorb(*)

      dimension orb_blip_tmp(norb,ndet),dorb_blip_tmp(3,norb,ndet),ddorb_blip_tmp(norb,ndet)

! Determine position in lattice coordinates
! Find vector in basis coordinates
! Note we add 1 to r_basis because interpol_orb expects a positive input


!wparker This was what was originally here
!     do 20 k=1,ndim
!       r_basis(k)=0
!       do 10 i=1,ndim
!  10     r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*x(i)
!  20   r_basis(k)=r_basis(k)-nint(r_basis(k))+1

!wparker Here is what came from the later version of this subroutine
!r_basis in the line after 10 is between  0 and 1 if r_basis in line 10 is > 0
!                            but between -1 and 0 if r_basis in line 10 is < 0
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
!wparker Get the values from Lagrange interpolation routine

      if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
         call interpol_orb(ngrid_orbx,ngrid_orby,ngrid_orbz,xi,yi,zi, &
     &                     orb,dorb,ddorb)

!wparker Added from the later version of this subroutine
               do 40 iorb=1,norb
                  orb(iorb)=orb(iorb)*isgn
                  do 30 k=1,ndim
   30                dorb(k,iorb)=dorb(k,iorb)*isgn
   40             ddorb(iorb)=ddorb(iorb)*isgn

      endif

      if(inum_orb.eq.6 .or. inum_orb.eq.-6) then
         if(ndet.gt.1) stop 'can only handle one determinant with blips'
         call bwfdet_main(x,1,1,1,orb_blip_tmp,dorb_blip_tmp, &
     &                  ddorb_blip_tmp)
         do iorb=1,norb
              orb(iorb)=orb_blip_tmp(iorb,1)
              dorb(:,iorb)=dorb_blip_tmp(:,iorb,1)
              ddorb(iorb)=ddorb_blip_tmp(iorb,1)
         enddo
      endif

      if(inum_orb.eq.8 .or. inum_orb.eq.-8) then
#ifdef NOEINSPLINE
          write(6,*) 'the code must be linked to einspline library'
          stop 'the code must be linked to einspline library'
#else
         call evaluate_bsplines_with_derivatives(x,orb,dorb,ddorb)
#endif

      endif

      return
      end

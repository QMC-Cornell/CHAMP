module jastrow_mod

  use all_tools_mod
  use electrons_mod
  use basis_mod
  use orbitals_mod
  use determinants_mod
  use montecarlo_mod

! Declaration of global variables and default values
  real(dp)                               :: fen
  real(dp), allocatable                  :: dfen_drn (:)
  real(dp)                               :: feen
  real(dp), allocatable                  :: dfeen_drn (:)
  real(dp), allocatable                  :: dist_en_scaled_wlk (:,:,:)
  real(dp), allocatable                  :: dist_en_scaled_deriv1_wlk (:,:,:)
  real(dp), allocatable                  :: dist_ee_scaled2_wlk (:,:,:)
  real(dp), allocatable                  :: dist_en_scaled2_wlk (:,:,:)
  real(dp), allocatable                  :: dist_en_scaled2_deriv1_wlk (:,:,:)
  real(dp), allocatable                  :: dist_ee_scaled2_deriv1_wlk (:,:,:)

  contains

! ==============================================================================
  function dist_scaled (dist, kappa)
! ------------------------------------------------------------------------------
! Description   : scaling function used for distances in e-n and e-e jastrow
!
! Created       : J. Toulouse, 29 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  real(dp), intent(in) :: dist
  real(dp), intent(in) :: kappa

! output
  real(dp) dist_scaled

! local
  character(len=max_string_len_rout), save :: lhere = 'dist_scaled'

! initialize
  dist_scaled = 0.d0

  select case(isc)
   case (2)
     if (kappa == 0.d0) then
      dist_scaled = dist
     else
      dist_scaled = (1.d0 - dexp (-kappa * dist)) / kappa
     endif
     return
   case (4)
     dist_scaled = dist / (1.d0 + kappa * dist)
     return
   case default
     call die (lhere, 'unknown case isc='+isc+'.')
  end select

  end function dist_scaled

! ==============================================================================
  function dist_scaled_deriv1 (dist, kappa)
! ------------------------------------------------------------------------------
! Description   : first derivative of dist_scaled
!
! Created       : J. Toulouse, 29 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  real(dp), intent(in) :: dist
  real(dp), intent(in) :: kappa

! output
  real(dp) dist_scaled_deriv1

! local
  character(len=max_string_len_rout), save :: lhere = 'dist_scaled_deriv1'

! initialize
  dist_scaled_deriv1 = 0.d0

  select case(isc)
   case (2)
     dist_scaled_deriv1 = dexp (-kappa * dist)
     return
   case (4)
     dist_scaled_deriv1 = 1.d0/(1.d0 + kappa * dist)**2
     return
   case default
     call die (lhere, 'unknown case isc='+isc+'.')
  end select

  end function dist_scaled_deriv1

! ==============================================================================
  function dist_scaled2 (dist, kappa)
! ------------------------------------------------------------------------------
! Description   : scaling function used for distances in e-e-n jastrow
!
! Created       : J. Toulouse, 30 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  real(dp), intent(in) :: dist
  real(dp), intent(in) :: kappa

! output
  real(dp) dist_scaled2

! local
  character(len=max_string_len_rout), save :: lhere = 'dist_scaled2'

! initialize
  dist_scaled2 = 0.d0

  select case(isc)
   case (2)
     dist_scaled2 = dexp (-kappa * dist)
     return
   case (4)
     dist_scaled2 = 1.d0 / (1.d0 + kappa * dist)
     return
   case default
     call die (lhere, 'unknown case isc='+isc+'.')
  end select

  end function dist_scaled2

! ==============================================================================
  function dist_scaled2_deriv1 (dist, kappa)
! ------------------------------------------------------------------------------
! Description   : first derivative of dist_scaled2
!
! Created       : J. Toulouse, 30 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  real(dp), intent(in) :: dist
  real(dp), intent(in) :: kappa

! output
  real(dp) dist_scaled2_deriv1

! local
  character(len=max_string_len_rout), save :: lhere = 'dist_scaled2_deriv1'

! initialize
  dist_scaled2_deriv1 = 0.d0

  select case(isc)
   case (2)
     dist_scaled2_deriv1 = -kappa * dexp (-kappa * dist)
     return
   case (4)
     dist_scaled2_deriv1 = -kappa / (1.d0 + kappa * dist)**2
     return
   case default
     call die (lhere, 'unknown case isc='+isc+'.')
  end select

  end function dist_scaled2_deriv1

! ==============================================================================
  subroutine dist_en_scaled_wlk_bld
! ------------------------------------------------------------------------------
! Description   : scaled electron-nuclei distances for e-n jastrow
!
! Created       : J. Toulouse, 29 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer  :: elec_i, cent_i, walk_i

! header
  if (header_exe) then

   call object_create ('dist_en_scaled_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('ncent')
   call object_needed ('dist_en_wlk')
   call object_needed ('scalek')

   return

  endif

! begin

! allocations
  call object_alloc ('dist_en_scaled_wlk', dist_en_scaled_wlk, nelec, ncent, nwalk)

  do walk_i = 1, nwalk
    do cent_i = 1, ncent
      do elec_i = 1, nelec
         dist_en_scaled_wlk (elec_i, cent_i, walk_i) = dist_scaled (dist_en_wlk (elec_i, cent_i, walk_i), scalek(1))
       enddo ! elec_i
    enddo ! cent_i
  enddo ! walk_i

  end subroutine dist_en_scaled_wlk_bld

! ==============================================================================
  subroutine dist_en_scaled_deriv1_wlk_bld
! ------------------------------------------------------------------------------
! Description   : first derivative of scaled electron-nuclei distances for e-n jastrow
!
! Created       : J. Toulouse, 29 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer  :: elec_i, cent_i, walk_i

! header
  if (header_exe) then

   call object_create ('dist_en_scaled_deriv1_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('ncent')
   call object_needed ('dist_en_wlk')
   call object_needed ('scalek')

   return

  endif

! begin

! allocations
  call object_alloc ('dist_en_scaled_deriv1_wlk', dist_en_scaled_deriv1_wlk, nelec, ncent, nwalk)

  do walk_i = 1, nwalk
    do cent_i = 1, ncent
      do elec_i = 1, nelec
         dist_en_scaled_deriv1_wlk (elec_i, cent_i, walk_i) = dist_scaled_deriv1 (dist_en_wlk (elec_i, cent_i, walk_i), scalek(1))
       enddo ! elec_i
    enddo ! cent_i
  enddo ! walk_i

  end subroutine dist_en_scaled_deriv1_wlk_bld

! ==============================================================================
  subroutine dist_en_scaled2_wlk_bld
! ------------------------------------------------------------------------------
! Description   : scaled electron-nuclei distances for e-e-n jastrow
!
! Created       : J. Toulouse, 30 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer  :: elec_i, cent_i, walk_i

! header
  if (header_exe) then

   call object_create ('dist_en_scaled2_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('ncent')
   call object_needed ('dist_en_wlk')
   call object_needed ('scalek')

   return

  endif

! begin

! allocations
  call object_alloc ('dist_en_scaled2_wlk', dist_en_scaled2_wlk, nelec, ncent, nwalk)

  do walk_i = 1, nwalk
    do cent_i = 1, ncent
      do elec_i = 1, nelec
         dist_en_scaled2_wlk (elec_i, cent_i, walk_i) = dist_scaled2 (dist_en_wlk (elec_i, cent_i, walk_i), scalek(1))
       enddo ! elec_i
    enddo ! cent_i
  enddo ! walk_i

  end subroutine dist_en_scaled2_wlk_bld

! ==============================================================================
  subroutine dist_ee_scaled2_wlk_bld
! ------------------------------------------------------------------------------
! Description   : scaled electron-electron distances for e-e-n jastrow
!
! Created       : J. Toulouse, 30 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer  :: elec_i, elec_j, walk_i

! header
  if (header_exe) then

   call object_create ('dist_ee_scaled2_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('scalek')

   return

  endif

! begin

! allocations
  call object_alloc ('dist_ee_scaled2_wlk', dist_ee_scaled2_wlk, nelec, nelec, nwalk)

  do walk_i = 1, nwalk
    do elec_j = 1, nelec
      do elec_i = elec_j+1, nelec
         dist_ee_scaled2_wlk (elec_i, elec_j, walk_i) = dist_scaled2 (dist_ee_wlk (elec_i, elec_j, walk_i), scalek(1))
         dist_ee_scaled2_wlk (elec_j, elec_i, walk_i) = dist_ee_scaled2_wlk (elec_i, elec_j, walk_i)
       enddo ! elec_i
    enddo ! cent_j
  enddo ! walk_i

  end subroutine dist_ee_scaled2_wlk_bld

! ==============================================================================
  subroutine dist_en_scaled2_deriv1_wlk_bld
! ------------------------------------------------------------------------------
! Description   : first derivative of scaled electron-nuclei distances for e-e-n jastrow
!
! Created       : J. Toulouse, 30 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer  :: elec_i, cent_i, walk_i

! header
  if (header_exe) then

   call object_create ('dist_en_scaled2_deriv1_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('ncent')
   call object_needed ('dist_en_wlk')
   call object_needed ('scalek')

   return

  endif

! begin

! allocations
  call object_alloc ('dist_en_scaled2_deriv1_wlk', dist_en_scaled2_deriv1_wlk, nelec, ncent, nwalk)

  do walk_i = 1, nwalk
    do cent_i = 1, ncent
      do elec_i = 1, nelec
         dist_en_scaled2_deriv1_wlk (elec_i, cent_i, walk_i) = dist_scaled2_deriv1 (dist_en_wlk (elec_i, cent_i, walk_i), scalek(1))
       enddo ! elec_i
    enddo ! cent_i
  enddo ! walk_i

  end subroutine dist_en_scaled2_deriv1_wlk_bld

! ==============================================================================
  subroutine dist_ee_scaled2_deriv1_wlk_bld
! ------------------------------------------------------------------------------
! Description   : derivatives of scaled electron-electron distances for e-e-n jastrow
!
! Created       : J. Toulouse, 31 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer  :: elec_i, elec_j, walk_i

! header
  if (header_exe) then

   call object_create ('dist_ee_scaled2_deriv1_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('scalek')

   return

  endif

! begin

! allocations
  call object_alloc ('dist_ee_scaled2_deriv1_wlk', dist_ee_scaled2_deriv1_wlk, nelec, nelec, nwalk)

  do walk_i = 1, nwalk
    do elec_j = 1, nelec
      do elec_i = elec_j+1, nelec
         dist_ee_scaled2_deriv1_wlk (elec_i, elec_j, walk_i) = dist_scaled2_deriv1 (dist_ee_wlk (elec_i, elec_j, walk_i), scalek(1))
         dist_ee_scaled2_deriv1_wlk (elec_j, elec_i, walk_i) = dist_ee_scaled2_deriv1_wlk (elec_i, elec_j, walk_i)
       enddo ! elec_i
    enddo ! cent_j
  enddo ! walk_i

  end subroutine dist_ee_scaled2_deriv1_wlk_bld

! ==============================================================================
  subroutine fen_bld
! ------------------------------------------------------------------------------
! Description   : Jastrow function fen for ijas=4 (used only for checkings)
!
! Created       : J. Toulouse, 30 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer cent_i, elec_i, cent_type_i, order_i

! header
  if (header_exe) then

   call object_create ('fen')

   call object_needed ('ncent')
   call object_needed ('iwctype')
   call object_needed ('nelec')
   call object_needed ('a4')
   call object_needed ('norda')
   call object_needed ('dist_en_scaled_wlk')

   return

  endif

! begin

! association
  call object_associate ('fen', fen)

! implemented only for ijas=4
  if (ijas /= 4) then
    call die (here, 'implemented only for ijas=4.')
  endif

  fen = 0.d0
  do cent_i = 1, ncent
    cent_type_i = iwctype (cent_i)
    do elec_i = 1, nelec
      fen = fen + a4(1,cent_type_i,1)*dist_en_scaled_wlk (elec_i, cent_i, 1)/(1.d0 + a4(2,cent_type_i,1)*dist_en_scaled_wlk (elec_i, cent_i, 1))
      do order_i = 2, norda
        fen = fen + a4(order_i+1,cent_type_i,1) * (dist_en_scaled_wlk (elec_i, cent_i, 1)**(order_i))
      enddo ! order_i
        fen = fen - asymp_jasa(cent_type_i,1)
    enddo ! elec_i
  enddo ! cent_i

  end subroutine fen_bld

! ==============================================================================
  subroutine dfen_drn_bld
! ------------------------------------------------------------------------------
! Description   : derivative of Jastrow function fen wrt nuclear coordinates
! Description   : for ijas=4
!
! Created       : J. Toulouse, 27 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer force_i, cent_i, dim_i, elec_i, cent_type_i, order_i
  real(dp) :: dfen_drn_scaled

! header
  if (header_exe) then

   call object_create ('dfen_drn')

   call object_needed ('forces_nb')
   call object_needed ('forces_cent')
   call object_needed ('forces_direct')
   call object_needed ('iwctype')
   call object_needed ('nelec')
   call object_needed ('a4')
   call object_needed ('norda')
   call object_needed ('dist_en_scaled_wlk')
   call object_needed ('dist_en_scaled_deriv1_wlk')
   call object_needed ('grd_dist_en')

   return

  endif

! begin

! allocation
  call object_alloc ('dfen_drn', dfen_drn, forces_nb)

! implemented only for ijas=4
  if (ijas /= 4) then
    call die (here, 'implemented only for ijas=4.')
  endif

  do force_i = 1, forces_nb
    dfen_drn (force_i) = 0.d0
    cent_i = forces_cent (force_i)
    cent_type_i = iwctype (cent_i)
    dim_i = forces_direct (force_i)
    do elec_i = 1, nelec
      dfen_drn_scaled = a4(1,cent_type_i,1)/(1.d0 + a4(2,cent_type_i,1)*dist_en_scaled_wlk (elec_i, cent_i, 1))**2
      do order_i = 2, norda
        dfen_drn_scaled = dfen_drn_scaled + a4(order_i+1,cent_type_i,1) * order_i * (dist_en_scaled_wlk (elec_i, cent_i, 1)**(order_i -1))
      enddo ! order_i
      dfen_drn (force_i) = dfen_drn (force_i) - dfen_drn_scaled * dist_en_scaled_deriv1_wlk (elec_i, cent_i, 1) * grd_dist_en (dim_i, elec_i, cent_i)
    enddo ! elec_i
  enddo ! force_i

  end subroutine dfen_drn_bld

! ==============================================================================
  subroutine feen_bld
! ------------------------------------------------------------------------------
! Description   : Jastrow function feen for ijas=4 (used only for checkings)
!
! Created       : J. Toulouse, 30 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer elec_i, elec_j, cent_i, cent_type_i
  integer ll, n, k, m, l, l_hi
  real(dp) uu, rri, rrj, tt

! header
  if (header_exe) then

   call object_create ('feen')

   call object_needed ('nelec')
   call object_needed ('ncent')
   call object_needed ('iwctype')
   call object_needed ('c')
   call object_needed ('nordc')
   call object_needed ('dist_ee_scaled2_wlk')
   call object_needed ('dist_en_scaled2_wlk')

   return

  endif

! begin

! association
  call object_associate ('feen', feen)

! implemented only for ijas=4
  if (ijas /= 4) then
    call die (here, 'implemented only for ijas=4.')
  endif

  feen = 0.d0
  do elec_i = 2, nelec
    do elec_j = 1, elec_i - 1

!     scaled ee distance
      uu = dist_ee_scaled2_wlk (elec_i, elec_j, 1)

      do cent_i = 1, ncent

!       scaled en distances
        rri = dist_en_scaled2_wlk (elec_i, cent_i, 1)
        rrj = dist_en_scaled2_wlk (elec_j, cent_i, 1)
        tt = rri * rrj

        cent_type_i = iwctype (cent_i)
        ll=0
        do n=2,nordc
          do k=n-1,0,-1
            if(k.eq.0) then
              l_hi=n-k-2
             else
              l_hi=n-k
            endif
            do l=l_hi,0,-1
              m=(n-k-l)/2
              if(2*m.eq.n-k-l) then
                ll=ll+1
                feen = feen + c(ll,cent_type_i,1) * (uu**k) * (rri**l+rrj**l) * (tt**m)
              endif
            enddo ! l
          enddo ! k
        enddo ! n
      enddo ! cent_i
    enddo ! elec_j
  enddo ! elec_i

  end subroutine feen_bld

! ==============================================================================
  subroutine dfeen_drn_bld
! ------------------------------------------------------------------------------
! Description   : derivative of Jastrow function feen wrt nuclear coordinates
! Description   : for ijas=4
! Description   : coding to be optimized?
!
! Created       : J. Toulouse, 31 Jul 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer elec_i, elec_j, cent_i, cent_type_i, force_i, dim_i
  integer ll, n, k, m, l, l_hi
  real(dp) uu, rri, rrj, tt, drri, drrj, dri_drn,  drj_drn

! header
  if (header_exe) then

   call object_create ('dfeen_drn')

   call object_needed ('forces_nb')
   call object_needed ('forces_cent')
   call object_needed ('forces_direct')
   call object_needed ('nelec')
   call object_needed ('ncent')
   call object_needed ('iwctype')
   call object_needed ('c')
   call object_needed ('nordc')
   call object_needed ('dist_ee_scaled2_wlk')
   call object_needed ('dist_en_scaled2_wlk')
   call object_needed ('dist_en_scaled2_deriv1_wlk')
   call object_needed ('grd_dist_en')

   return

  endif

! begin

! association
  call object_alloc ('dfeen_drn', dfeen_drn, forces_nb)

! implemented only for ijas=4
  if (ijas /= 4) then
    call die (here, 'implemented only for ijas=4.')
  endif

  do force_i = 1, forces_nb
    dfeen_drn (force_i) = 0.d0
    cent_i = forces_cent (force_i)
    cent_type_i = iwctype (cent_i)
    dim_i = forces_direct (force_i)

    do elec_i = 2, nelec
    
!     scaled en distance and derivatives
      rri = dist_en_scaled2_wlk (elec_i, cent_i, 1)
      drri = dist_en_scaled2_deriv1_wlk (elec_i, cent_i, 1)
      dri_drn = -grd_dist_en (dim_i, elec_i, cent_i)
    
      do elec_j = 1, elec_i - 1
    
!       scaled ee distance
        uu = dist_ee_scaled2_wlk (elec_i, elec_j, 1)
    
!       scaled en distance
        rrj = dist_en_scaled2_wlk (elec_j, cent_i, 1)
        drrj = dist_en_scaled2_deriv1_wlk (elec_j, cent_i, 1)
        drj_drn = -grd_dist_en (dim_i, elec_j, cent_i)
        tt = rri * rrj
    
        ll=0
        do n=2,nordc
          do k=n-1,0,-1
            if(k.eq.0) then
              l_hi=n-k-2
             else
              l_hi=n-k
            endif
            do l=l_hi,0,-1
              m=(n-k-l)/2
              if(2*m.eq.n-k-l) then
                ll=ll+1
                dfeen_drn (force_i) = dfeen_drn (force_i) + c(ll,cent_type_i,1) * (uu**k) * (                        &
                   ( (l*rri**(l-1)) * (tt**m) + (rri**l+rrj**l) * m*(rri**(m-1)) * (rrj**m) ) * drri * dri_drn       &
                 + ( (l*rrj**(l-1)) * (tt**m) + (rri**l+rrj**l) * m*(rrj**(m-1)) * (rri**m) ) * drrj * drj_drn  ) 
              endif
            enddo ! l
          enddo ! k
        enddo ! n
      enddo ! elec_j
    enddo ! elec_i
  enddo ! force_i

  end subroutine dfeen_drn_bld

end module jastrow_mod

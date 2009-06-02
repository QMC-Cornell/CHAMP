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

!===========================================================================
  subroutine jastrow_menu
!---------------------------------------------------------------------------
! Description : menu for jastrow
!
! Created     : J. Toulouse, 08 Apr 2009
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'jastrow_menu'
  integer isp, iparm, it
  real(dp) parm2min, cutjas_ee_tmp, cutjas_en_tmp
  real(dp),parameter :: eps=1.d-4


! begin
  write(6,*)
  write(6,'(a)') 'Beginning of jastrow menu --------------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for jastrow menu:'
   write(6,'(a)') 'basis'
   write(6,'(a)') '  ijas = [integer] type of Jastrow factor (default=4)'
   write(6,'(a)') '  isc = [integer] type of scaling function for coordinates (default=2)'
   write(6,'(a)') '  nspin1 = [integer] starting spin index (default=1)'
   write(6,'(a)') '  nspin2 = [integer] ending spin index (default=1)'
   write(6,'(a)') '  nord = [integer] order of the polynomial (default=5)'
   write(6,'(a)') '  fock = [integer] Fock terms (default=0)'
   write(6,'(a)') '  ianalyt_lap = [integer] analytic Laplacian of Jastrow?  (default=1)'
   write(6,'(a)') '  scalek = [real] scale factor for ijas>= 2 isc>= 2 (default=0.5)'
   write(6,'(a)') '  a21 = [real] some constant for ijas and isc>= 2 (default=0.)'
   write(6,'(a)') '  norda = [integer] order of e-n polynomial for ijas >= 4 (default=5)'
   write(6,'(a)') '  nordb = [integer] order of e-e polynomial for ijas >= 4 (default=5)'
   write(6,'(a)') '  nordc = [integer] order of e-e-n polynomial for ijas >= 4 (default=5)'
   write(6,'(a)') '  parameters ... end = e-n, e-e and e-e-n parameters'
   write(6,'(a)') 'end'
   write(6,*)


  case ('ijas')
   call get_next_value (ijas)

  case ('isc')
   call get_next_value (isc)
   call object_modified ('isc')

  case ('nspin1')
   call get_next_value (nspin1)

  case ('nspin2')
   call get_next_value (nspin2)

  case ('nord')
   call get_next_value (nord)

  case ('ifock')
   call get_next_value (ifock)

  case ('ianalyt_lap')
   call get_next_value (ianalyt_lap)

  case ('scalek')
   call object_provide ('nwf')
   call alloc ('scalek', scalek, nwf)
   call get_next_value (scalek(1))
   call object_modified ('scalek')

  case ('a21')
   call get_next_value (a21)

  case ('norda')
   call get_next_value (norda)
   call object_modified ('norda')

  case ('nordb')
   call get_next_value (nordb)
   call object_modified ('nordb')

  case ('nordc')
   call get_next_value (nordc)
   call object_modified ('nordc')

  case ('cutjas_en')
   call get_next_value (cutjas_en_tmp)

  case ('cutjas_ee')
   call get_next_value (cutjas_ee_tmp)

  case ('parameters')
   call jastrow_parameters_rd

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  write(6,'(a,i4)') ' type of Jastrow factor: ijas=',ijas
  write(6,'(a,i4)') ' type of scaling function: isc=',isc
  write(6,'(a,i4)') ' starting spin index: nspin1=',nspin1
  write(6,'(a,i4)') ' ending   spin index: nspin2=',nspin2
  write(6,'(a,i4)') ' order of polynomial: nord=',nord
  write(6,'(a,i4)') ' Fock terms: ifock=',ifock
  write(6,'(a,i4)') ' analytic Laplacian: ianalyt_lap=',ianalyt_lap

  call object_provide ('nloc')
  call object_provide ('ndn')

  if(ianalyt_lap.eq.0 .and. nloc.gt.0) stop 'Cannot have numerical Lap. with pseudopot'
  if(ianalyt_lap.eq.0 .and. iperiodic.gt.0) stop 'Cannot have numerical Lap. with periodic system: distances in jastrow_num not correct'
  if(ijas.ne.4 .and. iperiodic.gt.0) stop 'Only ijas=4 implemented for periodic systems'
  if(ijas.gt.6) stop 'only ijas=1,2,3,4,5,6 implemented'
  if(ifock.lt.0.or.ifock.gt.4) stop 'ifock must be between 0 and 4'
  if(ndn.eq.1.and.nspin2.eq.3) stop '1 spin down and nspin2=3'
  if((ijas.eq.4.or.ijas.eq.5).and.(isc.ne.2.and.isc.ne.4.and.isc.ne.6.and.isc.ne.7.and.    &
      isc.ne.8.and.isc.ne.10.and.isc.ne.12.and.isc.ne.14.and.isc.ne.16.and.isc.ne.17))     &
       stop 'if ijas=4 or 5, isc must be one of 2,4,6,7,8,10,12,14,16,17'
  if((ijas.eq.6).and.(isc.ne.6.and.isc.ne.7)) stop 'if ijas=6, isc must be 6 or 7'

  if(ijas.eq.3.and.nspin2.gt.1) stop 'ijas=3 and nspin2>1'


  if(ijas.eq.1) then
    write(6,'(a,f10.5)') ' Jastrow numerator =',cjas1(1)
    write(6,'(a,f10.5)') ' Jastrow denominator =',cjas2(1)
  elseif(ijas.eq.2) then
!    nparm_read=69
    write(6,'(a,f10.5)') ' scale factor: scalek=',scalek(1)
    write(6,'(a,f10.5)') ' a21=',a21
    do isp=nspin1,nspin2
       call object_provide ('ncent')
       if(ncent.gt.1.and.a1(2,isp,1).ne.zero) then
         write(6,'(a)') ' Warning: e-n cusp condition cannot be imposed for molecules with present weighted form of Jastrow'
       endif
       write(6,'(a,10f10.6)') ' e-n terms: a=',(a1(iparm,isp,1),iparm=1,nparm_read)
     enddo       
     do isp=nspin1,nspin2
       write(6,'(a,10f10.6)') ' e-e terms: b=',(a2(iparm,isp,1),iparm=1,nparm_read)
     enddo

   elseif(ijas.eq.3) then
!     nparm_read=2
!     nparmc_read=(nord**3+5*nord)/6+nord**2+nord
!     write(6,'(a,3i5)') ' nparm_read, nparmc_read=', nparm_read,nparmc_read
     if(isc.ge.2) then
       write(6,'(a,f10.5)') ' scale factor: scalek=',scalek(1)
       write(6,'(a,f10.5)') ' a21=',a21
     endif
     write(6,'(a,10f10.6)') ' e-n terms: a=',(a(iparm,1),iparm=1,nparm_read)
     do isp=nspin1,nspin2b
        write(6,'(a,10f10.6)') ' e-e terms: b=',(b(iparm,isp,1),iparm=1,nparm_read)
     enddo
     do it=1,nctype
        write(6,'(a,50f10.6)') ' e-e-n terms: c=',(c(iparm,it,1),iparm=1,nparmc_read)
     enddo
     if(ifock.gt.0) then
          do it=1,nctype
            write(6,'(a,10f10.6)') ' Fock terms f=',(fck(iparm,it,1),iparm=1,nfock)
          enddo
     endif
   elseif(ijas.ge.4.and.ijas.le.6) then
     if(ifock.gt.0) stop 'fock not yet implemented for ijas=4,5,6'
        write(6,'(a,i5)') ' order of e-n polynomial: norda=',norda
        write(6,'(a,i5)') ' order of e-e polynomial: nordb=',nordb
        write(6,'(a,i5)') ' order of e-e-n polynomial: nordc=',nordc
!        nparma_read=2+max(0,norda-1)
!        nparmb_read=2+max(0,nordb-1)
!        nparmc_read=nterms4(nordc)
!        write(6,'(a,3i5)') ' nparma_read,nparmb_read,nparmc_read=', nparma_read,nparmb_read,nparmc_read
        if(iperiodic.gt.0 .and. nordc.gt.0 .and. ijas .le. 3) stop 'J_een only implemented with ijas= 4,5,6'
        if(isc.ge.2) then
          write(6,'(a,f10.5)') ' scale factor: scalek=',scalek(1)
          write(6,'(a,f10.5)') ' a21=',a21
        endif
        if(isc.ne.8 .and. isc.ne.10) then
          parm2min=-scalek(1)
        else
          parm2min=-1.d0
        endif
        do it=1,nctype
           write(6,'(a,10f10.6)') ' e-n terms: a=',(a4(iparm,it,1),iparm=1,nparma_read)
           if(nparma_read.ge.2 .and. a4(2,it,1).lt.parm2min) then
               write(6,'(a)') ' Warning: a4(2,it,1) too low, Jastrow denom could become negative'
               stop 'a4(2,it,1) too low, Jastrow denom could become negative'
             else
           endif
        enddo     
        do isp=nspin1,nspin2b
          write(6,'(a,10f10.6)') ' e-e terms: b=',(b(iparm,isp,1),iparm=1,nparmb_read)
           if(nparmb_read.ge.2 .and. b(2,isp,1).lt.parm2min) then
             write(6,'(a)') ' Warning: b(2,isp,1) too low, Jastrow denom could become negative'
             stop 'b(2,isp,1) too low, Jastrow denom could become negative'
           endif
        enddo    
        do it=1,nctype
          write(6,'(a,50f10.6)') ' e-e-n terms: c=',(c(iparm,it,1),iparm=1,nparmc_read)
       enddo
! Note: Fock terms yet to be put in ijas=4,5,6.
   endif

! Read cutoff for Jastrow4,5,6 and call set_scale_dist to evaluate constants
! that need to be reset if scalek is being varied.
! If cutjas=0, then reset cutjas_en, cutjas_ee to infinity
! Warning: At present we are assuming that the same scalek is used
! for primary and secondary wavefns.  Otherwise c1_jas6i,c1_jas6,c2_jas6
! should be dimensioned to MWF
      if(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17) then
!        read(5,*) cutjas_en_tmp,cutjas_ee_tmp
        if(iperiodic.ne.0 .and. cutjas_en_tmp.gt.cutjas_en+eps) then
          write(6,'(''Warning: input cutjas > half shortest primitive cell lattice vector;cutjas_en reset from'',f9.5,'' to'',f9.5)') cutjas_en_tmp,cutjas_en
         else
          if(cutjas_en_tmp.lt.cutjas_en-eps) then
            write(6,'(''Warning: Could use larger cutjas_en='',f9.5,'' instead of the input value='',f9.5)') cutjas_en,cutjas_en_tmp
          endif
          write(6,'(''input cutjas_en='',d12.5)') cutjas_en_tmp
          cutjas_en=cutjas_en_tmp
        endif
        if(iperiodic.ne.0 .and. cutjas_ee_tmp.gt.cutjas_ee+eps) then
          write(6,'(''Warning: input cutjas > half shortest simulation cell lattice vector;cutjas_ee reset from'',f9.5,'' to'',f9.5)') cutjas_ee_tmp,cutjas_ee
         else
          if(cutjas_ee_tmp.lt.cutjas_ee-eps) then
            write(6,'(''Warning: Could use larger cutjas_ee='',f9.5,'' instead of the input value='',f9.5)') cutjas_ee,cutjas_ee_tmp
          endif
          write(6,'(''input cutjas_ee='',d12.5)') cutjas_ee_tmp
          cutjas_ee=cutjas_ee_tmp
        endif
        if(cutjas_en_tmp.le.0.d0) then
          write(6,'(''cutjas_en reset to infinity'')')
          cutjas_en=1.d99
        endif
        if(cutjas_ee_tmp.le.0.d0) then
          write(6,'(''cutjas_ee reset to infinity'')')
          cutjas_ee=1.d99
        endif
      endif
      call set_scale_dist(-1,1)

      if(ifock.gt.0) then
!     Setup for Chris' Fock
!         fflag=7

! Read pars for Chris's wf
!       call wfpars
        if(ifock.eq.4) then
        call die (lhere, 'fock terms need to be updated')
!          open(11, file =
!     &    '/afs/theory.cornell.edu/user/tc/cyrus/qmc/vmc/lob.dat')
!          rewind 11
!          read(11,*) (rlobx(i),rloby(i),i=1,nsplin)
!          call spline(rlobx,rloby,nsplin,0.d0,0.d0,rloby2)
        endif
      endif


  write(6,'(a)') 'End of jastrow menu --------------------------------------------------------------------------------------'

  end subroutine jastrow_menu

!===========================================================================
  subroutine jastrow_parameters_rd
!---------------------------------------------------------------------------
! Description : read Jastrow parameters (e-n, e-e, e-e-n)
!
! Created     : J. Toulouse, 09 Apr 2009
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'jastrow_parameters_rd'
  integer isp, iparm, it

  integer, external :: nterms4

! begin
  call object_provide ('nctype')

  nspin2b=iabs(nspin2)
  nocuspb=0
  if(nspin2.lt.0) then
    if(nspin2.eq.-1) nocuspb=1
    nspin2=1
  endif

  if(ijas.eq.1) then
    read(5,*) cjas1(1),cjas2(1)
  elseif(ijas.eq.2) then
    nparm_read=69
    call alloc ('a1', a1, nparm_read, nspin2-nspin1+1, nwf)
    do isp=nspin1,nspin2
       read(5,*) (a1(iparm,isp,1),iparm=1,nparm_read)
    enddo
    call alloc ('a2', a2, nparm_read, nspin2-nspin1+1, nwf)
    do isp=nspin1,nspin2
       read(5,*) (a2(iparm,isp,1),iparm=1,nparm_read)
    enddo
  elseif(ijas.eq.3) then
    nparm_read=2
    nparmc_read=(nord**3+5*nord)/6+nord**2+nord
!    write(6,'(a,3i5)') ' nparm_read, nparmc_read=', nparm_read,nparmc_read
    call alloc ('a', a, nparm_read, nwf)
    read(5,*) (a(iparm,1),iparm=1,nparm_read)
    call alloc ('b', b, nparm_read, nspin2b-nspin1+1,nwf)
    do isp=nspin1,nspin2b
       read(5,*) (b(iparm,isp,1),iparm=1,nparm_read)
    enddo
    call alloc ('c', c, nparmc_read, nctype, nwf)
    do it=1,nctype
       read(5,*) (c(iparm,it,1),iparm=1,nparmc_read)
    enddo
    if(ifock.gt.0) then
      nfock=9
      if(ifock.eq.2) nfock=15
         call alloc ('fck', fck, nfock, nctype, nwf)
         do it=1,nctype
           read(5,*) (fck(iparm,it,1),iparm=1,nfock)
         enddo
    endif
  elseif(ijas.ge.4.and.ijas.le.6) then
       nparma_read=2+max(0,norda-1)
       nparmb_read=2+max(0,nordb-1)
       nparmc_read=nterms4(nordc)
!       write(6,'(a,3i5)') ' nparma_read,nparmb_read,nparmc_read=', nparma_read,nparmb_read,nparmc_read
       call alloc ('a4', a4, nparma_read, nctype, nwf)
       do it=1,nctype
          read(5,*) (a4(iparm,it,1),iparm=1,nparma_read)
       enddo     
       call object_modified ('a4')
       call alloc ('b', b, nparmb_read, nspin2b-nspin1+1,nwf)
       do isp=nspin1,nspin2b
         read(5,*) (b(iparm,isp,1),iparm=1,nparmb_read)
       enddo    
       call alloc ('c', c, nparmc_read, nctype, nwf)
       do it=1,nctype
         read(5,*) (c(iparm,it,1),iparm=1,nparmc_read)
      enddo
   endif

   call read_up_to_end

   call object_modified ('nspin2b')
   call object_modified ('nparma_read')
   call object_modified ('nparmb_read')
   call object_modified ('nparmc_read')
   call object_modified ('b')
   call object_modified ('c')


  end subroutine jastrow_parameters_rd

! ==============================================================================
  function dist_scaled (dist, kappa)
! ------------------------------------------------------------------------------
! Description   : scaling function used for distances in e-n and e-e jastrow
!
! Created       : J. Toulouse, 29 Jul 2008
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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
  include 'modules.h'
  implicit none

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

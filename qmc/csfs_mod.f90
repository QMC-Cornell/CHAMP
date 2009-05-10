module csfs_mod

  use all_tools_mod
  use orbitals_mod
  use determinants_mod

! declaration of global variables and default values
  type (type_integer_row), allocatable :: det_unq_up_in_csf (:)
  type (type_integer_row), allocatable :: det_unq_dn_in_csf (:)
  type (type_real_row),    allocatable :: cdet_unq_in_csf (:)
  real(dp), allocatable                :: csf_prefac (:)

  real(dp)                  :: wfdet_ovlp
  real(dp), allocatable     :: dens_mat_wfdet_up (:,:)
  real(dp), allocatable     :: dens_mat_wfdet_dn (:,:)
  real(dp), allocatable     :: dens_mat_wfdet (:,:)

  real(dp), allocatable     :: csfs_ovlp (:,:)
  real(dp), allocatable     :: csfs_wfdet_ovlp (:)

  contains

!===========================================================================
  subroutine csfs_menu
!---------------------------------------------------------------------------
! Description : menu for CSFs
!
! Created     : J. Toulouse, 07 Apr 2009
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'csfs_menu'
  integer elec_i, elec_j, det_i, csf_i, det_in_csf_i
  real(dp), allocatable :: csf_coef_read (:)

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of csfs menu -----------------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for csfs menu:'
   write(6,'(a)') ' csfs'
   write(6,'(a)') '  determinants : list of orbitals occupied in spin-up and spin-down determinants'
   write(6,'(a)') '   1 2  1 2'
   write(6,'(a)') '   1 3  1 2'
   write(6,'(a)') '  end'
   write(6,'(a)') '  csf_coef 1.0 0.3 ... end : CSF coefficients'
   write(6,'(a)') '  dets_in_csfs ... end : determinants in CSFs'
   write(6,'(a)') ' end'
   write(6,*)

  case ('determinants')
    call determinants_rd 

  case ('csf_coef')
   call get_next_value_list ('csf_coef_read', csf_coef_read, ncsf)
   call object_modified ('nforce')
   call alloc ('csf_coef', csf_coef, ncsf, max(3,nforce))
   csf_coef (:,1) = csf_coef_read (:)
   call object_modified ('ncsf')
   call object_modified ('csf_coef')

  case ('dets_in_csfs')
   call dets_in_csfs_rd

  case ('end')
   exit
  
  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select
   
  enddo ! end loop over menu lines

  call object_provide ('nelec')
  call object_provide ('nup')
  call object_provide ('norb')
  call object_provide ('ndet')
  write(6,'(a,i5)') ' number of determinants = ',ndet
  write(6,'(a)') ' determinants have orbitals (spin-up | spin-down):'
  do det_i=1,ndet
     write(6,'(a,i5,a,100i4)',advance='no') ' det # ', det_i, ': ', (iworbd(elec_i,det_i),elec_i=1,nup) 
     write(6,'(a,100i4)') '  |', (iworbd(elec_i,det_i),elec_i=nup+1,nelec) 
     do elec_i=1,nelec
       if(iworbd(elec_i,det_i).gt.norb) stop 'iworbd(j,i) > norb'
     enddo
  enddo

  do det_i = 1, ndet
    do elec_i=2,nup
       do elec_j=1,elec_i-1
          if(iworbd(elec_i,det_i).eq.iworbd(elec_j,det_i)) then
           call die (lhere, 'a spin-up determinant has 2 identical orbitals')
          endif
       enddo
    enddo
    do elec_i=2,ndn
       do elec_j=1,elec_i-1
          if(iworbd(nup+elec_i,det_i).eq.iworbd(nup+elec_j,det_i)) then
           call die (lhere, 'a spin-down determinant has 2 identical orbitals')
          endif
       enddo
    enddo
   enddo  

! determine unique up and dn determinants
  call determinant_up_dn

  write(6,*)
  call object_provide ('ncsf')
  write(6,'(a,i5)') ' number of CSFs = ',ncsf
  write(6,'(a,200f10.6)') ' CSF coefficients = ',(csf_coef(csf_i,1),csf_i=1,ncsf)

  call object_provide ('ndet_in_csf')
  call object_provide ('iwdet_in_csf')
  call object_provide ('cdet_in_csf')
  do csf_i = 1, ncsf
    write(6,'(a,i3)') ' CSF # ',csf_i
    write(6,'(a,200i4)') ' determinants in CSF:',(iwdet_in_csf(det_in_csf_i,csf_i),det_in_csf_i=1,ndet_in_csf(csf_i))
    write(6,'(a,200f8.5)') ' coefficients:',(cdet_in_csf(det_in_csf_i,csf_i),det_in_csf_i=1,ndet_in_csf(csf_i))
  enddo

! number of orbitals to compute
  call object_provide ('orb_occ_last_in_wf_lab')
  norb = orb_occ_last_in_wf_lab
  write(6,'(a,i8)') ' Number of computed orbitals initialized to ', norb
  call object_modified ('norb')

  write(6,'(a)') 'End of csfs menu -----------------------------------------------------------------------------------------'

  end subroutine csfs_menu

!===========================================================================
  subroutine dets_in_csfs_rd
!---------------------------------------------------------------------------
! Description : read Slater determinants in CSFs
!
! Created     : J. Toulouse, 07 Apr 2009
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'dets_in_csfs_rd'
  character(len=50000) line
  integer iostat, csf_i, det_i, det_in_csf_i, ndet_in_csf_check, normalize_csf
  integer iwdet_in_csf_temp_nb, cdet_csf_temp_nb
  integer, allocatable :: iwdet_in_csf_temp (:)
  real(dp), allocatable :: cdet_csf_temp (:)
  real(dp) :: csf_norm

! begin
  call object_provide ('ncsf')
  call object_provide ('ndet')

  call alloc ('ndet_in_csf', ndet_in_csf, ncsf)

! If all the cdet_in_csf are inputted in integer format (no dots in those lines) then csf_coef are
! assumed to correspond to normalized CSF's and the cdet_in_csf are renormalized so that the
! CSF's are normalized.
! normalize_csf is reset to 0 if any of the cdet_in_csf's are in floating format.
  normalize_csf=1

  do csf_i = 1, ncsf

!  read list of determinants
   read(unit_input,'(a)',iostat=iostat) line
   if (iostat < 0) then
     call die (lhere, 'error while reading dets_in_csfs')
   endif
   
   ndet_in_csf (csf_i) = words_number_in_string (line)
   if (ndet_in_csf(csf_i) > MDET_CSF) then
     call die (lhere, 'ndet_in_csf='+ndet_in_csf(csf_i)+' > MDET_CSF='+MDET_CSF)
   endif

   call alloc ('iwdet_in_csf', iwdet_in_csf, maxval(ndet_in_csf), ncsf)
   call alloc ('cdet_in_csf', cdet_in_csf, maxval(ndet_in_csf), ncsf)

   read(line,*,iostat=iostat) (iwdet_in_csf(det_in_csf_i,csf_i),det_in_csf_i=1,ndet_in_csf(csf_i))
   if (iostat < 0) then
     call die (lhere, 'error while reading dets_in_csfs')
   endif

   do det_in_csf_i=1,ndet_in_csf(csf_i)
     if(iwdet_in_csf(det_in_csf_i,csf_i).gt.ndet) stop 'iwdet_in_csf(det_in_csf_i,csf_i) > ndet'
   enddo

!  read list of coefficients
   read(unit_input,'(a)',iostat=iostat) line
   if (iostat < 0) then
     call die (lhere, 'error while reading dets_in_csfs')
   endif
   ndet_in_csf_check = words_number_in_string (line)

   if(index(line,'.').ne.0) normalize_csf=0
   
   if (ndet_in_csf_check /= ndet_in_csf (csf_i)) then
    write(6,'(a,i5,a)') ' in CSF # ',csf_i,':'
    write(6,'(a,i4,a,i4)') ' the number of determinants = ' ,ndet_in_csf (csf_i),' is different from the number of coefficients =',ndet_in_csf_check
    call die (lhere, ' number of determinants /= number of coefficients in a CSF')
   endif

   read(line,*,iostat=iostat) (cdet_in_csf(det_in_csf_i,csf_i),det_in_csf_i=1,ndet_in_csf(csf_i))
   if (iostat < 0) then
     call die (lhere, 'error while reading dets_in_csfs')
   endif
   
!   call read_next_line_list ('iwdet_in_csf_temp', iwdet_in_csf_temp, iwdet_in_csf_temp_nb)
!   write(6,*) 'iwdet_in_csf_temp=',iwdet_in_csf_temp
!   call read_next_line_list ('cdet_csf_temp', cdet_csf_temp, cdet_csf_temp_nb)
!   write(6,*) 'cdet_csf_temp=',cdet_csf_temp
!   ndet_in_csf (csf_i) = iwdet_in_csf_temp_nb
!   iwdet_in_csf(:, csf_i) = iwdet_in_csf_temp (:)

  enddo

  read(unit_input,'(a)',iostat=iostat) line
!  write(6,*) 'line >',trim(line),'<'
  if (iostat < 0) then
     call die (lhere, 'error while reading dets_in_csfs')
  endif
  if (index(line,'end') == 0) then
     call die (lhere, 'end keyword expected')
  endif

  if(normalize_csf.eq.1) then
!   First calculate normalization and adjust csf_coef to correspond to that.
    write(6,'(a)') ' normalizing cdet_in_csf'
    do csf_i=1,ncsf
       csf_norm=0
       do det_in_csf_i=1,ndet_in_csf(csf_i)
          csf_norm=csf_norm+cdet_in_csf(det_in_csf_i,csf_i)**2
       enddo
       csf_norm=sqrt(csf_norm)
       do det_in_csf_i=1,ndet_in_csf(csf_i)
          cdet_in_csf(det_in_csf_i,csf_i)=cdet_in_csf(det_in_csf_i,csf_i)/csf_norm
       enddo
    enddo
  endif

  call sort_iworbd

  call object_modified ('iworbd')
  call object_modified ('ndet_in_csf')
  call object_modified ('iwdet_in_csf')
  call object_modified ('cdet_in_csf')


  end subroutine dets_in_csfs_rd

! ==============================================================================
  subroutine det_unq_in_csf_bld
! ------------------------------------------------------------------------------
! Description : each csf as a sorted list of unique determinants
!
! Created     : J. Toulouse, 14 Dec 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  integer csf_i, det_i, det_in_csf_i
  integer det_unq_up_cur, det_unq_dn_cur
  real(dp) cdet_unq_cur

! header
  if (header_exe) then

   call object_create ('det_unq_in_csf_nb')
   call object_create ('det_unq_up_in_csf')
   call object_create ('det_unq_dn_in_csf')
   call object_create ('cdet_unq_in_csf')
   call object_create ('csf_prefac')

   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('det_unq_up_in_csf', det_unq_up_in_csf, ncsf)
  call object_alloc ('det_unq_dn_in_csf', det_unq_dn_in_csf, ncsf)
  call object_alloc ('cdet_unq_in_csf', cdet_unq_in_csf, ncsf)
  call object_alloc ('csf_prefac', csf_prefac, ncsf)

  do csf_i = 1, ncsf

   do det_in_csf_i = 1, ndet_in_csf (csf_i)

!    index of original determinant
     det_i = iwdet_in_csf (det_in_csf_i, csf_i)

!    indexes of corresponding unique spin-up and spin-dn determinants
     det_unq_up_cur = det_to_det_unq_up (det_i)
     det_unq_dn_cur = det_to_det_unq_dn (det_i)

!    coefficients of unique determinants, with sign
     cdet_unq_cur = cdet_in_csf (det_in_csf_i, csf_i)

     call append (det_unq_up_in_csf (csf_i)%row, det_unq_up_cur)
     call append (det_unq_dn_in_csf (csf_i)%row, det_unq_dn_cur)
     call append (cdet_unq_in_csf (csf_i)%row, cdet_unq_cur)

   enddo ! det_in_csf_i

!  sort determinants in csf
   call sort (det_unq_up_in_csf (csf_i)%row (:), det_unq_dn_in_csf (csf_i)%row (:), cdet_unq_in_csf (csf_i)%row (:))

!  prefactor of csf = coefficient of first determinant
   csf_prefac (csf_i) = cdet_unq_in_csf (csf_i)%row (1)

!  normalize csf so that the coefficient of first determinant is 1
   cdet_unq_in_csf (csf_i)%row (:) = cdet_unq_in_csf (csf_i)%row (:) / csf_prefac (csf_i)

  enddo ! csf_i

!  write(6,'(2a)') trim(here), ': csfs in terms of sorted unique determinants:'
!  do csf_i = 1, ncsf
!   write(6,'(2a,i4,a)') trim(here), ': csf # ', csf_i,' :'
!   write(6,'(2a,100i3)') trim(here), ': det_unq_up_in_csf =', det_unq_up_in_csf (csf_i)%row (:)
!   write(6,'(2a,100i3)') trim(here), ': det_unq_dn_in_csf =', det_unq_dn_in_csf (csf_i)%row (:)
!   write(6,'(2a,100f7.3)')  trim(here), ': cdet_unq_in_csf =', cdet_unq_in_csf (csf_i)%row (:)
!  enddo

  end subroutine det_unq_in_csf_bld

! ==============================================================================
  subroutine wfdet_ovlp_bld_old
! ------------------------------------------------------------------------------
! Description   : overlap of the determinantal part of the wave function
! Description   : assuming orthnomality of orbitals
!
! Created       : J. Toulouse, 29 Nov 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  integer csf_i, csf_j
  integer det_i, det_j, det_unq_up_i, det_unq_dn_i, det_unq_up_j, det_unq_dn_j
  integer det_in_csf_i, det_in_csf_j

! header
  if (header_exe) then

   call object_create ('wfdet_ovlp')

   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('csf_coef')
   call object_needed ('orb_occ_in_det_unq_up')
   call object_needed ('orb_occ_in_det_unq_dn')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_associate ('wfdet_ovlp', wfdet_ovlp)

  wfdet_ovlp = 0.d0

  do csf_i = 1, ncsf
   do det_in_csf_i = 1, ndet_in_csf (csf_i)
    det_i = iwdet_in_csf (det_in_csf_i, csf_i)
    det_unq_up_i = det_to_det_unq_up (det_i)
    det_unq_dn_i = det_to_det_unq_dn (det_i)

     do csf_j = 1, ncsf
      do det_in_csf_j = 1, ndet_in_csf (csf_j)
       det_j = iwdet_in_csf (det_in_csf_j, csf_j)
       det_unq_up_j = det_to_det_unq_up (det_j)
       det_unq_dn_j = det_to_det_unq_dn (det_j)

        if (arrays_equal (orb_occ_in_det_unq_up (:, det_unq_up_i), orb_occ_in_det_unq_up (:, det_unq_up_j)) .and. &
            arrays_equal (orb_occ_in_det_unq_dn (:, det_unq_dn_i), orb_occ_in_det_unq_dn (:, det_unq_dn_j))) then

          wfdet_ovlp = wfdet_ovlp + csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i) * csf_coef (csf_j, 1) * cdet_in_csf (det_in_csf_j, csf_j)

        endif

       enddo ! det_in_csf_j
      enddo ! csf_j
    enddo ! det_in_csf_i
   enddo ! csf_i

!   write(6,*) trim(here),': wfdet_ovlp=',wfdet_ovlp

  end subroutine wfdet_ovlp_bld_old

! ==============================================================================
  subroutine csfs_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap of the CSFs  (without Jastrow)
! Description   : assuming orthonormality of orbitals
!
! Created       : J. Toulouse, 31 Mar 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  integer csf_i, csf_j
  integer det_i, det_j, det_unq_up_i, det_unq_dn_i, det_unq_up_j, det_unq_dn_j
  integer det_in_csf_i, det_in_csf_j

! header
  if (header_exe) then

   call object_create ('csfs_ovlp')

   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('orb_occ_in_det_unq_up')
   call object_needed ('orb_occ_in_det_unq_dn')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('csfs_ovlp', csfs_ovlp, ncsf, ncsf)

  csfs_ovlp (:,:) = 0.d0

  do csf_i = 1, ncsf
   do det_in_csf_i = 1, ndet_in_csf (csf_i)
    det_i = iwdet_in_csf (det_in_csf_i, csf_i)
    det_unq_up_i = det_to_det_unq_up (det_i)
    det_unq_dn_i = det_to_det_unq_dn (det_i)

     do csf_j = 1, ncsf
      do det_in_csf_j = 1, ndet_in_csf (csf_j)
       det_j = iwdet_in_csf (det_in_csf_j, csf_j)
       det_unq_up_j = det_to_det_unq_up (det_j)
       det_unq_dn_j = det_to_det_unq_dn (det_j)

        if (arrays_equal (orb_occ_in_det_unq_up (:, det_unq_up_i), orb_occ_in_det_unq_up (:, det_unq_up_j)) .and. &
            arrays_equal (orb_occ_in_det_unq_dn (:, det_unq_dn_i), orb_occ_in_det_unq_dn (:, det_unq_dn_j))) then

          csfs_ovlp (csf_i, csf_j) = csfs_ovlp (csf_i, csf_j) + cdet_in_csf (det_in_csf_i, csf_i) * cdet_in_csf (det_in_csf_j, csf_j)

        endif

       enddo ! det_in_csf_j
      enddo ! csf_j
    enddo ! det_in_csf_i
   enddo ! csf_i


  end subroutine csfs_ovlp_bld

! ==============================================================================
  subroutine csfs_wfdet_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap of the CSFs  with the multi-configuration wave function (without Jastrow)
! Description   : assuming orthonormality of orbitals
!
! Created       : J. Toulouse, 31 Mar 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  integer csf_i, csf_j

! header
  if (header_exe) then

   call object_create ('csfs_wfdet_ovlp')

   call object_needed ('ncsf')
   call object_needed ('csfs_ovlp')
   call object_needed ('csf_coef')

   return

  endif

! begin

! allocations
  call object_alloc ('csfs_wfdet_ovlp', csfs_wfdet_ovlp, ncsf)

  csfs_wfdet_ovlp (:) = 0.d0

  do csf_i = 1, ncsf
     do csf_j = 1, ncsf
       csfs_wfdet_ovlp (csf_i) = csfs_wfdet_ovlp (csf_i) + csf_coef (csf_j, 1) * csfs_ovlp (csf_i, csf_j)
     enddo ! csf_j
  enddo ! csf_i


  end subroutine csfs_wfdet_ovlp_bld

! ==============================================================================
  subroutine wfdet_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap of the determinantal part of the wave function
! Description   : assuming orthnomality of orbitals
!
! Created       : J. Toulouse, 31 Mar 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  integer csf_i

! header
  if (header_exe) then

   call object_create ('wfdet_ovlp')

   call object_needed ('ncsf')
   call object_needed ('csfs_wfdet_ovlp')
   call object_needed ('csf_coef')

   return

  endif

! begin

! allocations
  call object_associate ('wfdet_ovlp', wfdet_ovlp)

  wfdet_ovlp = 0.d0

  do csf_i = 1, ncsf
      wfdet_ovlp = wfdet_ovlp + csf_coef (csf_i, 1) * csfs_wfdet_ovlp (csf_i)
  enddo ! csf_i

  end subroutine wfdet_ovlp_bld

! ==============================================================================
  subroutine dens_mat_wfdet_bld
! ------------------------------------------------------------------------------
! Description   : density matrix of determinantal part of wave function
! Description   : assuming orthonormality of orbitals
!
! Created       : J. Toulouse, 29 Nov 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  integer orb_i, orb_j
  integer csf_i, csf_j
  integer det_i, det_unq_up_i, det_unq_dn_i
  integer det_j, det_unq_up_j, det_unq_dn_j
  integer det_in_csf_i, det_in_csf_j
  real(dp) dens_mat_wfdet_up_trace
  real(dp) dens_mat_wfdet_dn_trace
  real(dp) dens_mat_wfdet_trace

! header
  if (header_exe) then

   call object_create ('dens_mat_wfdet_up')
   call object_create ('dens_mat_wfdet_dn')
   call object_create ('dens_mat_wfdet')

   call object_needed ('nelec')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('orb_tot_nb')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('csf_coef')
   call object_needed ('orb_occ_in_adet_unq_up')
   call object_needed ('orb_occ_in_adet_unq_dn')
   call object_needed ('sgn_adet_unq_dn')
   call object_needed ('sgn_adet_unq_up')
   call object_needed ('wfdet_ovlp')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('dens_mat_wfdet_up', dens_mat_wfdet_up, orb_tot_nb, orb_tot_nb)
  call object_alloc ('dens_mat_wfdet_dn', dens_mat_wfdet_dn, orb_tot_nb, orb_tot_nb)
  call object_alloc ('dens_mat_wfdet', dens_mat_wfdet, orb_tot_nb, orb_tot_nb)

  dens_mat_wfdet_up = 0.d0
  dens_mat_wfdet_dn = 0.d0

  do orb_i = 1, orb_tot_nb
  do orb_j = 1, orb_tot_nb

  do csf_i = 1, ncsf
   do det_in_csf_i = 1, ndet_in_csf (csf_i)
    det_i = iwdet_in_csf (det_in_csf_i, csf_i)
    det_unq_up_i = det_to_det_unq_up (det_i)
    det_unq_dn_i = det_to_det_unq_dn (det_i)

     do csf_j = 1, ncsf
      do det_in_csf_j = 1, ndet_in_csf (csf_j)
       det_j = iwdet_in_csf (det_in_csf_j, csf_j)
       det_unq_up_j = det_to_det_unq_up (det_j)
       det_unq_dn_j = det_to_det_unq_dn (det_j)

!       spin-up contribution
        if (arrays_equal (orb_occ_in_adet_unq_up (:, orb_i, det_unq_up_i), orb_occ_in_adet_unq_up (:, orb_j, det_unq_up_j)) .and. &
            arrays_equal (orb_occ_in_det_unq_dn (:, det_unq_dn_i), orb_occ_in_det_unq_dn (:, det_unq_dn_j))) then

          dens_mat_wfdet_up (orb_i, orb_j) = dens_mat_wfdet_up (orb_i, orb_j) +                        &
                                          sgn_adet_unq_up (orb_i, det_unq_up_i) * sgn_adet_unq_up (orb_j, det_unq_up_j)      &
                                        * csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)      &
                                        * csf_coef (csf_j, 1) * cdet_in_csf (det_in_csf_j, csf_j)

        endif

!       spin-down contribution
        if (arrays_equal (orb_occ_in_adet_unq_dn (:, orb_i, det_unq_dn_i), orb_occ_in_adet_unq_dn (:, orb_j, det_unq_dn_j)) .and. &
            arrays_equal (orb_occ_in_det_unq_up (:, det_unq_up_i), orb_occ_in_det_unq_up (:, det_unq_up_j))) then

          dens_mat_wfdet_dn (orb_i, orb_j) = dens_mat_wfdet_dn (orb_i, orb_j) +                       &
                                          sgn_adet_unq_dn (orb_i, det_unq_dn_i) * sgn_adet_unq_dn (orb_j, det_unq_dn_j)     &
                                        * csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)     &
                                        * csf_coef (csf_j, 1) * cdet_in_csf (det_in_csf_j, csf_j)

        endif

       enddo ! det_in_csf_j
      enddo ! csf_j
    enddo ! det_in_csf_i
   enddo ! csf_i

   enddo ! orb_j
   enddo ! orb_i

   dens_mat_wfdet_up = dens_mat_wfdet_up/wfdet_ovlp
   dens_mat_wfdet_dn = dens_mat_wfdet_dn/wfdet_ovlp
   dens_mat_wfdet = dens_mat_wfdet_up + dens_mat_wfdet_dn

   write(6,*) trim(here),': dens_mat_wfdet=',dens_mat_wfdet

!  check trace up
   dens_mat_wfdet_up_trace = 0.d0
   do orb_i = 1, orb_tot_nb
    dens_mat_wfdet_up_trace = dens_mat_wfdet_up_trace + dens_mat_wfdet_up (orb_i, orb_i)
   enddo
   write(6,*) trim(here),': dens_mat_wfdet_up_trace=',dens_mat_wfdet_up_trace

   if (dabs(dens_mat_wfdet_up_trace-nup) > 10d-6) then
    write(6,*) trim(here),': dens_mat_wfdet_up_trace=',dens_mat_wfdet_up_trace,' /= nup=',nup
    call die (here)
   endif

!  check trace down
   dens_mat_wfdet_dn_trace = 0.d0
   do orb_i = 1, orb_tot_nb
    dens_mat_wfdet_dn_trace = dens_mat_wfdet_dn_trace + dens_mat_wfdet_dn (orb_i, orb_i)
   enddo
   write(6,*) trim(here),': dens_mat_wfdet_dn_trace=',dens_mat_wfdet_dn_trace

   if (dabs(dens_mat_wfdet_dn_trace-ndn) > 10d-6) then
    write(6,*) trim(here),': dens_mat_wfdet_dn_trace=',dens_mat_wfdet_dn_trace,' /= ndn=',ndn
    call die (here)
   endif

!  check trace
   dens_mat_wfdet_trace = 0.d0
   do orb_i = 1, orb_tot_nb
    dens_mat_wfdet_trace = dens_mat_wfdet_trace + dens_mat_wfdet (orb_i, orb_i)
   enddo
   write(6,*) trim(here),': dens_mat_wfdet_trace=',dens_mat_wfdet_trace

   if (dabs(dens_mat_wfdet_trace-nelec) > 10d-6) then
    write(6,*) trim(here),': dens_mat_wfdet_trace=',dens_mat_wfdet_trace,' /= nelec=',nelec
    call die (here)
   endif

  end subroutine dens_mat_wfdet_bld

end module csfs_mod

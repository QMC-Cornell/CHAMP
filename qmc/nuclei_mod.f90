module nuclei_mod

  use all_tools_mod

! Declaration of global variables and default values
  logical                   :: l_convert_from_angstrom_to_bohr = .false.
  real(dp), allocatable     :: dist_nn  (:,:)
  real(dp), allocatable     :: mass_nucl  (:)
  real(dp), allocatable     :: mass_nucl_center (:)
  real(dp)                  :: mass_nucl_total
  real(dp), allocatable     :: cent2 (:,:,:)
  real(dp), allocatable     :: cent_sav (:,:)
  real(dp), allocatable     :: cent_best (:,:)

  real(dp), allocatable     :: cent_ref(:,:)
  real(dp), allocatable     :: delc(:,:,:)
  real(dp), allocatable     :: pecentn(:)

  contains

!===========================================================================
  subroutine nuclei_menu
!---------------------------------------------------------------------------
! Description : menu for nuclei
!
! Created     : J. Toulouse, 30 Apr 2008
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'nuclei_menu'
  integer cent_i, dim_i, i, ic, ict, lpotp1_nb

! initialization
  l_convert_from_angstrom_to_bohr = .false.
  write(6,*)
  write(6,'(a)') 'Beginning of nuclei menu ---------------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') ' HELP for nuclei menu:'
   write(6,'(a)') '  nuclei'
   write(6,'(a)') '   ndim = [integer] number of dimensions (default: 3)'
   write(6,'(a)') '   nloc = [integer] type of external potential (default: 0, i.e. -Z/r)'
   write(6,'(a)') '   nquad = [integer] number of quadrature points for pseudopotential (default: 6)'
   write(6,'(a)') '   lpotp1 2 1 ... end : local components of the pseudopotential for each atom type'
   write(6,'(a)') '   geometry ... end: atom types, nuclear charges and cartesian coordinates. Example for H2O:'
   write(6,'(a)') '    geometry'
   write(6,'(a)') '     1 8.0   0.00000000   0.00000000   0.00000000'
   write(6,'(a)') '     2 1.0   0.00000000  -1.43121000   1.10861000'
   write(6,'(a)') '     2 1.0   0.00000000   1.43121000   1.10861000'
   write(6,'(a)') '    end'
   write(6,'(a)') '   convert_from_angstrom_to_bohr = [logical] : convert nuclear coordinates from Angstrom to Bohr units (default=false)'
   write(6,'(a)') '  end'
   write(6,*)

  case ('ndim')
   call get_next_value (ndim)
   call object_modified ('ndim')

  case ('nloc')
   call get_next_value (nloc)
   call object_modified ('nloc')

  case ('nquad')
   call get_next_value (nquad)
   call object_modified ('nquad')

  case ('lpotp1')
   call get_next_value_list ('lpotp1', lpotp1, lpotp1_nb)
   call object_modified ('lpotp1')

  case ('geometry')
   call geometry_rd

  case ('convert_from_angstrom_to_bohr')
   call get_next_value (l_convert_from_angstrom_to_bohr)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  if (use_parser) then
  write(6,'(a,i1)') ' system of dimension ', ndim
  if(ndim.ne.2.and.ndim.ne.3) stop 'ndim must be 2 or 3'
  if(ndim.eq.2.and.iperiodic.gt.0) stop 'ndim=2 not yet implemented for periodic systems'
  if(ndim.eq.2.and.imetro.ne.1.and.index(mode,'vmc').ne.0) stop 'imetro!=1 not yet implemented for ndim=2'

  call object_provide ('nctype')
  call object_provide ('ncent')
  write(6,'(a,i3)') ' number of atomic center types = ', nctype
  write(6,'(a,i5)') ' number of atomic centers = ', ncent
  call object_provide ('iwctype')
  do ic=1,ncent
    if(iwctype(ic).gt.nctype) stop 'iwctype(ic) > nctype'
  enddo

  call object_provide ('cent')
  call object_provide ('znuc')
  write(6,'(a)')' geometry:'
  write(6,'(a)')'                  type  charge            cartesian coordinates'
  do cent_i=1,ncent
    write(6,'(a,i4,a,i4,3x,f5.1,4x,3f12.6)') ' nucleus # ', cent_i,': ', iwctype(cent_i), znuc(iwctype(cent_i)), (cent(dim_i,cent_i),dim_i=1,ndim)
  enddo

! Store the number of centers of each center type in ncent_ctype()
  call alloc ('ncent_ctype', ncent_ctype, nctype)
  do ict=1,nctype
    ncent_ctype(ict)=0
  enddo
  do ic=1,ncent
    ncent_ctype(iwctype(ic))=ncent_ctype(iwctype(ic))+1
!    write(6,'(''Number of centers of each centertype:'',20i3)') (ncent_ctype(ict),ict=1,nctype)
  enddo
  call object_modified ('ncent_ctype')

  write(6,*)
  write(6,'(a,i5)') ' type of external potential: nloc=',nloc
  if (nloc > 0) then
    write(6,'(a)') ' this is a pseudopotential calculation'
  endif

  if(nloc.lt.-5 .or. nloc.gt.6) stop 'nloc must be between -5 and 6 inclusive'

  if(nloc.ge.2) then
    call object_provide ('nquad')
    write(6,'(a,i4)') ' number of quadrature points for pseudopotential = ', nquad
    if (nquad > MPS_QUAD) call die (lhere, 'nquad='+nquad+' > MPS_QUAD='+MPS_QUAD)
  elseif(nloc.eq.-1) then
    call die (lhere, 'nloc=-1 not yet implemented for new input')
  elseif(nloc.eq.-2) then
    call die (lhere, 'nloc=-2 not yet implemented for new input')
  elseif(nloc.eq.-4) then
    call die (lhere, 'nloc=-4 not yet implemented for new input')
  elseif(nloc.eq.-5) then
    call die (lhere, 'nloc=-5 not yet implemented for new input')
  endif

  if(iperiodic.ne.0) then
    call die (lhere, 'new input not yet implemented for periodic calculations.')
  endif

  if(nloc.eq.-3) then
    call die (lhere, 'new input not implemented for Jellium calculations.')
   else
    zconst = 0  ! normal case
  endif

! local components of the pseudopotential
  if(nloc.gt.0) then
    call object_provide ('lpotp1')
    call require (lhere, 'lpotp1_nb = nctype', lpotp1_nb == nctype)
    write(6,'(a,20i3)') ' local components of pseudopotential: ',lpotp1(1:nctype)
!JT    do i=1,nctype
!JT       if(lpotp1(i).gt.MPS_L) stop 'lpotp1(i) > MPS_L'
!JT    enddo
    MPS_L=maxval(lpotp1(1:nctype))  ! JT: set MPS_L from input
  endif

  if(nloc > 0) then
    call alloc ('vps', vps, nelec, ncent, MPS_L)
    call alloc ('npotd', npotd, nctype)
    if(nloc.eq.1) then
      call readps
    elseif(nloc.eq.2.or.nloc.eq.3) then
      call readps_tm
    elseif(nloc.eq.4 .or. nloc.eq.5) then
      call readps_champ
    elseif(nloc.eq.6) then
      call readps_gauss
    else
      stop 'nloc >= 7'
    endif
    do ict=1,nctype
      if(npotd(ict).ge.4 .and. nquad.lt.12) then
        nquad=12
        write(6,'(''Number of quadrature points for psp, nquad, reset to 12 because npotd='',i2)') npotd(ict)
      endif
      if(npotd(ict).ge.5 .and. nquad.lt.24) then
        nquad=24
        write(6,'(''Number of quadrature points for psp, nquad, reset to 24 because npotd='',i2)') npotd(ict)
      endif
      if(npotd(ict).ge.6) then
        write(6,'(a)') ' Warning: we are not ensuring the right number of quadrature points for npotd >=6'
      endif
     enddo
     if (nquad > MPS_QUAD) call die (lhere, 'nquad='+nquad+' > MPS_QUAD='+MPS_QUAD)
     call alloc ('xq0', xq0, MPS_QUAD)
     call alloc ('yq0', yq0, MPS_QUAD)
     call alloc ('zq0', zq0, MPS_QUAD)
     call alloc ('xq', xq, MPS_QUAD)
     call alloc ('yq', yq, MPS_QUAD)
     call alloc ('zq', zq, MPS_QUAD)
     call alloc ('wq', wq, MPS_QUAD)
     call gesqua(nquad,xq0,yq0,zq0,wq)
     if(ipr.ge.0) then
       write(6,'(a)') ' quadrature points for nonlocal pseudopotential:'
       do i=1,nquad
         write(6,'(a,4f10.5)') ' x,y,z,w : ',xq0(i),yq0(i),zq0(i),wq(i)
       enddo
     endif
   endif


  endif ! if use_parser

! convert unit of nuclear coordinates
  if (l_convert_from_angstrom_to_bohr) then
   write(6,*)
   write(6,'(a)') ' Converting nuclear coordinates from Angstrom to Bohr:'
   call object_provide ('ncent')
   call object_provide ('ndim')
   call object_provide ('cent')
   cent(1:ndim,1:ncent) = cent(1:ndim,1:ncent) * angstrom_to_bohr
   call object_modified ('cent')
   call object_provide ('znuc')
   call object_provide ('iwctype')
   write(6,'(a)')'                  type  charge            cartesian coordinates'
   do cent_i=1,ncent
    write(6,'(a,i4,a,i4,3x,f5.1,4x,3f12.6)') ' nucleus # ', cent_i,': ', iwctype(cent_i), znuc(iwctype(cent_i)), (cent(dim_i,cent_i),dim_i=1,ndim)
   enddo
!  recalculate nuclear potential energy
!   call object_provide ('znuc')
!   call object_provide ('iwctype')
!   call pot_nn(cent,znuc,iwctype,ncent,pecent)
!   write(6,'(a,f14.7)') ' recalculting nuclear potential energy: pecent=',pecent
!   call object_modified ('pecent')
  endif

! get nuclear potential energy
  call object_provide ('znuc')
  call object_provide ('iwctype')
  call object_provide ('ncent')
  call pot_nn(cent,znuc,iwctype,ncent,pecent)
  write(6,'(a,f14.7)') ' Nuclear potential energy =',pecent
  call object_modified ('pecent')

  write(6,'(a)') 'End of nuclei menu ---------------------------------------------------------------------------------------'

  end subroutine nuclei_menu

! ==============================================================================
  subroutine geometry_rd
! ------------------------------------------------------------------------------
! Description   : read types, charges and geometry of nuclei
!
! Created       : J. Toulouse, 07 Apr 2009
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'geometry_rd'
  character(len=500) line
  integer iostat, cent_i, dim_i
  real(dp) nuclear_charge

! begin
  call object_provide ('ndim')

  cent_i = 0

  do
   read(unit_input,'(a)',iostat=iostat) line
!   write(6,*) 'line >',trim(line),'<'

!  no next line found
   if(iostat < 0) then
     call die (lhere, 'error while reading geometry')
   endif

!  convert to lower case
   call upplow (line)

!  exit when 'end' is read
   if (index(line,'end') /= 0) exit

   cent_i = cent_i + 1
   call alloc ('cent', cent, 3, cent_i)
   call alloc ('iwctype', iwctype, cent_i)

   read(line,*,iostat=iostat) iwctype(cent_i), nuclear_charge, (cent(dim_i,cent_i),dim_i=1,ndim)
   if(iostat < 0) then
     call die (lhere, 'error while reading geometry')
   endif

   nctype = maxval(iwctype)
   call alloc ('znuc', znuc, nctype)
   znuc(iwctype(cent_i)) = nuclear_charge
  enddo

  ncent = cent_i
  call object_modified ('ncent')
  call object_modified ('nctype')
  call object_modified ('iwctype')
  call object_modified ('znuc')
  call object_modified ('cent')

  end subroutine geometry_rd

! ==============================================================================
  subroutine mass_nucl_bld
! ------------------------------------------------------------------------------
! Description   : atomic mass of each nuclear center
!
! Created       : J. Toulouse, 02 Feb 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer cent_i
  integer znuc_int
  real(dp) :: znuc_to_mass_nucl (10) = (/1, 4, 7, 9, 11, 12, 14, 16, 19, 20/)

! begin

! header
  if (header_exe) then

   call object_create ('mass_nucl')

   call object_needed ('ncent')
   call object_needed ('iwctype')
   call object_needed ('znuc')

   return

  endif

! allocations
  call object_alloc ('mass_nucl', mass_nucl, ncent)

  do cent_i = 1, ncent
   znuc_int = int(znuc(iwctype(cent_i)))
   if (znuc_int >= 1 .and. znuc_int <= 10) then
    mass_nucl (cent_i) = znuc_to_mass_nucl (znuc_int)
   else
    call die (here, 'atomic mass associated to nuclear charge ='+ znuc_int+' is unknown.')
   endif
  enddo ! cent_i


  end subroutine mass_nucl_bld

! ==============================================================================
  subroutine mass_nucl_total_bld
! ------------------------------------------------------------------------------
! Description   : total nuclear mass
!
! Created       : J. Toulouse, 03 Feb 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer cent_i

! begin

! header
  if (header_exe) then

   call object_create ('mass_nucl_total')

   call object_needed ('ncent')
   call object_needed ('mass_nucl')

   return

  endif

! allocations
  mass_nucl_total = 0.d0

  do cent_i = 1, ncent
    mass_nucl_total = mass_nucl_total + mass_nucl (cent_i)
  enddo ! cent_i

  end subroutine mass_nucl_total_bld

! ==============================================================================
  subroutine mass_nucl_center_bld
! ------------------------------------------------------------------------------
! Description   : calculate nuclear center of mass
!
! Created       : J. Toulouse, 03 Feb 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer dim_i, cent_i

! begin

! header
  if (header_exe) then

   call object_create ('mass_nucl_center')

   call object_needed ('ndim')
   call object_needed ('ncent')
   call object_needed ('cent')
   call object_needed ('mass_nucl_total')

   return

  endif

! allocations
  call object_alloc ('mass_nucl_center', mass_nucl_center, ndim)

  do dim_i = 1, ndim
   mass_nucl_center (dim_i) = 0.d0
   do cent_i = 1, ncent
    mass_nucl_center (dim_i) = mass_nucl_center (dim_i) + mass_nucl (cent_i) * cent (dim_i, cent_i)
   enddo ! cent_i
    mass_nucl_center (dim_i) = mass_nucl_center (dim_i) / mass_nucl_total
  enddo ! dim_i

!  write (6,*) "center of mass:", mass_nucl_center(:)

  end subroutine mass_nucl_center_bld

! ==============================================================================
  subroutine dist_nn_bld
! ------------------------------------------------------------------------------
! Description   : distance |Ri - Rj| between two nuclei
!
! Created       : J. Toulouse, 26 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer  cent_i, cent_j, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('dist_nn')

   call object_needed ('ncent')
   call object_needed ('ndim')
   call object_needed ('cent')

   return

  endif

! allocations
  call object_alloc ('dist_nn', dist_nn, ncent, ncent)
  dist_nn (:,:) = 0.d0

  do cent_i = 1, ncent
    do cent_j = cent_i+1, ncent
      do dim_i = 1, ndim
        dist_nn (cent_i, cent_j) = dist_nn (cent_i, cent_j) + (cent(dim_i,cent_i) - cent(dim_i,cent_j))**2
      enddo ! dim_i
      dist_nn (cent_i, cent_j) = dsqrt (dist_nn (cent_i, cent_j))
      dist_nn (cent_j, cent_i) = dist_nn (cent_i, cent_j)
    enddo ! cent_j
  enddo ! cent_i

  end subroutine dist_nn_bld

end module nuclei_mod

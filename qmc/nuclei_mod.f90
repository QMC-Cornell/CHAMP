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
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'nuclei_menu'
  integer cent_i, dim_i, i, ic, it

! initialization
  l_convert_from_angstrom_to_bohr = .false.
  if (use_parser) then
    nforce = 1
  endif

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
   write(6,'(a)') '   nloc = [integer] type of external potential (default: 0)'
   write(6,'(a)') '   nquad = [integer] number of quadrature points for pseudopotential (default: 0)'
   write(6,'(a)') '   nforce = [integer] number of geometries for correlated sampling (default: 1)'
   write(6,'(a)') '   ndim = [integer] number of dimensions (default: 3)'
   write(6,'(a)') '   geometry ... end: atom types, nuclear charges and cartesian coordinates. Example for H2O:'
   write(6,'(a)') '    geometry'
   write(6,'(a)') '     1 8.0   0.00000000   0.00000000   0.00000000'
   write(6,'(a)') '     2 1.0   0.00000000  -1.43121000   1.10861000'
   write(6,'(a)') '     2 1.0   0.00000000   1.43121000   1.10861000'
   write(6,'(a)') '    end'
   write(6,'(a)') '   convert_from_angstrom_to_bohr = [logical] : convert nuclear coordinates from Angstrom to Bohr units (default=false)'
   write(6,'(a)') '  end'
   write(6,*)

  case ('nloc')
   call get_next_value (nloc)
   call object_modified ('nloc')

  case ('nforce')
   call get_next_value (nforce)
   call object_modified ('nforce')
   write(6,'(a,i5)') ' number of geometries for correlated sampling calculation = ',nforce
   call require (lhere, 'nforce <= MFORCE', nforce <= MFORCE)

  case ('ndim')
   call get_next_value (ndim)
   call object_modified ('ndim')

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

  write(6,'(a,i5)') ' type of external potential: nloc=',nloc
  if(nloc.lt.-5 .or. nloc.gt.6) stop 'nloc must be between -5 and 6 inclusive'
  if(nloc.ge.2) then
     write(6,'(a,i4)') ' number of quadrature points for pseudopotential = ', nquad
     if (nquad.gt.MPS_QUAD) stop 'nquad > MPS_QUAD'

! warning quantum dots, wire,... stuff not added here 
!   .... and don't forget nloc.eq.-5 case whenever this is implemented
!  elseif(nloc.eq.-1) then
!        read(5,*) w0,bext,glande
!        we=dsqrt(w0*w0+0.25d0*bext*bext)
!        write(6,'(''spring const of dot pot., w0='',t31,f10.5)') w0
!        write(6,'(''applied magnetic field., bext='',t31,f10.5)') bext
!        write(6,'(''effective spring const., we='',t31,f10.5)') we
!        write(6,'(''Lande factor, glande='',t31,f10.5)') glande
!       elseif(nloc.eq.-2) then
!        read(5,*) p1,p2,p3,p4
!        write(6,'(''quartic dot pot. p1,p2,p3,p4='',t31,9f9.6)') p1,p2,p3,p4
!       elseif(nloc.eq.-4) then 
!        read(5,*) wire_w,wire_length,wire_potential_cutoff
!        we=wire_w  !  this is a quick fix:  needed for the subroutine basis_fns_2dgauss
!        write(6,'(''wire_w,wire_length,wire_potential_cutoff='',t31,9f9.6)') wire_w,wire_length,wire_potential_cutoff
   endif

! Warning: periodic stuff not added yet
!      if(iperiodic.ne.0) then
!
!c npoly is the polynomial order for short-range part
!        read(5,*) npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big
!        write(6,'(/,''Npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big'',2i4,9f8.2)')
!     &   npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big
!        if(npoly.ne.7) then
!          write(6,'(''present version works best with npoly=7'')')
!          stop 'present version works best with npoly=7'
!        endif
!
!        ncoef=npoly+1
!        if(ncoef.gt.NCOEFX) stop 'ncoef gt NCOEFX'
!
!        read(5,*) alattice
!        do 10 i=1,ndim
!          read(5,*) (rlatt(k,i),k=1,ndim)
!          do 10 k=1,ndim
!   10       rlatt(k,i)=rlatt(k,i)*alattice
!
!        write(6,'(/,''Lattice basis vectors'',3(/,3f10.6))')
!     &   ((rlatt(k,j),k=1,ndim),j=1,ndim)
!
!c read the dimensions of the simulation 'cube'
!        do 20 i=1,ndim
!          read(5,*) (rlatt_sim(k,i),k=1,ndim)
!          do 20 k=1,ndim
!   20       rlatt_sim(k,i)=rlatt_sim(k,i)*alattice
!
!        write(6,'(/,''Simulation lattice basis vectors'',3(/,3f10.6))')
!    &   ((rlatt_sim(k,j),k=1,ndim),j=1,ndim)
!
!c read k-shift for generating k-vector lattice
!        read(5,*) (rkvec_shift_latt(k),k=1,ndim)
!        do 22 k=1,ndim
!   22     if(rkvec_shift_latt(k).ne.0.d0 .and. rkvec_shift_latt(k).ne..5d0)
!     &    stop 'rkvec_shift_latt components must be 0 or 1/2 to have real orbs'
!
!      endif

      call object_provide ('nctype')
      call object_provide ('ncent')
      write(6,'(a,i3)') ' number of atomic center types = ', nctype
      write(6,'(a,i5)') ' number of atomic centers = ', ncent

      call object_provide ('iwctype')
      do ic=1,ncent
        if(iwctype(ic).gt.nctype) stop 'iwctype(ic) > nctype'
      enddo

!     warning Jellium stuff not added yet
      if(nloc.eq.-3) then ! Jellium RM
!!MS Jellium sphere plus charge at center
!        dn_background = nelec - znuc(1)
!        rs_jel = 1.d0
!        radius_b = (dn_background*(rs_jel)**3)**(1.d0/3.d0)
!        zconst = 20 !* 27Aug06
       else
        zconst = 0  ! normal case
      endif

! warning pseudo stuff not added yet
!c Read in which is the local component of the potential
!      if(nloc.gt.0) then
!        read(5,*) (lpotp1(i),i=1,nctype)
!        write(6,'(''lpotp1='',t31,20i3,(20i3))') (lpotp1(i),i=1,nctype)
!        do 35 i=1,nctype
!   35     if(lpotp1(i).gt.MPS_L) stop 'lpotp1(i) > MPS_L'
!      endif

      call object_provide ('cent')
      call object_provide ('znuc')
      write(6,'(a)')' geometry:'
      if(iperiodic.eq.0) then
       write(6,'(a)')'                  type  charge            cartesian coordinates'
       do cent_i=1,ncent
        write(6,'(a,i4,a,i4,3x,f5.1,4x,3f12.6)') ' nucleus # ', cent_i,': ', iwctype(cent_i), znuc(iwctype(cent_i)), (cent(dim_i,cent_i),dim_i=1,ndim)
       enddo
      endif


! Convert center positions from primitive lattice vector units to cartesian coordinates
!     if(iperiodic.ne.0) then
!       write(6,'(/,''center positions in primitive lattice vector units and in cartesian coordinates'')')
!       do 66 ic=1,ncent
!         do 62 k=1,ndim
!  62       cent_tmp(k)=cent(k,ic)
!         do 65 k=1,ndim
!           cent(k,ic)=0
!           do 65 i=1,ndim
!  65         cent(k,ic)=cent(k,ic)+cent_tmp(i)*rlatt(k,i)
!  66     write(6,'(''center'',i4,1x,''('',3f9.5,'')'',1x,''('',3f9.5,'')'')') ic, (cent_tmp(k),k=1,ndim),(cent(k,ic),k=1,ndim)
!     endif
!     write(6,*)

!     if(nloc.gt.0) then
!       write(6,'(/,''pseudopotential calculation'')')
!       if(nloc.eq.1) then
!         call readps
!        elseif(nloc.eq.2.or.nloc.eq.3) then
!         call readps_tm
!        elseif(nloc.eq.4 .or. nloc.eq.5) then
!         call readps_champ
!        elseif(nloc.eq.6) then
!         call readps_gauss
!        else
!         stop 'nloc >= 7'
!       endif
!       do 67 ict=1,nctype
!         if(npotd(ict).ge.4 .and. nquad.lt.12) then
!           nquad=12
!           write(6,'(''Number of quadrature points, nquad, reset to 12 because npotd='',i2)') npotd(ict)
!         endif
!         if(npotd(ict).ge.5 .and. nquad.lt.24) then
!           nquad=24
!           write(6,'(''Number of quadrature points, nquad, reset to 24 because npotd='',i2)') npotd(ict)
!         endif
!         if(npotd(ict).ge.6) write(6,'(''Warning: We are not ensuring the right number of quadrature points for npotd >=6'')')
!  67   continue
!       call gesqua(nquad,xq0,yq0,zq0,wq)
!        if(ipr.ge.0) then
!          write(6,'(''Quadrature points for nonlocal pseudopotential'')')
!          do 68 i=1,nquad
!   68       write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)
!        endif
!      endif
!
!      if(iperiodic.ne.0) call set_ewald 
!
! Compute total nuclear charge and compare to number of electrons
! Warn if not equal, stop if they differ by more than 2.
!      znuc_tot=0
!      do 69 ic=1,ncent
!        ict=iwctype(ic)
!   69   znuc_tot=znuc_tot+znuc(ict)
!      if(iperiodic.ne.0) znuc_tot=znuc_tot*vcell_sim/vcell
!      if(znuc_tot.ne.dfloat(nelec)) write(6,'(''znuc_tot='',f6.1,'' != nelec='',i4)') znuc_tot,nelec
!JT      if(abs(znuc_tot-dfloat(nelec)).gt.3) stop 'abs(znuc_tot - nelec) > 3'

!      if(nloc.ne.-3) then ! RM
!        if(abs(znuc_tot-dfloat(nelec)).gt.6) stop 'abs(znuc_tot - nelec) > 6'
!      endif

! TEMPORARY: Warning: we are not calling readforce and only using one geometry
!     if(index(mode,'fit').ne.0) then 
!       nforce=1
!       nwftype=1
!       iwftype(1)=1
!     endif

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
   do cent_i = 1,ncent
    write(6,'(a,i4,a,3f8.5)') ' center # ',cent_i,' : ',cent(1:ndim,cent_i)
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
  include 'modules.h'
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
  include 'modules.h'
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
  include 'modules.h'
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
  include 'modules.h'
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
  include 'modules.h'
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

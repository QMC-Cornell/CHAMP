module correlated_sampling_mod

  use all_tools_mod
  use nuclei_mod

! Declaration of global variables and default values
  real(dp)                  :: alfstr

  contains

!===========================================================================
  subroutine correlated_sampling_menu
!---------------------------------------------------------------------------
! Description : menu for correlated sampling
!
! Created     : J. Toulouse, 09 Mar 2010
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'correlated_sampling_menu'
  integer force_i, cent_i, cent_j, dim_i
  real(dp) rsq, rcm, rsq1

! initialization
  nforce = 1
  istrech = 2
  alfstr = 4.d0

  write(6,*)
  write(6,'(a)') 'Beginning of correlated_sampling menu --------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') ' HELP for correlated_sampling menu:'
   write(6,'(a)') '  correlated_sampling'
   write(6,'(a)') '   nforce = [integer] number of geometries for correlated sampling (default: 1)'
   write(6,'(a)') '   istrech = [integer] type of wrap-space transformation (default: 2)'
   write(6,'(a)') '   alfstr = [real] parameter of wrap transformation (default: 4.)'
   write(6,'(a)') '   geometry_displacements ... end: Example for H2:'
   write(6,'(a)') '    geometry_displacements'
   write(6,'(a)') '     0.00000000   0.00000000   0.00000000'
   write(6,'(a)') '     0.00000000   0.00000000   0.00000000'
   write(6,'(a)') '     0.00000000   0.00000000   0.00000000'
   write(6,'(a)') '     0.00000000   0.00000000   0.00100000'
   write(6,'(a)') '    end'
   write(6,'(a)') '  end'
   write(6,*)

  case ('nforce')
   call get_next_value (nforce)
   call object_modified ('nforce')
   write(6,'(a,i5)') ' number of geometries for correlated sampling calculation: nforce= ',nforce

  case ('istrech')
   call get_next_value (istrech)
   call object_modified ('istrech')
   write(6,'(a,i5)') ' type of wrap-space transformation: istrech= ',istrech

  case ('alfstr')
   call get_next_value (alfstr)
   call object_modified ('alfstr')
   write(6,'(a,f12.6)') ' parameter of wrap-space transformation: alfstr= ',alfstr

  case ('geometry_displacements')
   call geometry_displacements_rd

  case ('nwftype')
   call get_next_value (nwftype)
   call object_modified ('nwftype')
   write(6,'(a,i5)') ' number of wave functions: nwftype=',nwftype

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  call object_provide ('nforce')
  call object_provide ('ncent')
  call object_provide ('ndim')

  if (nwftype /= 1) then
    call die (lhere, 'nwftype /=1 is not yet implemented for new input')
  else
   call alloc ('iwftype', iwftype, nforce)
   do force_i = 1, nforce
     iwftype (force_i) = 1
   enddo
   call object_modified ('iwftype')
  endif

  if (nwftype > nforce) then
   call die (lhere, 'nwftype > nforce')
  endif
  if (iwftype(1) /= 1) then
   call die (lhere, 'iwftype(1) /= 1')
  endif

  if (l_mode_dmc) then
!     if(index(mode,'dmc').ne.0) then
!        read(3,*) nwprod,itausec
!        write(6,'(''nwprod,itausec='',2i4)') nwprod,itausec
!        if(nwprod.lt.1) stop 'nwprod must be 1 or more'
!        if(itausec.ne.0.and.itausec.ne.1) stop 'itausec must be 0 or 1'
!        call object_modified ('nwprod')
!      endif
    call die (lhere, 'correlated sampling not yet implemented for new input')
  endif

  call alloc ('cent_ref', cent_ref, 3, ncent)
  cent_ref (1:ndim,1:ncent) = cent (1:ndim,1:ncent)

  call alloc ('pecentn', pecentn, nforce)
  do force_i = 1, nforce
      call pot_nn(cent(1:ndim,1:ncent)+delc(1:ndim,1:ncent,force_i),znuc,iwctype,ncent,pecentn(force_i))
  enddo
  write(6,'(a,100f10.5)') ' Nuclear-nuclear potential energies = ', pecentn(1:nforce)

  call alloc ('deltot', deltot, nforce)
  do force_i = 1,nforce
    deltot(force_i)=0.d0
    rsq=0.d0
    do cent_j=1,ncent
      do dim_i=1,ndim
        rcm=0.d0
        do cent_i=1,ncent
           rcm=rcm+delc(dim_i,cent_i,force_i)
           rsq=rsq+(cent(dim_i,cent_i)+delc(dim_i,cent_i,force_i)-cent(dim_i,cent_j)-delc(dim_i,cent_j,force_i))**2
        enddo !cent_i
        rcm=rcm/ncent
        deltot(force_i)=deltot(force_i)+(delc(dim_i,cent_j,force_i)-rcm)**2
      enddo ! dim_i
    enddo ! cent_j
    if(force_i == 1) rsq1=rsq
!   Warning: TEMPORARY: multiplication by ncent right for diatomics
    deltot(force_i)=sign(dsqrt(deltot(force_i)*ncent),rsq-rsq1)
!   Warning: this is a temporary fix to put deltot=1 if one is doing excited
!   states rather than forces
    if(deltot(force_i) == 0.d0) deltot(force_i)=1
  enddo ! force_i

  write(6,'(a)') 'End of correlated_sampling menu --------------------------------------------------------------------------'

  end subroutine correlated_sampling_menu

! ==============================================================================
  subroutine geometry_displacements_rd
! ------------------------------------------------------------------------------
! Description   : read geometry displacements
!
! Created       : J. Toulouse, 07 Apr 2009
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'geometry_displacements_rd'
  integer force_i, cent_i

! begin
  call object_provide ('nforce')
  call object_provide ('ndim')
  call object_provide ('ncent')

  call alloc ('delc', delc, 3, ncent, nforce)
  do force_i = 1, nforce
     do cent_i = 1, ncent
       read(5,*) delc (1:ndim, cent_i, force_i)
     enddo
  enddo

  call read_up_to_end

  call object_modified ('delc')

  end subroutine geometry_displacements_rd

end module correlated_sampling_mod

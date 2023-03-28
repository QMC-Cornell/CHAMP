subroutine orb_params
! Written by Cyrus Umrigar
! Figure out number of variational linear and basis exponent parameters, nparml, nparme
! and which those parameters are, iwbasi, iworb, iwbase.
! Also, figure out the number of orbital constraints, necn, and which those contraints are, iebasi, ieorb
! for imposing symmetries while optimizing orbital parameters
! The total number of LCAO parameters is norb*nbasis, where norb is the number of orbs that are occupied in at least one det.
! The number of nonredundant LCAO parameters is smaller because:
! 1) One (or more if there is symmetry) parameter per orb is redundant because normalization does not matter
! 2) For all-electron calculations, the number of parameters is reduced by the number of atoms on which each orbital has
!    a nonzero value, because of imposing cusp conditions on those orbitals and atoms.  This is done in cuspco and cusorb.
! 3) Some parameters are zero by symmetry
! 4) Spatial symmetry also results in equalities among nonzero LCAO coefs.  For Abelian point groups these equalities are only within
!    a given orbital, but for non-Abelian point groups they are also between pairs of degenerate orbitals.
! The number of parameters we get here is not the same as what we get from the standard method implemented by Julien
! which uses the number of
! core (C), active (A), and virtual (V) orbitals of each symmetry.  That gives:
! NC*NA + NC*NV + NA*NV + NA*(NA-1)/2 parameters for each symmetry.
! The number in this routine can be larger, but it can also be smaller because here we use that orbitals in a
! 2-dimensional representation are not independent, whereas Julien does not take that into account.
! This routine relies on equality of coefs. to figure out symmetry.  This equality may be destroyed by a prior non-fit optimization.
! So, the fit optimization should be done on the vanilla input from the quantum chemistry program.  Then ieorb, iebasi printed out can be
! put back in the input and used to do a fit optimization and/or restore the symmetry.
! Note that finding which orbs are closed, active and virtual is not really necessary in this routine.
! necn -> necoef

  use all_tools_mod
  use atom_mod
  use coefs_mod
  use optim_mod
  use numbas_mod
  use lbas_mod
  use basis1_mod
  use basis2_mod
  use basis_mod
  use dets_mod
  use dorb_mod
  use orbitals_mod
  use contrl_opt_mod
  use basis_mod, only : which_analytical_basis

  implicit real*8(a-h,o-z)

  integer, external :: findloc2
  integer coef_max_loc(1)
  character*8, allocatable :: orbital_type(:)
  character*80 fmt
  logical up_occ, dn_occ, closed, virtual, bas_exp_equal
  logical, allocatable :: orb_closed_sym(:,:),coef_varied(:,:),first_basis_fn(:)

  parameter(eps=1.d-6)

! LCAO coefs.
  if(necn < 0) then

    call object_provide ('orb_sym_lab')
    write(6,'(''orb_sym_lab='',1000a10)') orb_sym_lab

!   First we find which orbs are closed, active and virtual, but we use only which are virtual
!   Find closed-shell orbitals (not used because we need to know closed for each symmetry)

    allocate(orbital_type(norb))
    orbital_type='active'

    do iorb=1,norb
      closed=.true.
      do idet=1,ndet
        up_occ=.false.
        dn_occ=.false.
        do ielec=1,nup
          if(iworbd(ielec,idet).eq.iorb) then
            up_occ=.true.
            exit
          endif
        enddo
        do ielec=1,ndn
          if(iworbd(ielec+nup,idet).eq.iorb) then
            dn_occ=.true.
            exit
          endif
        enddo
!       write(6,'(i4,2l)') iorb,up_occ,dn_occ
        if(.not. up_occ .or. .not. dn_occ) then
          closed=.false.
!         write(6,'(i4,3l)') iorb,up_occ,dn_occ,closed
          exit
        endif
      enddo ! idet
!     write(6,*)
      if(closed) orbital_type(iorb)='closed'
    enddo ! iorb

!   Find virtual orbitals
    do iorb=1,norb
      virtual=.true.
      do idet=1,ndet
        up_occ=.false.
        dn_occ=.false.
        do ielec=1,nup
          if(iworbd(ielec,idet).eq.iorb) then
            up_occ=.true.
            exit
          endif
        enddo
        do ielec=1,ndn
          if(iworbd(ielec+nup,idet).eq.iorb) then
            dn_occ=.true.
            exit
          endif
        enddo
        if(up_occ .or. dn_occ) then
          virtual=.false.
          exit
        endif
      enddo ! idet
      if(virtual) orbital_type(iorb)='virtual'
    enddo ! iorb

    write(6,'(i4,x,a)') (iorb,orbital_type(iorb),iorb=1,norb)

!   orb_closed_sym(iorb,jorb) is true if iorb is present in every up/dn det where jorb is present.
!   This allows us to use the invariance property of determinants to reduce the number of free parameters.
    allocate(orb_closed_sym(norb,norb))
    allocate(coef_varied(nbasis,norb))
    orb_closed_sym=.true.
    coef_varied=.true.
    
    do iorb=1,norb
      do jorb=1,norb
        if(orb_sym_lab(jorb) /= orb_sym_lab(iorb)) then
          orb_closed_sym(iorb,jorb)=.false. ! This should never be needed
        else
          do idet=1,ndet
            if((findloc2(iworbd(1:nup,idet),iorb,nup) == 0 .and. findloc2(iworbd(1:nup,idet),jorb,nup) /= 0) .or. &
               (findloc2(iworbd(nup+1:nelec,idet),iorb,ndn) == 0 .and. findloc2(iworbd(nup+1:nelec,idet),jorb,ndn) /= 0)) then ! If jorb present but iorb not present
               orb_closed_sym(iorb,jorb)=.false.
               exit
            endif
          enddo
        endif
      enddo
    enddo

    write(6,'(''orb_closed_sym:'')')
    do iorb=1,norb
      write(6,'(10000l2)') (orb_closed_sym(iorb,jorb),jorb=1,norb)
    enddo

!   Figure out if basis functions are the first basis function on any atom.  These are used for imposing cusp conditions and so are not varied.
    allocate(first_basis_fn(nbasis))
    write(6,'(''imnbas='',1000i4)') imnbas ; call systemflush(6)
    first_basis_fn=.false.
    do ibasis=1,nbasis
      do icent=1,ncent
        if(imnbas(icent).eq.ibasis) then
          write(6,'(''coef('',i4,'',:) not varied because of e-n cusp condition'')') ibasis
          first_basis_fn(ibasis)=.true.
          coef_varied(ibasis,:)=.false.
        endif
      enddo
    enddo
    write(6,'(''first_basis_fn'',1000l2)') first_basis_fn; call systemflush(6)

! If it is the largest abs coef (omitting first basis fn on each atom for Slater basis) or zero do not vary it
    do iorb=1,norb
      coef_max_loc=maxloc(abs(coef(:,iorb,1)),MASK=.not.first_basis_fn)
      do ibasis=1,nbasis
        if(ibasis.eq.coef_max_loc(1) .or. abs(coef(ibasis,iorb,1)).le.eps) then ! If it is the largest coef or zero do not vary it
          if(ibasis.eq.coef_max_loc(1)) write(6,'(''coef('',i4,'','',i4,'') not varied because normalization is irrelevant'')') ibasis,iorb
          if(abs(coef(ibasis,iorb,1)).le.eps) write(6,'(''coef('',i4,'','',i4,'') not varied because coef is zero'')') ibasis,iorb
          coef_varied(ibasis,iorb)=.false.
          cycle
        endif
      enddo
    enddo

!   If orb_closed_sym(iorb,jorb) is true, then one additional parameter in jorb is redundant.
!   One could add iorb to jorb to zero out that element of jorb, but that is optional.
    do iorb=1,norb
!     coef_max_loc=maxloc(abs(coef(:,iorb,1)),MASK=.not.first_basis_fn)
!     write(6,'(''iorb, coef_max_loc='',2i4)') iorb, coef_max_loc(1)
      do jorb=1,norb
        if(iorb.eq.jorb .or. orb_sym_lab(jorb) /= orb_sym_lab(iorb) .or. .not. orb_closed_sym(iorb,jorb)) cycle
        do jbasis=1,nbasis
          if(abs(coef(jbasis,jorb,1)).ge.eps .and. coef_varied(jbasis,jorb)) then
            write(6,'(''coef('',i4,'','',i4,'') not varied because of invariance of determinant'')') jbasis,jorb
            coef_varied(jbasis,jorb)=.false.
            exit
          endif
        enddo
      enddo
    enddo

!   allocate(iworb(norb*(norb-1)/2))
!   allocate(iwbasi(norb*(norb-1)/2))
    write(6,'(''size(iworb),norb*(norb-1)/2)'',2i6)') size(iworb),norb*(norb-1)/2
    call alloc('iwbasi', iwbasi, norb*nbasis)
    call alloc('iworb', iworb, norb*nbasis)
    call alloc ('iebasi', iebasi, 2, norb*nbasis)
    call alloc ('ieorb', ieorb, 2, norb*nbasis)


!   write(6,'(''norb, nbasis='',2i5)') norb, nbasis
    write(6,'(''orb_sym_lab='',1000(a,x))') (trim(orb_sym_lab(iorb)),iorb=1,norb)

!   Figure out which orbital coefs are to others, ieorb, iebasi.  After that figure out which orbital coefs are to be varied iworb, iwbasi.
    write(6,'(''nparml,nparm='',2i5)') nparml,nparm
    nparml_sav=nparml
    nparm_sav=nparm
    nparml=0
    necn=0
    do 25 iorb=1,norb
      if(orbital_type(iorb)=='virtual') cycle
      coef_max_loc=maxloc(abs(coef(:,iorb,1)),MASK=.not.first_basis_fn)
      do 20 ibasis=1,nbasis
        if(abs(coef(ibasis,iorb,1)).lt.eps) cycle
!       For 1d irreps, coefs of a given orb may be equal and for 2d irreps coefs. of different orbs may be equal too.
        do 15 jorb=1,iorb
          if(jorb==iorb) then
            nbasis2=ibasis-1
          else
            nbasis2=nbasis
          endif
          do 10 jbasis=1,nbasis2
            if(zex(ibasis,1).eq.zex(jbasis,1) .and. (jorb==iorb .or. &
              (orb_sym_lab(iorb)(1:1)=='e' .and. orb_sym_lab(jorb)(1:1)=='e') .or. &
              (orb_sym_lab(iorb)(1:1)=='p' .and. orb_sym_lab(jorb)(1:1)=='p') .or. &
              (orb_sym_lab(iorb)(1:1)=='d' .and. orb_sym_lab(jorb)(1:1)=='d') .or. &
              (orb_sym_lab(iorb)(1:1)=='f' .and. orb_sym_lab(jorb)(1:1)=='f') .or. &
              (orb_sym_lab(iorb)(1:1)=='g' .and. orb_sym_lab(jorb)(1:1)=='g') .or. &
              (orb_sym_lab(iorb)(1:1)=='h' .and. orb_sym_lab(jorb)(1:1)=='h') .or. &
              (orb_sym_lab(iorb)(1:1)=='E' .and. orb_sym_lab(jorb)(1:1)=='E') .or. &
              (orb_sym_lab(iorb)(1:1)=='P' .and. orb_sym_lab(jorb)(1:1)=='P') .or. &
              (orb_sym_lab(iorb)(1:1)=='D' .and. orb_sym_lab(jorb)(1:1)=='D') .or. &
              (orb_sym_lab(iorb)(1:1)=='F' .and. orb_sym_lab(jorb)(1:1)=='F') .or. &
              (orb_sym_lab(iorb)(1:1)=='G' .and. orb_sym_lab(jorb)(1:1)=='G') .or. &
              (orb_sym_lab(iorb)(1:1)=='H' .and. orb_sym_lab(jorb)(1:1)=='H'))) then
              if(abs(abs(coef(ibasis,iorb,1))-abs(coef(jbasis,jorb,1))).le.eps .and. abs(coef(ibasis,iorb,1)).gt.eps ) then ! coefs equal or opposite and nonzero
                necn=necn+1
                if(necn.gt.size(ieorb,2)) stop 'necn.gt.size(ieorb,2)'
                coef_varied(ibasis,iorb)=.false.
                iebasi(1,necn)=ibasis
                iebasi(2,necn)=jbasis
                ieorb(1,necn)=iorb
                if(abs(coef(ibasis,iorb,1)-coef(jbasis,jorb,1)).le.eps) then ! coefs equal
                  ieorb(2,necn)=jorb
                  goto 20
                elseif(abs(coef(ibasis,iorb,1)+coef(jbasis,jorb,1)).le.eps) then ! coefs opposite
                  ieorb(2,necn)=-jorb
                  goto 20
                endif
              endif
            endif
10        continue
15      continue
        if(coef_varied(ibasis,iorb)) then
          nparml=nparml+1  ! Since this coef is not equal to another, make it a variational parameter
          if(nparml.gt.size(iworb)) then
            write(6,'(''nparml,size(iworb)'',2i5)') nparml,size(iworb)
            call systemflush(6)
            stop 'nparml.gt.size(iworb)'
          endif
          iwbasi(nparml)=ibasis
          iworb(nparml)=iorb
        endif
20    continue
25  continue

    nparm=nparm_sav+nparml-nparml_sav
    write(6,'(''nparm, nparml, necn='',3i5)') nparm, nparml, necn
!   write(6,'(''2 size(iebasi), size(ieorb),size(iwbasi),size(iworb),'',9i5)') size(iebasi), size(ieorb),size(iwbasi),size(iworb)
    call alloc ('iebasi', iebasi, 2, necn)
    call alloc ('ieorb', ieorb, 2, necn)
    call alloc ('iwbasi', iwbasi, nparml)
    call alloc ('iworb', iworb, nparml)
!   write(6,'(''3 size(iebasi), size(ieorb),size(iwbasi),size(iworb),'',9i5)') size(iebasi), size(ieorb),size(iwbasi),size(iworb)


!   write(6,'(''nparml='',i5)') nparml
!   write(6,'(''iworb,iwbasi='',/,(2i4,x))') (iworb(iparml),iwbasi(iparml),iparml=1,nparml)
!   write(6,'(''necn='',i5)') necn
!   write(6,'(''ieorb,iebasi='',/,2(2i4,x),x)') ((ieorb(i,iecn),iebasi(i,iecn),i=1,2),iecn=1,necn)
!   write(6,'(''lin. coefs of orbs (ieorb,iebasi) set equal='',5(2(2i3,2x),2x))') ((ieorb(i,j),iebasi(i,j),i=1,2),j=1,necn)

    write(6,'(/,2i5,'' nparml,necn'')') nparml,necn
    write(fmt,'(''(''i5,''(2i3,x),a)'')') nparml
    if(nparml.gt.0) write(6,fmt) (iworb(ibasis),iwbasi(ibasis),ibasis=1,nparml),' (iworb(ibasis),iwbasi(ibasis),ibasis=1,nparml)'
    write(fmt,'(''(''i5,''(2(2i3,x),x),a)'')') necn
    if(necn.gt.0) write(6,fmt) ((ieorb(iorb,ibasis),iebasi(iorb,ibasis),iorb=1,2),ibasis=1,necn),' ((ieorb(iorb,ibasis),iebasi(iorb,ibasis),iorb=1,2),ibasis=1,necn)'

  endif ! necn<0

! Basis set exponents
  if(nebase < 0) then
    call alloc ('iwbase', iwbase, nbasis)
    call alloc ('iebase', iebase, 2, nbasis)

    nparme_sav=nparme
    nparme=0
    nebase=0
    do ibasis=1,nbasis
      bas_exp_equal=.false.
      do jbasis=1,ibasis-1
        if(zex(jbasis,1).eq.zex(ibasis,1)) then
          nebase=nebase+1
          iebase(1,nebase)=ibasis
          iebase(2,nebase)=jbasis
          bas_exp_equal=.true.
          exit
        endif
      enddo
      if(.not.bas_exp_equal) then
        nparme=nparme+1  ! Since this coef is not equal to another, make it a variational parameter
        iwbase(nparme)=ibasis
      endif
    enddo

    nparm=nparm+nparme-nparme_sav

    call alloc ('iwbase', iwbase, nparme)
    call alloc ('iebase', iebase, 2, nebase)

    write(6,'(/,2i4,'' nparme,nebase'')') nparme,nebase
    write(fmt,'(''(''i3,''i3,a)'')') nparme
    write(6,fmt) (iwbase(ibasis),ibasis=1,nparme),' (iwbase(ibasis),ibasis=1,nparme)'
    write(fmt,'(''(''i3,''(2i3,2x),a,/)'')') nebase
    write(6,fmt) ((iebase(iorb,ibasis),iorb=1,2),ibasis=1,nebase),' ((iebase(iorb,ibasis),iorb=1,2),ibasis=1,nebase)'
  endif

end subroutine orb_params
!-------------------------------------------------------------------------------------------------------
integer function findloc2(iworbd,iorb,n)

integer, intent(in) :: iworbd(n),iorb,n

findloc2=0
do i=1,size(iworbd)
  if(iworbd(i).eq.iorb) then
    findloc2=i
    exit
  endif
enddo

end function findloc2 

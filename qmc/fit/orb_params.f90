subroutine orb_params
! Written by Cyrus Umrigar
! Figure out number of variational linear and basis exponent parameters, nparml, nparme
! and which those parameters are, iwbasi, iworb, iwbase.
! Also, figure out the number of orbital constraints, necn, and which those contraints are, iebasi, ieorb
! for imposing symmetries while optimizing orbital parameters
! The number of parameters we get here is not the same as what we get from the standard method implemented by Julien
! which uses the number of
! core (C), active (A), and virtual (V) orbitals of each symmetry.  That gives:
! NC*NA + NC*NV + NA*NV + NA*(NA-1)/2 parameters for each symmetry.
! The number in this routine can be larger, but it can also be smaller because here we use that orbitals in a
! 2-dimensional representation are not independent, whereas Julien does not take that into account.
! The number could be further reduced by imposing e-N cusp conditions on orbitals that are nonzero at nuclei.
! This is done in cuspco and cusorb for orbitals for which lo=0.
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

  implicit real*8(a-h,o-z)

  integer, external :: findloc2
  integer coef_max_loc(1)
  character*8, allocatable :: orbital_type(:)
  character*80 fmt
  logical up_occ, dn_occ, closed, virtual
  logical, allocatable :: orb_closed_sym(:,:),first_basis_fn(:)

  parameter(eps=1.d-6)

! LCAO coefs.
  if(necn < 0) then

!   JT: use basis_fns_name instead of lbasis
    call object_provide ('basis_fns_name')
    write(6,'(''basis_fns_name'',1000a10)') basis_fns_name

    call object_provide ('orb_sym_lab')
    write(6,'(''orb_sym_lab='',1000a10)') orb_sym_lab

!   nparml=0
!   necn=0
!   do 20 k=1,norb
!     do 20 i=1,nbasis
!       do 15 l=1,norb
!         do 10 j=1,nbasis
!           if((l-1)*nbasis+j .ge. (k-1)*nbasis+i) goto 15
!           indexi=index(basis_fns_name(i),'_',.true.) !JT
!           indexj=index(basis_fns_name(j),'_',.true.) !JT
!           if(basis_fns_name(i)(indexi+1:indexi+2).eq.basis_fns_name(j)(indexj+1:indexj+2) .and. zex(i,1).eq.zex(j,1)) then
!             if( ((i.ne.j).or.(k.ne.l)) .and. abs(coef(i,k,1)-coef(j,l,1)).le.eps .and. abs(coef(i,k,1)).gt.eps ) then
!               necn=necn+1
!               call alloc ('iebasi', iebasi, 2, necn)
!               call alloc ('ieorb', ieorb, 2, necn)
!               iebasi(1,necn)=i
!               iebasi(2,necn)=j
!               ieorb(1,necn)=k
!               ieorb(2,necn)=l
!               goto 20
!             elseif( ((i.ne.j).or.(k.ne.l)) .and. abs(coef(i,k,1)+coef(j,l,1)).le.eps .and. abs(coef(i,k,1)).gt.eps ) then
!               necn=necn+1
!               call alloc ('iebasi', iebasi, 2, necn)
!               call alloc ('ieorb', ieorb, 2, necn)
!               iebasi(1,necn)=i
!               iebasi(2,necn)=j
!               ieorb(1,necn)=k
!               ieorb(2,necn)=-l
!               goto 20
!             endif
!           endif
!    10   continue
!    15 continue
!    20 continue

    allocate(orbital_type(norb))
    orbital_type='active'

!   First we find which orbs are closed, active and virtual, but we use only which are virtual
!   Find closed-shell orbitals (not used because we need to know closed for each symmetry)
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

!   Find which orbs are closed for each symmetry
!   orb_closed_sym(iorb,jorb) is true if iorb is present in every up/dn det where jorb is present.
    allocate(orb_closed_sym(norb,norb))
    orb_closed_sym=.true.
    
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

    write(6,'(''In orb_param, coef='')')
    do iorb=1,norb
      write(6,'(1000es12.4)') (coef(ibasis,iorb,1),ibasis=1,nbasis)
    enddo

!   Figure out if basis functions are the first basis function on any atom.  These are used for imposiing cusp conditions.
    allocate(first_basis_fn(nbasis))
    write(6,'(''imnbas='',1000i4)') imnbas ; call systemflush(6)

    first_basis_fn=.false.
    do ibasis=1,nbasis
      do icent=1,ncent
        if(imnbas(icent).eq.ibasis) then
          first_basis_fn(ibasis)=.true.
        endif
      enddo
    enddo
    write(6,'(''first_basis_fn'',1000l2)') first_basis_fn; call systemflush(6)


!   If orb_closed_sym(iorb,jorb) is true, then add iorb to jorb to zero out the element of jorb corresponding to the largest absolute element of iorb.
    do iorb=1,norb
      coef_max_loc=maxloc(abs(coef(:,iorb,1)),MASK=.not.first_basis_fn)
      write(6,'(''iorb, coef_max_loc='',2i4)') iorb, coef_max_loc(1)
      write(6,'(''orb(iorb)'',1000f10.6)') coef(:,iorb,1)
      do jorb=1,norb
        if(iorb.eq.jorb .or. orb_sym_lab(jorb) /= orb_sym_lab(iorb) .or. .not. orb_closed_sym(iorb,jorb)) cycle
!       coef(:,iorb,1)=coef(:,iorb,1)/(coef_max_loc(1),iorb) ! This changes the csf_coefs., so maybe omit
        if(abs(coef(coef_max_loc(1),jorb,1)).le.eps) then
          write(6,'(''iorb,jorb,orb_sym_lab(iorb),orb_sym_lab(jorb)'',2i5,2a)') iorb,jorb,orb_sym_lab(iorb),orb_sym_lab(jorb)
          write(6,'(''Warning, coef('',2i4,'') is already zero, rotated by orb'',i4)') coef_max_loc(1),jorb,iorb
        endif
        write(6,'(''adding orb'',i4,'' to orb'',i4,'' to zero out element'',i4)') iorb, jorb, coef_max_loc(1)
        coef(:,jorb,1)=coef(:,jorb,1) - coef(:,iorb,1)*coef(coef_max_loc(1),jorb,1)/coef(coef_max_loc(1),iorb,1)
      enddo
    enddo

    write(6,'(''2'')') ; call systemflush(6)
    do iorb=1,norb
      write(6,'(1000es12.4)') (coef(ibasis,iorb,1),ibasis=1,nbasis)
    enddo

!   allocate(iworb(norb*(norb-1)/2))
!   allocate(iwbasi(norb*(norb-1)/2))
    write(6,'(''size(iworb),norb*(norb-1)/2)'',2i6)') size(iworb),norb*(norb-1)/2
    call alloc('iwbasi', iwbasi, norb*nbasis)
    call alloc('iworb', iworb, norb*nbasis)
    call alloc ('iebasi', iebasi, 2, norb*nbasis)
    call alloc ('ieorb', ieorb, 2, norb*nbasis)
!   allocate(ieorb(2,norb*(norb-1)/2))
!   allocate(iebasi(2,norb*(norb-1)/2))

    write(6,'(''1 size(iebasi), size(ieorb),size(iwbasi),size(iworb),'',9i5)') size(iebasi), size(ieorb),size(iwbasi),size(iworb)

    write(6,'(''norb, nbasis='',2i5)') norb, nbasis
!   write(6,'(''orbital_type='',(a,x))') orbital_type
    write(6,'(''orb_sym_lab='',1000(a,x))') (trim(orb_sym_lab(iorb)),iorb=1,norb)

!   Figure out which orbital coefs are to be varied iworb, iwbasi, and which are equal to others, ieorb, iebasi
    write(6,'(''nparml,nparm='',2i5)') nparml,nparm
    nparml_sav=nparml
    nparm_sav=nparm
    nparml=0
    necn=0
    do 25 iorb=1,norb
      if(orbital_type(iorb)=='virtual') cycle
      coef_max_loc=maxloc(abs(coef(:,iorb,1)),MASK=.not.first_basis_fn)
      do 20 ibasis=1,nbasis
        if(ibasis==coef_max_loc(1)) cycle ! Since normalization does not matter, do not vary largest abs coef
        if(abs(coef(ibasis,iorb,1)).lt.eps) cycle ! If starting coef is very small, we assume it is not a variational parameter
        do 15 jorb=1,iorb
          if(jorb==iorb) then
            nbasis2=ibasis-1
          else
            nbasis2=nbasis
          endif
          do 10 jbasis=1,nbasis2
            if(zex(ibasis,1).eq.zex(jbasis,1) .and. (jorb==iorb .or. &
              (index(orb_sym_lab(iorb),'e').ne.0 .and. index(orb_sym_lab(jorb),'e').ne.0) .or. &
              (index(orb_sym_lab(iorb),'p').ne.0 .and. index(orb_sym_lab(jorb),'p').ne.0) .or. &
              (index(orb_sym_lab(iorb),'d').ne.0 .and. index(orb_sym_lab(jorb),'d').ne.0) .or. &
              (index(orb_sym_lab(iorb),'f').ne.0 .and. index(orb_sym_lab(jorb),'f').ne.0) .or. &
              (index(orb_sym_lab(iorb),'g').ne.0 .and. index(orb_sym_lab(jorb),'g').ne.0) .or. &
              (index(orb_sym_lab(iorb),'h').ne.0 .and. index(orb_sym_lab(jorb),'h').ne.0) .or. &
              (index(orb_sym_lab(iorb),'E').ne.0 .and. index(orb_sym_lab(jorb),'E').ne.0) .or. &
              (index(orb_sym_lab(iorb),'P').ne.0 .and. index(orb_sym_lab(jorb),'P').ne.0) .or. &
              (index(orb_sym_lab(iorb),'D').ne.0 .and. index(orb_sym_lab(jorb),'D').ne.0) .or. &
              (index(orb_sym_lab(iorb),'F').ne.0 .and. index(orb_sym_lab(jorb),'F').ne.0) .or. &
              (index(orb_sym_lab(iorb),'G').ne.0 .and. index(orb_sym_lab(jorb),'G').ne.0) .or. &
              (index(orb_sym_lab(iorb),'H').ne.0 .and. index(orb_sym_lab(jorb),'H').ne.0))) then
              if(abs(coef(ibasis,iorb,1)-coef(jbasis,jorb,1)).le.eps .and. abs(coef(ibasis,iorb,1)).gt.eps ) then
                necn=necn+1
                if(necn.gt.size(ieorb,2)) stop 'necn.gt.size(ieorb,2)'
                iebasi(1,necn)=ibasis
                iebasi(2,necn)=jbasis
                ieorb(1,necn)=iorb
                ieorb(2,necn)=jorb
                write(6,'(''1 iorb, ibasis, jorb, jbasis, coef'',4i4,f10.6)') iorb, ibasis, jorb, jbasis, coef(ibasis,iorb,1)
                goto 20
              elseif(abs(coef(ibasis,iorb,1)+coef(jbasis,jorb,1)).le.eps .and. abs(coef(ibasis,iorb,1)).gt.eps ) then
                necn=necn+1
                if(necn.gt.size(ieorb,2)) stop 'necn.gt.size(ieorb,2)'
                iebasi(1,necn)=ibasis
                iebasi(2,necn)=jbasis
                ieorb(1,necn)=iorb
                ieorb(2,necn)=-jorb
                write(6,'(''2 iorb, ibasis, jorb, jbasis, coef'',4i4,f10.6)') iorb, ibasis, jorb, jbasis, coef(ibasis,iorb,1)
                goto 20
              endif
            endif
10        continue
15      continue
        nparml=nparml+1  ! Since this coef is not equal to another, make it a variational parameter
        if(nparml.gt.size(iworb)) then
          write(6,'(''nparml,size(iworb)'',2i5)') nparml,size(iworb)
          call systemflush(6)
          stop 'nparml.gt.size(iworb)'
        endif
        iwbasi(nparml)=ibasis
        iworb(nparml)=iorb
20    continue
25   continue

    nparm=nparm_sav+nparml-nparml_sav
    write(6,'(''nparm, nparml, necn='',3i5)') nparm, nparml, necn
    write(6,'(''2 size(iebasi), size(ieorb),size(iwbasi),size(iworb),'',9i5)') size(iebasi), size(ieorb),size(iwbasi),size(iworb)
    call alloc ('iebasi', iebasi, 2, necn)
    call alloc ('ieorb', ieorb, 2, necn)
    call alloc ('iwbasi', iwbasi, nparml)
    call alloc ('iworb', iworb, nparml)
    write(6,'(''3 size(iebasi), size(ieorb),size(iwbasi),size(iworb),'',9i5)') size(iebasi), size(ieorb),size(iwbasi),size(iworb)


    write(6,'(''nparml='',i5)') nparml
    write(6,'(''iworb,iwbasi='',/,(2i4,x))') (iworb(iparml),iwbasi(iparml),iparml=1,nparml)
    write(6,'(''necn='',i5)') necn
!   write(6,'(''ieorb,iebasi='',/,2(2i4,x),x)') ((ieorb(i,iecn),iebasi(i,iecn),i=1,2),iecn=1,necn)
    write(6,'(''lin. coefs of orbs (ieorb,iebasi) set equal='',5(2(2i3,2x),2x))') ((ieorb(i,j),iebasi(i,j),i=1,2),j=1,necn)

    write(6,'(2i4,'' nparml,necn'')') nparml,necn
    write(fmt,'(''(''i3,''(2i3,x),a)'')') nparml
    if(nparml.gt.0) write(6,fmt) (iworb(ibasis),iwbasi(ibasis),ibasis=1,nparml),' (iworb(ibasis),iwbasi(ibasis),ibasis=1,nparml)'
    write(fmt,'(''(''i3,''(2(2i3,x),x),a)'')') necn
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
      do jbasis=1,ibasis-1
        if(zex(jbasis,1).eq.zex(ibasis,1)) then
          nebase=nebase+1
          iebase(1,nebase)=ibasis
          iebase(2,nebase)=jbasis
          exit
        endif
      enddo
      nparme=nparme+1
      iwbase(nparme)=ibasis
    enddo

    nparm=nparm+nparme-nparme_sav

    call alloc ('iwbase', iwbase, nparme)
    call alloc ('iebase', iebase, 2, nebase)

    write(6,'(2i4,'' nparme,nebase'')') nparme,nebase
    write(fmt,'(''(''i3,''i3,a)'')') nparme
    write(6,fmt) (iwbase(ibasis),ibasis=1,nparme),' (iwbase(ibasis),ibasis=1,nparme)'
    write(fmt,'(''(''i3,''(2i3,2x),a)'')') nebase
    write(6,fmt) ((iebase(iorb,ibasis),iorb=1,2),ibasis=1,nebase),' ((iebase(iorb,ibasis),iorb=1,2),ibasis=1,nebase)'
    write(6,'(''expon. (iebase) set equal='',10(2i3,2x))') ((iebase(i,j),i=1,2),j=1,nebase)
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

module deriv_fast_mod
! Written by Tyler Anderson
  use all_tools_mod    !MED tools
  use wfsec_mod        !iwf
  use orbitals_mod     !orb_tot_nb
  use psi_mod          !sum_lap_lnj 
  use atom_mod         !pe_cent
  use dim_mod          !ndim
  use const_mod        !hb, nelec
  use vj_mod           !vj
  use dets_mod         !nup, ndn, ncsf, csf_coef, cdet_in_csf, ndet_in_csf
                       !& iwdet_in_csf
  use orb_mod          !orb, dorb, ddorb
  use dorb_mod         !ndetup, ndetdn, iwdorbdn, iwdorbup, iwdetup, iwdetdn
  use eloc_mod         !pe_ee, pe_en, eloc_pot
  use gamma_mod        !gup, gdn, psp_nonloc_pot
  use forcepar_mod     !nwf
  use optim_mod        !iwcsf
  use pseudo_mod       !nloc

  integer, allocatable, target :: ref_orbs_up(:), ref_orbs_dn(:)
  integer, allocatable, target :: irefup(:), irefdn(:)
  integer, allocatable, target :: ivirup(:), ivirdn(:)
  integer, target              :: ind_refup, ind_refdn !index of ref dets

  real(dp), target              :: deta_up, deta_dn     !det(A)
  real(dp), allocatable, target :: aiup(:,:), aidn(:,:) !A^-1
  real(dp), allocatable, target :: atup(:,:), atdn(:,:) !A^tilde
  real(dp), allocatable, target ::  tup(:,:),  tdn(:,:) !A^-1 A^tilde
  real(dp), allocatable, target :: dtup(:,:), dtdn(:,:) !elocal A^-1 A^tilde

  real(dp)                      :: chi
  real(dp), allocatable, target :: yup(:,:), ydn(:,:)
  real(dp), allocatable, target :: dyup(:,:), dydn(:,:)
  real(dp), allocatable, target :: dgup(:,:), dgdn(:,:)
  real(dp), allocatable, target :: dalpha_up(:,:), dalpha_dn(:,:)

  real(dp), allocatable :: elocal_orb_up(:,:), elocal_orb_dn(:,:)
  real(dp)              :: elocal_fast
  real(dp), allocatable :: kin_orb(:,:)
  real(dp), allocatable :: pot_orb(:,:)
  real(dp)              :: kin_fast, pot_fast

  integer                       :: kmaxup, kmaxdn, kmax
  integer , allocatable, target :: lab_refup(:,:), lab_refdn(:,:)
  integer , allocatable, target :: lab_excup(:,:), lab_excdn(:,:)
  integer , allocatable, target :: excup(:,:), excdn(:,:)
  integer , allocatable, target :: refup(:,:), refdn(:,:)
  integer , allocatable, target :: sgnup(:), sgndn(:)
  integer , allocatable, target :: ordup(:), orddn(:)
  real(dp), allocatable, target :: detup(:), detdn(:)
  real(dp), allocatable, target :: invup(:,:), invdn(:,:)
  real(dp), allocatable, target :: deloc_up(:), deloc_dn(:)

  real(dp), allocatable, target :: cdet_in_wf(:,:)

contains

  subroutine xorbs(orbs1, orbs2, norbs, res, n)
!Given two arrays of orbitals orbs1 and orbs2, returns the orbitals
!in orbs1 but not in orbs2.
!e.g.:
!(in)   orbs1 = (/ 1, 3 /)
!(in)   orbs2 = (/ 1, 2 /)
!(out)  res   = (/ 3 /)
!(out)  n     = 1
    implicit none

    integer, intent(in)  :: orbs1(:), orbs2(:), norbs
    integer, intent(out) :: res(:), n

    integer :: i

    n=0
    do i=1, norbs
      if (any(orbs2 .eq. orbs1(i))) cycle
      n=n+1
      res(n) = orbs1(i)
    enddo
  end subroutine

  subroutine index_of(list1, list2, len1, len2, inds)
!Given two arrays list1 and list2, returns the indices of elements of
!list1 in list2. If an element of list1 is not in list2, return -1
!e.g.:
!(in)   list1 = (/ 1, 2, 3, 4 /)
!(in)   list2 = (/ 2, 4 /) 
!(out)  inds  = (/ -1, 1, -1, 2 /)
    implicit none

    integer, intent(in)  :: list1(:), list2(:), len1, len2
    integer, intent(out) :: inds(:)

    integer :: i, j

    inds = -1
    do i=1, len1
      do j=1, len2
        if (list1(i) .eq. list2(j)) then
          inds(i) = j
          exit
        endif
      enddo
    enddo
  end subroutine

  integer function nswaps(exc, ref, len)
!Given exc, an array of orbitals of a particular excited Slater matrix, 
!and ref, an array of orbitals of the reference Slater matrix, return the
!number of swaps required s.t. the bottom k x k submatrix of (A^-1 A_n) 
!is alpha_n.
!e.g.:
!(in)   exc = (/ 1, 2, 3 /)
!(in)   ref = (/ 1, 3, 4 /)
!(out)  nswaps = 1 (In (A^-1 A_n), the first and second columns are 
!columns of I_3. However, the second column is the third column of I_3. 
!Must swap 2nd and 3rd rows of (A^-1 A_n) to account for this.)
    implicit none

    integer, intent(in) :: exc(:), ref(:), len

    integer :: nswp, n, i, inds(len)

    call index_of(exc, ref, len, len, inds)
    n = 1
    nswp = 0
    do i=1, len
      if ((inds(i) .ne. -1)) then
        nswp = nswp + (i - n) !add col swaps
        nswp = nswp + (inds(i) - n) !add row swaps
        n = n + 1
      endif
    enddo
    nswaps = nswp
  end function

  subroutine cdet_in_wf_bld
    implicit none

    real(dp) :: ccsf, cdet
    integer  :: icsf, i, idet, iwf

    if (header_exe) then
      call object_create('cdet_in_wf')
      call object_needed('nwf')
      call object_needed('ncsf')
      call object_needed('csf_coef')
      call object_needed('ndet_in_csf')
      call object_needed('iwdet_in_csf')
      call object_needed('cdet_in_csf')
      return
    endif

    cdet_in_wf(:,:) = 0d0
    do iwf=1, nwf
      do icsf=1, ncsf
        ccsf = csf_coef(icsf, iwf)
        do i=1, ndet_in_csf(icsf)
          idet = iwdet_in_csf(i, icsf)
          cdet =  cdet_in_csf(i, icsf)
          cdet_in_wf(idet, iwf) = cdet_in_wf(idet, iwf) + cdet*ccsf
        enddo !det
      enddo !csf
    enddo !wf
  end subroutine

  subroutine elocal_orb_bld
    implicit none
    
    integer     :: i, j
    real(dp)    :: ddjas
    
    if (header_exe) then
      call object_create('elocal_orb_up')
      call object_create('elocal_orb_dn')
      call object_needed('sum_lap_lnj')
      call object_needed('orb')
      call object_needed('dorb')
      call object_needed('ddorb')
      call object_needed('eloc_pot_loc')
      call object_needed('psp_nonloc_orb')
      return
    endif
    
    ddjas = (sum_lap_lnj + sum(vj*vj))/nelec

    do j=1, orb_tot_nb
      do i=1, nelec
        kin_orb(i,j) = &
        -hb*(ddorb(i,j) + 2*sum(vj(:,i)*dorb(:,i,j)) + ddjas*orb(i,j))

        !Although not strictly necessary, I include the potential
        !in the elocal_orb matrix by spreading the total potential
        !energy among all electrons.
        pot_orb(i,j) = eloc_pot_loc*orb(i,j)/nelec
        if (nloc.gt.0) pot_orb(i,j) = pot_orb(i,j) + psp_nonloc_orb(i,j)
      enddo
    enddo

    elocal_orb_up = kin_orb(1:nup      ,:) + pot_orb(1:nup      ,:)
    elocal_orb_dn = kin_orb(nup+1:nelec,:) + pot_orb(nup+1:nelec,:)
  end subroutine

  subroutine occ_bld
    implicit none

    integer :: iup, idn, iel, orb
    integer :: occup_bld(orb_tot_nb), occdn_bld(orb_tot_nb)

    if (header_exe) then
      call object_create('occup')
      call object_create('occdn')
      call object_create('noccup')
      call object_create('noccdn')
      call object_needed('orb_tot_nb')
      call object_needed('nup')
      call object_needed('ndn')
      call object_needed('ndetup')
      call object_needed('ndetdn')
      call object_needed('iworbdup')
      call object_needed('iworbddn')
      return
    endif

    noccup=0
    do iup=1,ndetup
      do iel=1,nup
        orb=iworbdup(iel,iup)
        if (.not.any(occup_bld(1:noccup).eq.orb)) then
          noccup=noccup+1
          occup_bld(noccup)=orb
        endif
      enddo
    enddo

    noccdn=0
    do idn=1,ndetdn
      do iel=1,ndn
        orb=iworbddn(iel,idn)
        if (.not.any(occdn_bld(1:noccdn).eq.orb)) then
          noccdn=noccdn+1
          occdn_bld(noccdn)=orb
        endif
      enddo
    enddo

    allocate(occup(noccup), occdn(noccdn))
    occup = occup_bld(1:noccup)
    occdn = occdn_bld(1:noccdn)
    call sort(occup)
    call sort(occdn)
  end subroutine

  subroutine excit_bld
!Allocate memory and calculate arrays that remain constant throughout.
    implicit none

    integer :: iup, idn, ne, nswp, ksum
    !type(AlphaMatrix), pointer :: aup, adn

    if (header_exe) then
      call object_create('excitation_info')
      call object_needed('orb_tot_nb')
      call object_needed('nwf')
      call object_needed('nup')
      call object_needed('ndn')
      call object_needed('iworbdup')
      call object_needed('iworbddn')
      call object_needed('occup')
      call object_needed('occdn')
      call object_needed('noccup')
      call object_needed('noccdn')
      return
    endif

    !orbital indexing
    allocate(ref_orbs_up(nup), ref_orbs_dn(ndn)) !reference orbs
    allocate(irefup(nup), irefdn(ndn))
    allocate(ivirup(noccup-nup), ivirdn(noccdn-ndn))
    allocate(ioccup(orb_tot_nb), ioccdn(orb_tot_nb))

    !1b observables
    allocate(aiup(nup,nup), aidn(ndn,ndn))
    allocate(atup(nup,orb_tot_nb), atdn(ndn,orb_tot_nb))
    allocate( tup(nup,orb_tot_nb),  tdn(ndn,orb_tot_nb))
    allocate(dtup(nup,orb_tot_nb), dtdn(ndn,orb_tot_nb))

    !wf matrices 
    allocate(yup(noccup,nup), ydn(noccdn,ndn))
    allocate(dyup(noccup,nup), dydn(noccdn,ndn))
    allocate(gup(noccup,nup), gdn(noccdn,ndn))
    allocate(dgup(noccup,nup), dgdn(noccdn,ndn))
    allocate(dalpha_up(3,ndetup), dalpha_dn(3,ndetdn))
    allocate(cdet_in_wf(ndet, nwf))
    allocate(invup(nup*nup,ndetup), invdn(ndn*ndn,ndetdn))
    allocate(lab_refup(nup,ndetup), lab_refdn(ndn,ndetdn))
    allocate(lab_excup(nup,ndetup), lab_excdn(ndn,ndetdn))
    allocate(excup(nup,ndetup), excdn(ndn,ndetdn))
    allocate(refup(nup,ndetup), refdn(ndn,ndetdn))
    allocate(detup(ndetup)    , detdn(ndetdn))
    allocate(sgnup(ndetup)    , sgndn(ndetdn))
    allocate(ordup(ndetup)    , orddn(ndetdn))
    allocate(deloc_up(ndetup), deloc_dn(ndetdn))

    !elocal
    allocate(elocal_orb_up(nup,orb_tot_nb), elocal_orb_dn(ndn,orb_tot_nb))
    allocate(kin_orb(nelec,orb_tot_nb), pot_orb(nelec,orb_tot_nb))
    
    !nonlocal pseudopotential
    allocate(psp_nonloc_orb(nelec,orb_tot_nb))
    psp_nonloc_orb(:,:) = 0d0
    call object_modified_by_index(psp_nonloc_orb_index)

    !Choose reference det to be the first
    ind_refup = 1 
    ind_refdn = 1
    !Find & store virtual orbitals
    ref_orbs_up = iworbdup(:,ind_refup)
    ref_orbs_dn = iworbddn(:,ind_refdn)
    call index_of(ref_orbs_up, occup, nup, noccup, irefup)
    call index_of(ref_orbs_dn, occdn, ndn, noccdn, irefdn)
    call xorbs((/ (iup, iup=1,noccup) /), irefup, noccup, ivirup, ne)
    call xorbs((/ (idn, idn=1,noccdn) /), irefdn, noccdn, ivirdn, ne)
    call index_of((/ (iup, iup=1,orb_tot_nb) /), occup, orb_tot_nb, &
                  noccup, ioccup)
    call index_of((/ (idn, idn=1,orb_tot_nb) /), occdn, orb_tot_nb, &
                  noccdn, ioccdn)

    kmaxup = 0
    do iup=1,ndetup
      !Find & store ref_orbs, exc_orbs, and iref_orbs
      call xorbs(ref_orbs_up, iworbdup(:,iup), nup, lab_refup(:,iup), ne)
      call xorbs(iworbdup(:,iup), ref_orbs_up, nup, lab_excup(:,iup), ne)
      call index_of(lab_refup(:,iup), ref_orbs_up, ne,    nup, refup(:,iup))
      call index_of(lab_excup(:,iup),       occup, ne, noccup, excup(:,iup))

      nswp = nswaps(iworbdup(:,iup), ref_orbs_up, nup)
      sgnup(iup) = 1-2*mod(nswp, 2)
      ordup(iup) = ne

      kmaxup = max(kmaxup, ne)
    enddo
    
    kmaxdn = 0
    do idn=1,ndetdn
      !Find & store ref_orbs, exc_orbs, and iref_orbs
      call xorbs(ref_orbs_dn, iworbddn(:,idn), ndn, lab_refdn(:,idn), ne)
      call xorbs(iworbddn(:,idn), ref_orbs_dn, ndn, lab_excdn(:,idn), ne)
      call index_of(lab_refdn(:,idn), ref_orbs_dn, ne,    ndn, refdn(:,idn))
      call index_of(lab_excdn(:,idn),       occdn, ne, noccdn, excdn(:,idn))

      nswp = nswaps(iworbddn(:,idn), ref_orbs_dn, ndn)
      sgndn(idn) = 1-2*mod(nswp, 2)
      orddn(idn) = ne

      kmaxdn = max(kmaxdn, ne)
    enddo
    kmax = max(kmaxup, kmaxdn)
  end subroutine

  subroutine t_bld
!Build ai = A^-1 & t = A^-1 A^tilde 
!Done once per MC step
    implicit none

    integer :: i

    if (header_exe) then
      call object_create('tup')
      call object_create('tdn')
      call object_create('aiup')
      call object_create('aidn')
      call object_create('atup')
      call object_create('atdn')
      call object_create('deta_up')
      call object_create('deta_dn')
      call object_needed('excitation_info')
      call object_needed('orb')
      return
    endif

    !TODO: avoid calculating first Nelec cols of t? (just id matrix)
    atup(:,:) = orb(1:nup,:)
    aiup(:,:) = orb(1:nup,ref_orbs_up)
    call matinv(aiup, nup, deta_up)
    tup(:,:) = matmul(aiup, atup)

    atdn(:,:) = orb(nup+1:nelec,:)
    aidn(:,:) = orb(nup+1:nelec,ref_orbs_dn)
    call matinv(aidn, ndn, deta_dn)
    tdn(:,:) = matmul(aidn, atdn)
  end subroutine

  subroutine dt_bld
!Build dt = -A^-1 B t + A^-1 B
    implicit none

    real(dp) :: aib_up(nup, orb_tot_nb), aib_dn(ndn, orb_tot_nb)

    integer :: i

    if (header_exe) then
      call object_create('dtup')
      call object_create('dtdn')
      call object_needed('elocal_orb_up')
      call object_needed('elocal_orb_dn')
      call object_needed('aiup')
      call object_needed('aidn')
      call object_needed('tup')
      call object_needed('tdn')
      return
    endif

    aib_up = matmul(aiup, elocal_orb_up)
    dtup = aib_up - matmul(aib_up(:,ref_orbs_up), tup)

    aib_dn = matmul(aidn, elocal_orb_dn)
    dtdn = aib_dn - matmul(aib_dn(:,ref_orbs_dn), tdn)
  end subroutine

  subroutine alpha_bld
!Build alpha^-1 and det(alpha)
!Done once per MC step
    implicit none

    integer :: iup, idn, kup, kdn, i, j
    real(dp) :: di, dd

    if (header_exe) then
      call object_create('alpha_up')
      call object_create('alpha_dn')
      call object_needed('tup')
      call object_needed('tdn')
      call object_needed('ndetup')
      call object_needed('ndetdn')
      return
    endif

    do iup=1, ndetup
      kup = ordup(iup)
      do i=1,kup
        do j=1,kup
          invup(j+kup*(i-1),iup) = tup(refup(j,iup),lab_excup(i,iup))
        enddo
      enddo

      call matinv(invup(:,iup), kup, detup(iup))
    enddo

    do idn=1, ndetdn
      kdn = orddn(idn)
      do i=1,kdn
        do j=1,kdn
          invdn(j+kdn*(i-1),idn) = tdn(refdn(j,idn),lab_excdn(i,idn))
        enddo
      enddo

      call matinv(invdn(:,idn), kdn, detdn(idn))
    enddo
  end subroutine

  subroutine alpha_deloc_bld
    implicit none

    integer           :: iup, idn, kup, kdn, i, j
    integer , pointer :: exc(:), ref(:)
    real(dp), pointer :: deloc, inv(:)

    real(dp) :: dd
    
    if (header_exe) then
      call object_create('deloc_up')
      call object_create('deloc_dn')
      call object_needed('dtup')
      call object_needed('dtdn')
      call object_needed('alpha_up')
      call object_needed('alpha_dn')
      return
    endif

    do iup=1, ndetup
      kup =  ordup(iup)
      ref => refup(:,iup)
      exc => lab_excup(:,iup)
      inv => invup(:,iup)
      deloc => deloc_up(iup)

      deloc = 0d0
      do j=1,kup
        do i=1,kup
          deloc = deloc + inv(j+kup*(i-1))*dtup(ref(i), exc(j))
        enddo
      enddo
    enddo

    do idn=1, ndetdn
      kdn =  orddn(idn)
      ref => refdn(:,idn)
      exc => lab_excdn(:,idn)
      inv => invdn(:,idn)
      deloc => deloc_dn(idn)

      deloc = 0d0
      do j=1,kdn
        do i=1,kdn
          deloc = deloc + inv(j+kdn*(i-1))*dtdn(ref(i), exc(j))
        enddo
      enddo
    enddo
  end subroutine

  subroutine y_bld
!Build chi and y
!Done once per MC step
    implicit none

    if (header_exe) then
      call object_create('chi')
      call object_create('yup')
      call object_create('ydn')
      call object_needed('alpha_up')
      call object_needed('alpha_dn')
      call object_needed('ncsf')
      call object_needed('csf_coef')
      call object_needed('ndet_in_csf')
      call object_needed('cdet_in_csf')
      call object_needed('iwdet_in_csf')
      call object_needed('iwf')
      call object_needed('ndet')
      call object_needed('cdet_in_wf')
      call object_needed('iwdetup')
      call object_needed('iwdetdn')
      return
    endif

    call eval_y(invup, invdn, detup, detdn, chi, yup, ydn)
  end subroutine y_bld

  subroutine eval_y(invup, invdn, detup, detdn, chi, yup, ydn)
    implicit none
 
    real(dp), intent(in) , target :: invup(:,:), invdn(:,:)
    real(dp), intent(in) , target :: detup(:)  , detdn(:)
    real(dp), intent(out)         :: chi, yup(:,:), ydn(:,:)

    real(dp) :: ccsf, cdet, coef, chii
    integer  :: icsf, id, i, j, idet, iup, idn, sgn, k

    real(dp) :: cup(ndetup), cdn(ndetdn)

    chi = 0d0
    cup = 0d0
    cdn = 0d0
    do idet=1, ndet
      cdet = cdet_in_wf(idet, iwf) 
      iup  = iwdetup(idet)
      idn  = iwdetdn(idet)
      sgn  = sgnup(iup)*sgndn(idn)
      coef = sgn*cdet*detup(iup)*detdn(idn)

      chi = chi + coef

      cup(iup) = cup(iup) + coef
      cdn(idn) = cdn(idn) + coef
    enddo

    yup = 0d0
    do iup=1, ndetup
      k = ordup(iup)
      do i=1,k
        do j=1,k
          yup(excup(j,iup),refup(i,iup)) = &
          yup(excup(j,iup),refup(i,iup)) + cup(iup)*invup(j+k*(i-1),iup)
        enddo
      enddo
    enddo !iup

    ydn = 0d0
    do idn=1, ndetdn
      k = orddn(idn)
      do i=1,k
        do j=1,k
          ydn(excdn(j,idn),refdn(i,idn)) = &
          ydn(excdn(j,idn),refdn(i,idn)) + cdn(idn)*invdn(j+k*(i-1),idn)
        enddo
      enddo
    enddo !idn

    chii = 1d0/chi
    yup = yup*chii
    ydn = ydn*chii
  end subroutine

  subroutine g_bld
    implicit none

    real(dp) :: yaiup(noccup, nup), yaidn(noccdn, ndn)

    real(dp) :: tyaiup(nup,nup), tyaidn(ndn,ndn)
    integer :: i, ii

    if (header_exe) then
      call object_create('gup')
      call object_create('gdn')
      call object_needed('yup')
      call object_needed('ydn')
      return
    endif

    yaiup = matmul(yup, aiup)
    tyaiup = matmul(tup(:,occup), yaiup)
    gup(irefup,:) = aiup - tyaiup
    gup(ivirup,:) = yaiup(ivirup,:)
     
    yaidn = matmul(ydn, aidn)
    tyaidn = matmul(tdn(:,occdn), yaidn)
    gdn(irefdn,:) = aidn - tyaidn
    gdn(ivirdn,:) = yaidn(ivirdn,:)
  end subroutine

  subroutine fmatmul(prod, mat1, mat2, d)
    implicit none

    real(dp), intent(out) :: prod(:)
    real(dp), intent(in)  :: mat1(:), mat2(:)
    integer , intent(in)  :: d
    
    integer  :: i,j,k
    
    do i=1,d
      do j=1,d
        prod(j+d*(i-1)) = 0d0
        do k=1,d
          prod(j+d*(i-1)) = &
          prod(j+d*(i-1)) + &
          mat1(j+d*(k-1))*mat2(k+d*(i-1))
        enddo
      enddo
    enddo
  end subroutine

  subroutine dy_bld
    implicit none

    real(dp)          :: ccsf, cdet, coef, esum, deloc
    real(dp)          :: fdt(kmax*kmax), idti(kmax*kmax), dti(kmax*kmax)
    real(dp)          :: chii 
    integer           :: icsf, id, idet, iup, idn, sgn, i, j, k
    integer , pointer :: ref(:), exc(:)
    real(dp), pointer :: inv(:)

    if (header_exe) then
      call object_create('dyup')
      call object_create('dydn')
      call object_needed('yup')
      call object_needed('ydn')
      call object_needed('dtup')
      call object_needed('dtdn')
      call object_needed('deloc_up')
      call object_needed('deloc_dn')
      return
    endif
    
    dyup = 0d0
    dydn = 0d0
    do icsf=1, ncsf
      ccsf = csf_coef(icsf, iwf)
      do id=1, ndet_in_csf(icsf)
        idet = iwdet_in_csf(id, icsf)
        cdet =  cdet_in_csf(id, icsf)
        iup  = iwdetup(idet)
        idn  = iwdetdn(idet)
        sgn  = sgnup(iup)*sgndn(idn)
        coef = sgn*cdet*ccsf*detup(iup)*detdn(idn)
        deloc = deloc_up(iup) + deloc_dn(idn)

        k   =  ordup(iup)
        exc => lab_excup(:,iup)
        ref => refup(:,iup)
        inv => invup(:,iup)

        do i=1,k
          do j=1,k
            fdt(j+k*(i-1)) = dtup(ref(j), exc(i))
          enddo
        enddo

        call fmatmul( dti, fdt, inv, k)
        call fmatmul(idti, inv, dti, k)
        exc => excup(:,iup)
        do i=1,k
          do j=1,k
            dyup(exc(j),ref(i)) = &
            dyup(exc(j),ref(i)) + &
            coef*(deloc*inv(j+k*(i-1)) - idti(j+k*(i-1)))
          enddo
        enddo

        k   =  orddn(idn)
        exc => lab_excdn(:,idn)
        ref => refdn(:,idn)
        inv => invdn(:,idn)

        do i=1,k
          do j=1,k
            fdt(j+k*(i-1)) = dtdn(ref(j), exc(i))
          enddo
        enddo

        call fmatmul( dti, fdt, inv, k)
        call fmatmul(idti, inv, dti, k)
        exc => excdn(:,idn)
        do i=1,k
          do j=1,k
            dydn(exc(j),ref(i)) = &
            dydn(exc(j),ref(i)) + &
            coef*(deloc*inv(j+k*(i-1)) - idti(j+k*(i-1)))
          enddo
        enddo
      enddo !idet
    enddo !icsf

    chii = 1d0/chi
    esum = sum(yup*transpose(dtup(:,occup))) &
         + sum(ydn*transpose(dtdn(:,occdn)))
    dyup = dyup*chii - esum*yup
    dydn = dydn*chii - esum*ydn
  end subroutine

  subroutine dg_bld
    implicit none

    real(dp) :: aux1u(nup,nup), aux1d(ndn,ndn)
    real(dp) :: aux2u(nup,nup), aux2d(ndn,ndn)
    real(dp) :: aux3u(noccup,nup), aux3d(noccdn,ndn)
    real(dp) :: aux4u(noccup,nup), aux4d(noccdn,ndn)

    integer :: i

    if (header_exe) then
      call object_create('dgup')
      call object_create('dgdn')
      call object_needed('gup')
      call object_needed('gdn')
      call object_needed('dyup')  
      call object_needed('dydn')  
      return
    endif

    aux1u = matmul(tup(:,occup), dyup) 
    aux2u = matmul(dtup(:,occup), yup)
    aux3u = matmul(gup, elocal_orb_up(:,ref_orbs_up))
    aux4u = 0d0
    aux4u(irefup,:) = aux1u + aux2u
    aux4u = dyup - aux3u - aux4u
    dgup = matmul(aux4u, aiup)

!    dgup = 0d0
!    dgup(irefup,:) = matmul(tup(:,occup), dyup) + matmul(dtup(:,occup), yup)
!    dgup = matmul(-matmul(gup, elocal_orb_up(:,ref_orbs_up)) + dyup - dgup, aiup) 

    aux1d = matmul(tdn(:,occdn), dydn) 
    aux2d = matmul(dtdn(:,occdn), ydn)
    aux3d = matmul(gdn, elocal_orb_dn(:,ref_orbs_dn))
    aux4d = 0d0
    aux4d(irefdn,:) = aux1d + aux2d
    aux4d = dydn - aux3d - aux4d
    dgdn = matmul(aux4d, aidn)

!    dgdn = 0d0
!    dgdn(irefdn,:) = matmul(tdn(:,occdn), dydn) + matmul(dtdn(:,occdn), ydn)
!    dgdn = matmul(-matmul(gdn, elocal_orb_dn(:,ref_orbs_dn)) + dydn - dgdn, aidn) 
  end subroutine

  subroutine csf_derivs_bld
    use slater_mod !deti_det
    use delocc_mod       !denergy
    implicit none

    real(dp) :: etot, cdet, coef, deloc, chii
    integer  :: iparm, icsf, i, idet, iup, idn, sgn

    if (header_exe) then
      call object_create('denergy')
      call object_create('deti_det')
      call object_needed('nparmcsf')
      call object_needed('deloc_up')
      call object_needed('deloc_dn')
      call object_needed('dtup')
      call object_needed('dtdn')
      call object_needed('yup')
      call object_needed('ydn')
      call object_needed('chi')
      return
    endif

    chii = 1d0/chi
    etot=sum(yup*transpose(dtup(:,occup)))+sum(ydn*transpose(dtdn(:,occdn)))
    do iparm=1, nparmcsf
      icsf=iwcsf(iparm)
      deti_det(iparm)=0d0
      denergy(iparm) =0d0
      do i=1, ndet_in_csf(icsf)
        idet = iwdet_in_csf(i, icsf)
        cdet =  cdet_in_csf(i, icsf)
        iup  = iwdetup(idet)
        idn  = iwdetdn(idet)
        sgn  = sgnup(iup)*sgndn(idn)

        coef  = sgn*cdet*detup(iup)*detdn(idn)*chii
        deloc = deloc_up(iup) + deloc_dn(idn)

        deti_det(iparm) = deti_det(iparm)+coef
        denergy (iparm) =  denergy(iparm)+coef*(deloc - etot)
      enddo !idet in csf
    enddo !iparm
  end subroutine

  subroutine elocal_fast_bld
    implicit none

    if (header_exe) then
      call object_create('elocal_fast')
      call object_needed('gup')
      call object_needed('gdn')
      call object_needed('elocal_orb_up')
      call object_needed('elocal_orb_dn')
      return
    endif

    elocal_fast = sum(gup*transpose(elocal_orb_up(:,occup))) &
                + sum(gdn*transpose(elocal_orb_dn(:,occdn)))
  end subroutine

  subroutine psp_nonloc_pot_bld
    implicit none

    if (header_exe) then
      call object_create('psp_nonloc_pot')
      call object_needed('gup')
      call object_needed('gdn')
      call object_needed('psp_nonloc_orb')
      return
    endif

    psp_nonloc_pot = sum(gup*transpose(psp_nonloc_orb(1:nup,occup))) &
                   + sum(gdn*transpose(psp_nonloc_orb(nup+1:nelec,occdn)))
  end subroutine
end module deriv_fast_mod

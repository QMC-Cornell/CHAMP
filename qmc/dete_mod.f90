module dete_mod
  use deriv_fast_mod
  use coefs_mod
  use orbe_mod

  !tmp values for sweep
!  real(dp), target              :: deta_upn, deta_dnn
!  real(dp), allocatable, target :: aiupn(:,:), aidnn(:,:)
!  real(dp), allocatable, target :: tupn(:,:), tdnn(:,:)
!  real(dp), allocatable, target :: yupn(:,:), ydnn(:,:)
!  real(dp), target              :: chin
!
!  real(dp), allocatable, target :: detupn(:), detdnn(:)
!  real(dp), allocatable, target :: invupn(:,:), invdnn(:,:)

  real(dp), pointer :: deta_upn, deta_dnn
  real(dp), pointer :: aiupn(:,:), aidnn(:,:)
  real(dp), pointer :: tupn(:,:), tdnn(:,:)
  real(dp), pointer :: yupn(:,:), ydnn(:,:)
  real(dp), pointer :: chin

  real(dp), pointer :: detupn(:), detdnn(:)
  real(dp), pointer :: invupn(:,:), invdnn(:,:)

contains

  subroutine dete_init(iel)
    implicit none
    
    integer, intent(in) :: iel

    integer :: idet, ord, i, j

!    if (.not.allocated(aiupn))  allocate(aiupn(nup,nup))
!    if (.not.allocated(tupn))   allocate(tupn(nup,noccup))
!    if (.not.allocated(aidnn))  allocate(aidnn(ndn,ndn))
!    if (.not.allocated(tdnn))   allocate(tdnn(ndn,noccdn))
!    if (.not.allocated(invdnn)) allocate(invdnn(ndn*ndn,ndetdn))
!    if (.not.allocated(invupn)) allocate(invupn(nup*nup,ndetup))
!    if (.not.allocated(detdnn)) allocate(detdnn(ndetdn))
!    if (.not.allocated(detupn)) allocate(detupn(ndetup))

    if (iel.le.nup) then
      aiupn = aiup
      deta_upn = deta_up
      tupn = tup(:,occup)
      detupn = detup
      do idet=1,ndetup
        ord = ordup(idet)
        do i=1,ord
          do j=1,ord
            invupn(j+ord*(i-1),idet) = invup(j+ord*(i-1),idet)
          enddo
        enddo
      enddo
    else
      aidnn = aidn
      deta_dnn = deta_dn
      tdnn = tdn(:,occdn)
      detdnn = detdn
      do idet=1,ndetdn
        ord = orddn(idet)
        do i=1,ord
          do j=1,ord
            invdnn(j+ord*(i-1),idet) = invdn(j+ord*(i-1),idet)
          enddo
        enddo
      enddo
    endif
    chin = chi
    yupn = yup
    ydnn = ydn
  end subroutine

  subroutine dete_save(iel)
    implicit none

    integer, intent(in) :: iel

    integer :: idet, ord, i, j

    if (iel.le.nup) then
      aiup = aiupn
      deta_up = deta_upn
      tup(:,occup) = tupn
      detup = detupn
      do idet=1,ndetup
        ord = ordup(idet)
        do i=1,ord
          do j=1,ord
            invup(j+ord*(i-1),idet) = invupn(j+ord*(i-1),idet)
          enddo
        enddo
      enddo
    else
      aidn = aidnn
      deta_dn = deta_dnn
      tdn(:,occdn) = tdnn
      detdn = detdnn
      do idet=1,ndetdn
        ord = orddn(idet)
        do i=1,ord
          do j=1,ord
            invdn(j+ord*(i-1),idet) = invdnn(j+ord*(i-1),idet)
          enddo
        enddo
      enddo
    endif
    chi = chin
    yup = yupn
    ydn = ydnn
  end subroutine

  subroutine dete_update(iel)
    implicit none

    integer, intent(in) :: iel

    integer, pointer           :: ref(:), exc(:), ref_orbs(:)
    integer                    :: offset, nel, ndet, i, j, k, ii, jj, nocc
    real(dp), pointer          :: ai(:,:), t(:,:), deta
!    real(dp)                   :: v(orb_tot_nb), vi2(orb_tot_nb), vaiu, dd, rat
    real(dp)                   :: vaiu, dd, rat
    real(dp), allocatable      :: aiu(:), vai(:), ui(:), vi1(:), v(:), vi2(:)
    real(dp), pointer          :: dets(:), invs(:,:), det, inv(:)
    integer , pointer          :: refs(:,:), excs(:,:), ord(:), occ(:)

    call dete_init(iel)

!assign a few aliases
    if (iel.le.nup) then
      offset   =  0
      nel      =  nup
      ndet     =  ndetup
      ref_orbs => ref_orbs_up
      ai       => aiupn
      deta     => deta_upn
      t        => tupn
      refs     => refup
      excs     => excup
      ord      => ordup
      dets     => detupn
      invs     => invupn
      occ      => occup
      nocc     =  noccup
      allocate(ui(nup), vi1(nup), v(noccup), vi2(noccup))
    else
      offset   =  nup
      nel      =  ndn
      ndet     =  ndetdn
      ref_orbs => ref_orbs_dn
      ai       => aidnn
      deta     => deta_dnn
      t        => tdnn
      refs     => refdn
      excs     => excdn
      ord      => orddn
      dets     => detdnn
      invs     => invdnn
      occ      => occdn
      nocc     =  noccdn
      allocate(ui(ndn), vi1(ndn), v(noccdn), vi2(noccdn))
    endif

!needed for Sherman-Morrison
!A' = A + u v
!A_inv' = A_inv - (ai u v ai)/(1 + v ai u)
!define ui := -ai u / (1 + vaiu)
!define vi := v ai
!A_inv' = A_inv + ui vi
    rat = dot_product(orbe(ref_orbs), ai(:,iel-offset))
    ui  = -ai(:,iel-offset)/rat
    vi1 = matmul(orbe(ref_orbs), ai)
    vi1(iel-offset) = vi1(iel-offset) - 1d0
    vi2 = matmul(orbe(ref_orbs),  t) - orbe(occ)

    do i=1,nel
      do j=1,nel
        ai(j,i) = ai(j,i) + ui(j)*vi1(i) 
      enddo 
    enddo 
 
    deta = deta*rat
 
    !update t w/ Sherman-Morrison 
    do i=1,nocc
      do j=1,nel
        t(j,i) = t(j,i) + ui(j)*vi2(i) 
      enddo 
    enddo 

    do i=1,ndet
      k   =   ord(i)
      ref => refs(1:k,i)
      exc => excs(1:k,i)
      inv => invs(1:k*k,i) 
      det => dets(i)

      do ii=1,k
        do jj=1,k
          inv(jj+k*(ii-1)) = t(ref(jj),exc(ii))
        enddo
      enddo

      call matinv(inv, k, det)
!      call matinve(inv, k, det, ui(ref), vi2(exc))
    enddo  

    if (iel.le.nup) then
      call eval_y(invupn, invdn, detupn, detdn, chin, yupn, ydnn)
    else
      call eval_y(invup, invdnn, detup, detdnn, chin, yupn, ydnn)
    endif
  end subroutine

  subroutine matinve(inv, n, det, u, v) 
!TODO: Old routine using Sherman-Morrison to calculate the rank-1 update
!to an alpha matrix. Has numerical problems (division by (1 + viu)).
!description: updates inverse and det of a matrix M (inv) if M' = M + uv^T  
!w/ Sherman-Morrison 
    implicit none 
 
    integer , intent(in)    :: n 
    real(dp), intent(in)    :: u(n), v(n) 
    real(dp), intent(inout) :: inv(n*n), det 
 
    integer  :: i, j
    real(dp) :: iu(n), vi(n), viu

    if (n.eq.0) return

    do i=1,n
      do j=1,n
        iu(j) = iu(j) + inv(j+n*(i-1))*u(i)
        vi(i) = vi(i) + v(j)*inv(j+n*(i-1))
      enddo
      viu = viu + vi(i)*u(i)
    enddo
 
    do i=1,n 
      do j=1,n 
        inv(j+n*(i-1)) = inv(j+n*(i-1)) - iu(j)*vi(i)/(1 + viu) 
      enddo 
    enddo 
 
    det = det*(1 + viu) 
  end subroutine

  subroutine eval_grad(iel, dorb, ai, tmat, y, grad)
    implicit none

    integer               :: iel
    real(dp), intent(in)  :: dorb(:,:), ai(:,:), tmat(:,:), y(:,:)
    real(dp), intent(out) :: grad(3)
    
    integer :: nel, offset, i, j, k, nocc
    real(dp), allocatable :: v(:)
    
    if (iel.le.nup) then
      allocate(v(noccup))
      do k=1,ndim
        grad(k) = sum(ai(:,iel)*dorb(k,ref_orbs_up))

        if (size(tmat,2).eq.orb_tot_nb) then
          v = dorb(k,occup) - matmul(dorb(k,ref_orbs_up), tmat(:,occup))
        else
          v = dorb(k,occup) - matmul(dorb(k,ref_orbs_up), tmat)
        endif
        do i=1,nup
          do j=1,noccup
              grad(k) = grad(k) + y(j,i)*ai(i,iel)*v(j)
          enddo
        enddo
      enddo
    else
      allocate(v(noccdn))
      do k=1,ndim
        grad(k) = sum(ai(:,iel-nup)*dorb(k,ref_orbs_dn))

        if (size(tmat,2).eq.orb_tot_nb) then
          v = dorb(k,occdn) - matmul(dorb(k,ref_orbs_dn), tmat(:,occdn))
        else
          v = dorb(k,occdn) - matmul(dorb(k,ref_orbs_dn), tmat)
        endif
        do i=1,ndn
          do j=1,noccdn
              grad(k) = grad(k) + y(j,i)*ai(i,iel-nup)*v(j)
          enddo
        enddo
      enddo
    endif
  end subroutine

  subroutine eval_lapl(iel, ddorb, ai, tmat, y, lapl)
    implicit none

    integer               :: iel
    real(dp), intent(in)  :: ddorb(:), ai(:,:), tmat(:,:), y(:,:)
    real(dp), intent(out) :: lapl
    
    integer :: i, j
    real(dp), allocatable :: v(:)
    
    if (iel.le.nup) then
      allocate(v(noccup))
      lapl = sum(ai(:,iel)*ddorb(ref_orbs_up))

      if (size(tmat,2).eq.orb_tot_nb) then
        v = ddorb(occup) - matmul(ddorb(ref_orbs_up), tmat(:,occup))
      else
        v = ddorb(occup) - matmul(ddorb(ref_orbs_up), tmat)
      endif
      do i=1,nup
        do j=1,noccup
          lapl = lapl + y(j,i)*ai(i,iel)*v(j)
        enddo
      enddo
    else
      allocate(v(noccdn))
      lapl = sum(ai(:,iel-nup)*ddorb(ref_orbs_dn))

      if (size(tmat,2).eq.orb_tot_nb) then
        v = ddorb(occdn) - matmul(ddorb(ref_orbs_dn), tmat(:,occdn))
      else
        v = ddorb(occdn) - matmul(ddorb(ref_orbs_dn), tmat)
      endif
      do i=1,ndn
        do j=1,noccdn
          lapl = lapl + y(j,i)*ai(i,iel-nup)*v(j)
        enddo
      enddo
    endif
  end subroutine

end module dete_mod

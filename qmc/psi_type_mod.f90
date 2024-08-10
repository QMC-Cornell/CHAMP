module psi_type_mod
    use types_mod, only: dp
    implicit none

    type psi_t
        real(dp)                :: det
        real(dp)                :: jas
        real(dp), allocatable   :: grad(:,:)
        real(dp)                :: lapl
        real(dp), allocatable   :: hess(:,:,:)
        real(dp)                :: eloc
        real(dp)                :: ekinpb
        real(dp)                :: ekinjf
        real(dp)                :: epot
        real(dp)                :: epot_ee

        !jastrow
        real(dp)                :: fsum, d2
        real(dp), allocatable   :: fj(:,:), fs(:,:), fij(:,:,:), d2ij(:,:)
        real(dp), allocatable   :: lapj(:), lapjij(:,:)

        !determinant
        real(dp), allocatable   :: orb(:,:), dorb(:,:,:), ddorb(:,:)
        real(dp), allocatable   :: aiup(:,:), aidn(:,:)
        real(dp), allocatable   :: tup(:,:), tdn(:,:)
        real(dp), allocatable   :: detup(:), detdn(:) 
        real(dp), allocatable   :: invup(:,:), invdn(:,:)
        real(dp), allocatable   :: yup(:,:), ydn(:,:)
        real(dp)                :: deta_up, deta_dn, chi
        
        !other
        real(dp), allocatable   :: rvec_en(:,:,:), r_en(:,:)
        real(dp), allocatable   :: rvec_ee(:,:), r_ee(:)
        real(dp), allocatable   :: quadr(:,:), quadx(:,:,:)
        integer,  allocatable   :: iwfragelec(:)
        real(dp), allocatable   :: enefrag(:)
    contains
        procedure               :: init => psi_init
    end type

contains

    subroutine psi_init(psi)
        use dim_mod, only: ndim
        use const_mod, only: nelec
        use dorb_mod, only: ndetup, ndetdn
        use gamma_mod, only: noccup, noccdn
        use dets_mod, only: nup, ndn
        use orbitals_mod, only: orb_tot_nb
        use qua_mod, only: nquad
        use atom_mod, only: ncent
        use fragments_mod, only: nfrag
        use all_tools_mod
        implicit none
        class(psi_t), intent(out)    :: psi


        call object_provide('noccup')
        call object_provide('noccdn')

        allocate(psi%grad(ndim, nelec))
        allocate(psi%hess(ndim, ndim, nelec))

        !jastrow
        allocate(psi%fj(ndim, nelec))
        allocate(psi%fs(nelec, nelec))
        allocate(psi%fij(ndim, nelec, nelec))
        allocate(psi%d2ij(nelec, nelec))
        allocate(psi%lapj(nelec))
        allocate(psi%lapjij(nelec, nelec))

        !determinant
        allocate(psi%orb(nelec, orb_tot_nb))
        allocate(psi%dorb(ndim, nelec, orb_tot_nb))
        allocate(psi%ddorb(nelec, orb_tot_nb))
        allocate(psi%aiup(nup,nup), psi%aidn(ndn,ndn))                   
        allocate(psi%tup(nup,orb_tot_nb),  psi%tdn(ndn,orb_tot_nb))             
        allocate(psi%detup(ndetup), psi%detdn(ndetdn))
        allocate(psi%invup(nup*nup,ndetup), psi%invdn(ndn*ndn,ndetdn))
        allocate(psi%yup(noccup,nup), psi%ydn(noccdn,ndn))

        !other
        allocate(psi%rvec_en(ndim, nelec, ncent), psi%r_en(nelec, ncent))
        allocate(psi%rvec_ee(ndim, nelec*(nelec-1)/2), psi%r_ee(nelec*(nelec-1)/2))
        allocate(psi%quadr(ncent*nquad,nelec), psi%quadx(ndim,ncent*nquad,nelec))
        allocate(psi%iwfragelec(nelec))
        allocate(psi%enefrag(nfrag+1))
    end subroutine psi_init

end module psi_type_mod

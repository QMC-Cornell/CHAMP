module walker_type_mod
    use constants_mod, only: dp
    use psi_utils_mod
    implicit none

!    real(dp) :: wmin = 1d0/2
!    real(dp) :: wmax = 2.0d0

    type walker_t
        real(dp)                :: weight
        real(dp), allocatable   :: x(:,:)
        type(psi_t)             :: psi
        integer                 :: age
    contains
        procedure, pass         :: walker_init1 
        procedure, pass         :: walker_init2
        generic                 :: init => walker_init1, walker_init2
    end type

contains

    subroutine walker_init1(walker)
        use dim_mod, only: ndim
        use const_mod, only: nelec
        implicit none
        class(walker_t), intent(inout)  :: walker

        walker%weight=0
        allocate(walker%x(ndim, nelec))
        call walker%psi%init
        walker%age=0
    end subroutine walker_init1

    subroutine walker_init2(walker,x,weight)
        use dim_mod, only: ndim
        use const_mod, only: nelec
        use variables_mod, only: l_mode_dmc
        implicit none
        class(walker_t), intent(inout)  :: walker
        real(dp), intent(in)            :: x(:,:)
        real(dp), intent(in)            :: weight
        logical :: tmp

        walker%weight=weight
        allocate(walker%x(ndim, nelec))
        walker%x=x
        call walker%psi%init
!        tmp=l_mode_dmc
!        l_mode_dmc=.FALSE.
        call psi_at(walker%x,walker%psi)
!        l_mode_dmc=tmp
        walker%age=0
    end subroutine walker_init2

!    function split_joined(walkers)
!        implicit none
!        type(walker_t), intent(in)  :: walkers(:)
!        type(walker_t), allocatable :: split_joined(:)
!        integer, allocatable        :: is(:)
!        real(dp), allocatable       :: ws(:)
!        integer                     :: itmp, ic, iw
!        integer                     :: nwalker_max, nwalker, ncopy
!        real(dp)                    :: wtmp, wtot, wt, rannyu
!        real(dp)                    :: wtot2
!
!        wtot=0
!        do iw=1,size(walkers)
!            wtot=wtot+walkers(iw)%weight
!        enddo
!        nwalker_max=ceiling(wtot/wmin)
!        allocate(is(nwalker_max))
!        allocate(ws(nwalker_max))
!
!        nwalker=0
!        wtmp=0
!        do iw=1,size(walkers)
!            wt=walkers(iw)%weight
!            if (wt.GT.wmax) then !split
!                ncopy=max(1, floor(wt))
!                do ic=1,ncopy
!                    call append_walker(iw,wt/ncopy)
!                enddo
!            else if (wt.LT.wmin) then !join
!                wtmp=wtmp+wt
!                if (rannyu(0).LT.wt/wtmp) itmp=iw
!                if (wtmp.GT.wmin) then
!                    call append_walker(itmp,wtmp)
!                    wtmp=0
!                endif
!            else
!                call append_walker(iw,wt)
!            endif
!        enddo
!
!        if (wtmp.GT.0) call append_walker(itmp,wtmp)
!
!        allocate(split_joined(nwalker))
!        do iw=1,nwalker
!            split_joined(iw)=walkers(is(iw))
!            split_joined(iw)%weight=ws(iw)
!        enddo
!    contains
!        subroutine append_walker(iw,wt)
!            implicit none
!            integer, intent(in) :: iw
!            real(dp), intent(in) :: wt
!
!            nwalker=nwalker+1
!            is(nwalker)=iw
!            ws(nwalker)=wt
!        end subroutine append_walker
!
!    end function split_joined
end module walker_type_mod

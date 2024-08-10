module stats_type_mod
!#if defined(MPI)
!    use mpi
!#endif
    use constants_mod!, only: dp
    implicit none

    type stats_t
        integer                 :: nest

        real(dp), allocatable    :: est_s1av(:) 
        real(dp), allocatable    :: est_bsav(:) 
        real(dp), allocatable    :: est_bssq(:) 
        real(dp), allocatable    :: est_tsav(:) 
        real(dp), allocatable    :: est_tssq(:) 
        real(dp), allocatable    :: est_tbsq(:) 

        real(dp), allocatable    :: wgt_s1av(:)
        real(dp), allocatable    :: wgt_bsav(:)
        real(dp), allocatable    :: wgt_bssq(:)
        real(dp), allocatable    :: wgt_tsav(:)
        real(dp), allocatable    :: wgt_tssq(:)
        real(dp), allocatable    :: wgt_tbsq(:)
    contains
        procedure               :: init             => stats_init
        procedure               :: clear            => stats_clear
        procedure               :: accuest1         => stats_accuest1
        procedure               :: end_step         => stats_end_step
        procedure               :: accuest          => stats_accuest
        procedure               :: end_block        => stats_end_block
        procedure               :: calculate        => stats_calculate
        procedure               :: avg              => stats_avg
        procedure               :: block_avg        => stats_block_avg
        procedure               :: nstep_eff        => stats_nstep_eff
        procedure               :: nblck_eff        => stats_nblck_eff
        procedure               :: sigma            => stats_sigma
    end type

contains
    
    subroutine stats_init(stats, nest)
        implicit none
        class(stats_t), intent(out)  :: stats
        integer, intent(in)         :: nest

        stats%nest = nest

        allocate(stats%est_s1av(nest)) 
        allocate(stats%est_bsav(nest)) 
        allocate(stats%est_bssq(nest)) 
        allocate(stats%est_tsav(nest)) 
        allocate(stats%est_tssq(nest)) 
        allocate(stats%est_tbsq(nest)) 

        allocate(stats%wgt_s1av(nest)) 
        allocate(stats%wgt_bsav(nest)) 
        allocate(stats%wgt_bssq(nest)) 
        allocate(stats%wgt_tsav(nest)) 
        allocate(stats%wgt_tssq(nest)) 
        allocate(stats%wgt_tbsq(nest)) 

        stats%est_s1av = 0
        stats%est_bsav = 0
        stats%est_bssq = 0
        stats%est_tsav = 0
        stats%est_tssq = 0
        stats%est_tbsq = 0

        stats%wgt_s1av = 0
        stats%wgt_bsav = 0
        stats%wgt_bssq = 0
        stats%wgt_tsav = 0
        stats%wgt_tssq = 0
        stats%wgt_tbsq = 0
    end subroutine stats_init

    subroutine stats_clear(stats)
        implicit none
        class(stats_t), intent(inout)  :: stats

        stats%est_s1av = 0
        stats%est_bsav = 0
        stats%est_bssq = 0
        stats%est_tsav = 0
        stats%est_tssq = 0
        stats%est_tbsq = 0

        stats%wgt_s1av = 0
        stats%wgt_bsav = 0
        stats%wgt_bssq = 0
        stats%wgt_tsav = 0
        stats%wgt_tssq = 0
        stats%wgt_tbsq = 0
    end subroutine stats_clear

    subroutine stats_accuest1(stats, iest, wgt, val)
        implicit none
        class(stats_t), intent(inout)   :: stats
        integer, intent(in)             :: iest
        real(dp), intent(in)            :: wgt, val

        stats%wgt_s1av(iest) = stats%wgt_s1av(iest) + wgt
        stats%est_s1av(iest) = stats%est_s1av(iest) + wgt*val
    end subroutine stats_accuest1

    subroutine stats_end_step(stats)
        use variables_mod, only: l_mode_dmc_mov1_mpi2
        implicit none
        class(stats_t), intent(inout)   :: stats
        real(dp)                        :: wgt, val
        integer                         :: i, IERROR

        if (l_mode_dmc_mov1_mpi2) then
            call MPI_Allreduce(MPI_IN_PLACE,stats%wgt_s1av,size(stats%wgt_s1av),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
            call MPI_Allreduce(MPI_IN_PLACE,stats%est_s1av,size(stats%est_s1av),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        endif

        do i=1,stats%nest
            wgt=stats%wgt_s1av(i)
            if (wgt.GT.0) then
                val=stats%est_s1av(i)/wgt
                call stats % accuest(i, wgt, val)
            endif
        enddo
        stats%wgt_s1av=0
        stats%est_s1av=0
    end subroutine stats_end_step

    subroutine stats_accuest(stats, iest, wgt, val)
        implicit none
        class(stats_t), intent(inout)   :: stats
        integer, intent(in)             :: iest
        real(dp), intent(in)            :: wgt, val

        stats%wgt_bsav(iest) = stats%wgt_bsav(iest) + wgt
        stats%wgt_bssq(iest) = stats%wgt_bssq(iest) + wgt**2
        stats%est_bsav(iest) = stats%est_bsav(iest) + wgt*val
        stats%est_bssq(iest) = stats%est_bssq(iest) + wgt*val**2
    end subroutine stats_accuest

!    subroutine stats_end_block(stats, bavg, bwgt, avg, sig, err, tcorr)
    subroutine stats_end_block(stats, bavg, bwgt, avg, err, err1)
        use mpi_mod, only: idtask
        use iterat_mod, only: ipass
        use stats_index_mod
        use variables_mod, only: l_mode_dmc_mov1_mpi1
        implicit none
        class(stats_t), intent(inout)   :: stats
!        real(dp), intent(out)           :: bavg(:), bwgt(:), avg(:), sig(:), err(:), tcorr(:)
        real(dp), intent(out)           :: bavg(:), bwgt(:), avg(:), err(:), err1(:)
!        real(dp)                        :: err1(size(err))
        real(dp)                        :: wbsq(size(stats%wgt_bsav))
        real(dp)                        :: ebsq(size(stats%est_bsav))
        integer                         :: IERROR
        integer   :: unit

        wbsq=stats%wgt_bsav**2
        ebsq=stats%est_bsav**2/stats%wgt_bsav

#if defined(MPI)
        if (l_mode_dmc_mov1_mpi1) then
            call MPI_Allreduce(MPI_IN_PLACE,wbsq          ,size(wbsq)          ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
            call MPI_Allreduce(MPI_IN_PLACE,ebsq          ,size(ebsq)          ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
            call MPI_Allreduce(MPI_IN_PLACE,stats%wgt_bsav,size(stats%wgt_bsav),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
            call MPI_Allreduce(MPI_IN_PLACE,stats%wgt_bssq,size(stats%wgt_bssq),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
            call MPI_Allreduce(MPI_IN_PLACE,stats%est_bsav,size(stats%est_bsav),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
            call MPI_Allreduce(MPI_IN_PLACE,stats%est_bssq,size(stats%est_bssq),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        endif
#endif

        bwgt = stats%wgt_bsav
        bavg = stats%est_bsav/stats%wgt_bsav

        stats%est_tsav = stats%est_tsav+stats%est_bsav
        stats%est_tssq = stats%est_tssq+stats%est_bssq
        stats%est_tbsq = stats%est_tbsq+ebsq
        stats%est_bsav = 0
        stats%est_bssq = 0

        stats%wgt_tsav = stats%wgt_tsav+stats%wgt_bsav
        stats%wgt_tssq = stats%wgt_tssq+stats%wgt_bssq
        stats%wgt_tbsq = stats%wgt_tbsq+wbsq
        stats%wgt_bsav = 0
        stats%wgt_bssq = 0

        call stats % calculate(avg, err, err1)
    end subroutine stats_end_block

    subroutine stats_calculate(stats, avg, err, err1)
        implicit none
        class(stats_t), intent(in)  :: stats
        real(dp), intent(inout)     :: avg(:), err(:), err1(:)
        real(dp)                    :: wsav(size(avg)), wssq(size(avg)), wbsq(size(avg))
        real(dp)                    :: essq(size(avg)), ebsq(size(avg))                     
        real(dp)                    :: nstep_eff(size(avg)), nblck_eff(size(avg)), sig_block(size(avg))

        if (any(size(avg).NE.(/ size(err), size(err1) /))) then
            stop 'inconsistent size in stats_calculate'
        endif

        wsav = stats%wgt_tsav + stats%wgt_bsav
        wssq = stats%wgt_tssq + stats%wgt_bssq
        wbsq = stats%wgt_tbsq + stats%wgt_bsav**2

        avg = (stats%est_tsav + stats%est_bsav)/wsav
        essq = stats%est_tssq + stats%est_bssq
        ebsq = stats%est_tbsq
        if (any(stats%wgt_bsav.GT.0d0)) ebsq = ebsq + stats%est_bsav**2/stats%wgt_bsav

        nstep_eff = wsav**2/wssq
        nblck_eff = wsav**2/wbsq

!        if (any(nblck_eff.LE.1).OR.any((ebsq/wsav-avg**2).LT.0)) then
        if (any(nblck_eff.LE.1)) then
            write(6,'(''Warning: undefined block sigma in stats_calculate. Setting it to zero.'')')
            err = 0 !TODO: tmp just so fpe-trap doesn't kill the program
        else
            err = dsqrt(1/(nblck_eff-1))*dsqrt(ebsq/wsav-avg**2)
        endif
!        if (any(nstep_eff.LE.1).OR.any((essq/wsav-avg**2).LT.0)) then
        if (any(nstep_eff.LE.1)) then
            write(6,'(''Warning: undefined sigma in stats_calculate. Setting it to zero.'')')
            err1 = 0 !TODO: tmp just so fpe-trap doesn't kill the program
        else
            err1 = dsqrt(1/(nstep_eff-1))*dsqrt(essq/wsav-avg**2)
        endif
    end subroutine stats_calculate

    real(dp) function stats_block_avg(stats, iest)
        implicit none
        class(stats_t), intent(in)  :: stats
        integer, intent(in)         :: iest

        stats_block_avg=stats%est_bsav(iest)/stats%wgt_bsav(iest)
    end function stats_block_avg

    function stats_avg(stats, iest) result(avg)
        implicit none
        class(stats_t), intent(in)  :: stats
        integer, intent(in)         :: iest
        real(dp)                    :: avg, cum, wgt

        cum=stats%est_tsav(iest)+stats%est_bsav(iest)
        wgt=stats%wgt_tsav(iest)+stats%wgt_bsav(iest)
        avg=cum/wgt
    end function stats_avg

    function stats_nstep_eff(stats, iest) result(nstep_eff)
        implicit none
        class(stats_t), intent(in)  :: stats
        integer, intent(in)         :: iest
        real(dp)                    :: nstep_eff, wsav, wssq

        wsav = stats%wgt_tsav(iest) + stats%wgt_bsav(iest)
        wssq = stats%wgt_tssq(iest) + stats%wgt_bssq(iest)
        nstep_eff = wsav**2/wssq
    end function stats_nstep_eff

    function stats_nblck_eff(stats, iest) result(nblck_eff)
        implicit none
        class(stats_t), intent(in)  :: stats
        integer, intent(in)         :: iest
        real(dp)                    :: nblck_eff, wsav, wbsq

        wsav = stats%wgt_tsav(iest) + stats%wgt_bsav(iest)
        wbsq = stats%wgt_tbsq(iest) + stats%wgt_bsav(iest)**2
        nblck_eff = wsav**2/wbsq
    end function stats_nblck_eff

    function stats_sigma(stats, iest) result(sigma)
        implicit none
        class(stats_t), intent(in)  :: stats
        integer, intent(in)         :: iest
        real(dp)                    :: sigma, wsav, avg, essq

        wsav = stats%wgt_tsav(iest) + stats%wgt_bsav(iest)
        avg = (stats%est_tsav(iest) + stats%est_bsav(iest))/wsav
        essq = stats%est_tssq(iest) + stats%est_bssq(iest)
        sigma = dsqrt(max(0d0, essq/wsav-avg**2))
    end function stats_sigma

end module stats_type_mod

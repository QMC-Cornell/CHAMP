!    subroutine split_join3(walkers, nwalk)
!    subroutine split_join3(walkers)
!        use walkers_t_mod, only: walker_t
!        type(walker_t), intent(inout)   :: walkers(:)
!!        integer, intent(inout)          :: nwalk
!        return
!    end subroutine split_join3

    subroutine test_subroutine1(num1, array, num2)
        real(8) :: array(:), num1, num2
        if (size(array).GT.0) write(6,'(f12.4)') array(1)
        write(6,'(2f12.4)') num1, num2
    end subroutine test_subroutine1

    subroutine test_subroutine2(num)
        real(8) :: num
        write(6,'(f12.4)') num
    end subroutine test_subroutine2

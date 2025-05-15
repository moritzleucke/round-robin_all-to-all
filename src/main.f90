program test
    use round_robin
    use mpi_env
    implicit none

    integer :: n_ops, i, n_ranks
    integer, dimension(:), allocatable :: schedule

    call init_mpi()

    print *, 'Hello from rank', mpi_rank, 'of', mpi_size


    n_ranks = 14

    do i = 0, n_ranks-1
        call get_round_robin_turnament_schedule(n_ranks, i, n_ops, schedule)
        print *, schedule(:)
    end do


    call finalize_mpi()
end program test
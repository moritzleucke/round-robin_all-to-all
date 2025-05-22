program test
    use round_robin
    use mpi_env
    use distributed_sum, only: sum_and_redistribute
    implicit none


    call init_mpi()
    call main()
    call finalize_mpi()

    !n_ranks = 14
    !do i = 0, n_ranks-1
    !    call get_round_robin_turnament_schedule(n_ranks, i, n_ops, schedule)
    !    print *, schedule(:)
    !end do

    ! TODO: write the sum + redistribution



contains

    subroutine main()

        integer, parameter :: n_dim1_glob = 20
        integer, parameter :: n_dim2_glob = 20
        integer, dimension(:), allocatable :: idx_list_1, idx_list_2
        integer, dimension(:), allocatable :: end_idx_list_1, end_idx_list_2
        integer :: n_dim1_loc, n_dim2_loc
        integer :: end_n_dim1_loc, end_n_dim2_loc
        real(kind=8), dimension(:, :), allocatable :: loc_mat
        real(kind=8), dimension(:, :), allocatable :: end_mat

        call build_local_matrix(n_dim1_glob, n_dim2_glob, n_dim1_loc, n_dim2_loc, &
                                  idx_list_1, idx_list_2, loc_mat)
        call get_end_datalayout(n_dim1_glob, n_dim2_glob, end_n_dim1_loc, end_n_dim2_loc, &
                                  end_idx_list_1, end_idx_list_2)
        call sum_and_redistribute(idx_list_1, idx_list_2, loc_mat, &
                                  end_idx_list_1, end_idx_list_2, end_mat)

    end subroutine main



    subroutine build_local_matrix(n_dim1_glob, n_dim2_glob, n_dim1_loc, n_dim2_loc, &
                                  idx_list_1, idx_list_2, loc_mat)
        integer, intent(in) :: n_dim1_glob, n_dim2_glob
        integer, intent(out) :: n_dim1_loc, n_dim2_loc
        integer, dimension(:), allocatable, intent(out) :: idx_list_1
        integer, dimension(:), allocatable, intent(out) :: idx_list_2
        real(kind=8), dimension(:, :), allocatable, intent(out) :: loc_mat

        ! internal variables
        integer :: i_dim1, i_dim2, counter_1, counter_2

        n_dim1_loc = n_dim1_glob - mpi_size +1
        n_dim2_loc = n_dim2_glob - mpi_size +1
        allocate(loc_mat(n_dim1_loc, n_dim2_loc))
        allocate(idx_list_1(n_dim1_loc))
        allocate(idx_list_2(n_dim2_loc))
        !print *, mpi_rank, n_dim1_loc, n_dim2_loc

        counter_1 = 0
        do i_dim1 = mpi_rank+1, n_dim1_glob - (mpi_size - mpi_rank) +1
            counter_1 = counter_1 + 1
            idx_list_1(counter_1) = i_dim1
            
            counter_2 = 0
            do i_dim2 = mpi_rank+1, n_dim2_glob - (mpi_size - mpi_rank) +1
                counter_2 = counter_2 + 1
                idx_list_2(counter_2) = i_dim2

                loc_mat(counter_1, counter_2) = 0.0d0 + mpi_rank
            end do
        end do

    end subroutine build_local_matrix


    subroutine get_end_datalayout(n_dim1_glob, n_dim2_glob, n_dim1_loc, n_dim2_loc, &
                                  idx_list_1, idx_list_2)
        integer, intent(in) :: n_dim1_glob, n_dim2_glob
        integer, intent(out) :: n_dim1_loc, n_dim2_loc
        integer, dimension(:), allocatable, intent(out) :: idx_list_1
        integer, dimension(:), allocatable, intent(out) :: idx_list_2

        ! internal variables
        integer :: i_dim1, i_dim2, counter_1

        ! 1D cyclic distribution

        n_dim1_loc = 0
        do i_dim1 = 1, n_dim1_glob
            if (mod(i_dim1, mpi_size) /= mpi_rank) cycle
            n_dim1_loc = n_dim1_loc + 1
        end do
        n_dim2_loc = n_dim2_glob

        allocate(idx_list_1(n_dim1_loc))
        allocate(idx_list_2(n_dim2_loc))

        counter_1 = 0
        do i_dim1 = 1, n_dim1_glob
            if (mod(i_dim1, mpi_size) /= mpi_rank) cycle
            counter_1 = counter_1 + 1
            idx_list_1(counter_1) = i_dim1
        end do

        do i_dim2 = 1, n_dim2_glob
            idx_list_2(i_dim2) = i_dim2
        end do

    end subroutine get_end_datalayout


end program test

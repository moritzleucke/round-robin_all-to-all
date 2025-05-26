module distributed_sum

    use mpi, only: MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_ALLREDUCE, MPI_MAX, MPI_STATUS_SIZE, MPI_Sendrecv
    use mpi_env, only: mpi_rank, mpi_size
    use round_robin, only: get_round_robin_turnament_schedule

    implicit none

    private
    public :: sum_and_redistribute

    type :: idx_list_package
        integer :: n_dim_1
        integer :: n_dim_2
        integer, dimension(:), pointer :: idx_list_1 => null()
        integer, dimension(:), pointer :: idx_list_2 => null()
        integer, dimension(:), pointer :: package => null()
    end type idx_list_package

    ! max package size in real(kind=8) bit units
    integer, parameter :: MAX_PACKAGE_LEN = 1000000

contains

    subroutine sum_and_redistribute(idx_list_1, idx_list_2, loc_mat, end_idx_list_1, end_idx_list_2, end_mat)
        integer, dimension(:), intent(in) :: idx_list_1
        integer, dimension(:), intent(in) :: idx_list_2
        real(kind=8), dimension(:, :), intent(in) :: loc_mat
        integer, dimension(:), intent(in) :: end_idx_list_1
        integer, dimension(:), intent(in) :: end_idx_list_2
        real(kind=8), dimension(:, :), allocatable, intent(out) :: end_mat

        ! internal variables
        integer :: loc_n_dim_1, loc_n_dim_2
        integer :: end_n_dim_1, end_n_dim_2
        integer :: n_comm, i_comm
        integer :: i_error
        integer :: max_end_package_len, my_max
        integer :: end_package_len
        integer :: len_offering, max_oponent_offering
        integer :: oponent
        integer, dimension(:), allocatable, target :: end_package
        integer, dimension(:), allocatable, target :: my_offering
        integer, dimension(:), allocatable, target :: recv_end_package
        integer, dimension(:), allocatable, target :: oponent_offering
        integer, dimension(:), allocatable :: comm_schedule
        integer, dimension(MPI_STATUS_SIZE) :: status
        real(kind=8), dimension(MAX_PACKAGE_LEN) :: send_buffer
        real(kind=8), dimension(MAX_PACKAGE_LEN) :: recv_buffer
        type(idx_list_package) :: oponent_end
        type(idx_list_package) :: oponent_offer

        loc_n_dim_1 = size(idx_list_1)
        loc_n_dim_2 = size(idx_list_2)
        end_n_dim_1 = size(end_idx_list_1)
        end_n_dim_2 = size(end_idx_list_2)

        ! consistency check

        ! allocate and initialize the end_mat
        if (allocated(end_mat)) deallocate (end_mat)
        allocate (end_mat(end_n_dim_1, end_n_dim_2))
        end_mat(:, :) = 0.0d0

        ! get communication schedule
        call get_round_robin_turnament_schedule(mpi_size, mpi_rank, n_comm, comm_schedule)

        ! prepare the package for communicating the end matrix layout
        call pack_idx_lists(end_idx_list_1, end_idx_list_2, end_package_len, end_package)

        ! get the maximum idx_list size of end data layout of all processors
        call MPI_Allreduce(end_package_len, max_end_package_len, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, i_error)
        allocate (recv_end_package(max_end_package_len))

        ! offering of oponent idices can only be so large as my end data layout
        allocate(oponent_offering(end_n_dim_1 + end_n_dim_2 + 2))

        ! communicate with every process
        do i_comm = 1, n_comm

            oponent = comm_schedule(i_comm)
            ! cycle if ghost oponent
            if (oponent == -1) cycle

            ! communicate the end data layout
            call MPI_Sendrecv(end_package, end_package_len, MPI_INTEGER, oponent, i_comm, &
                              recv_end_package, max_end_package_len, MPI_INTEGER, oponent, i_comm, &
                              MPI_COMM_WORLD, status, i_error)

            call unpack_idx_lists(recv_end_package, oponent_end)

            ! compute the sendbuffer layout
            call get_sendbuffer_layout(oponent_end, idx_list_1, idx_list_2, len_offering, my_offering)

            ! comminicate the send layout
            call MPI_Sendrecv(my_offering, len_offering, MPI_INTEGER, oponent, i_comm, &
                              oponent_offering, max_oponent_offering, MPI_INTEGER, oponent, i_comm, &
                              MPI_COMM_WORLD, status, i_error)
            call unpack_idx_lists(oponent_offering, oponent_offer)

            ! send buffer in chunks

        end do

        deallocate (comm_schedule)
        deallocate (end_package)
    end subroutine sum_and_redistribute



    subroutine get_sendbuffer_layout(oponent_end, idx_list_1, idx_list_2, len_offering, my_offering)
        type(idx_list_package), intent(in) :: oponent_end
        integer, dimension(:), intent(in) :: idx_list_1
        integer, dimension(:), intent(in) :: idx_list_2
        integer, intent(out) :: len_offering
        integer, dimension(:), allocatable, intent(out) :: my_offering

        ! internal variables
        integer, dimension(:), allocatable :: intersection_dim_1
        integer, dimension(:), allocatable :: intersection_dim_2
        integer :: n_intersect_1, n_intersect_2

        ! get overlapping indices
        call find_list_intersection(oponent_end%n_dim_1, oponent_end%idx_list_1, &
                                    size(idx_list_1), idx_list_1, &
                                    n_intersect_1, intersection_dim_1)
        call find_list_intersection(oponent_end%n_dim_2, oponent_end%idx_list_2, &
                                    size(idx_list_2), idx_list_2, &
                                    n_intersect_2, intersection_dim_2)

        ! pack the send buffer
        call pack_idx_lists(intersection_dim_1, intersection_dim_2, len_offering, my_offering)

    end subroutine get_sendbuffer_layout



    !> @brief find the intersection of two index lists
    !!
    !! @param[in] n_1 -- number of indices in the first list
    !! @param[in] idx_list_1 -- first list of indices
    !! @param[in] n_2 -- number of indices in the second list
    !! @param[in] idx_list_2 -- second list of indices
    !! @param[out] n_intersect -- number of indices in the intersection
    !! @param[out] intersection -- idices that are in both lists
    subroutine find_list_intersection(n_1, idx_list_1, n_2, idx_list_2, n_intersect, intersection)
        integer, intent(in) :: n_1
        integer, dimension(:), intent(in) :: idx_list_1
        integer, intent(in) :: n_2
        integer, dimension(:), intent(in) :: idx_list_2
        integer, intent(out) :: n_intersect
        integer, dimension(:), allocatable, intent(out) :: intersection

        ! internal variables
        integer :: i, j

        ! get the number of overlapping indices
        n_intersect = 0
        do i = 1, n_1
            do j = 1, n_2
                if (idx_list_1(i) == idx_list_2(j)) then
                    n_intersect = n_intersect + 1
                end if
            end do
        end do

        allocate (intersection(n_intersect))

        n_intersect = 0
        do i = 1, n_1
            do j = 1, n_2
                if (idx_list_1(i) == idx_list_2(j)) then
                    n_intersect = n_intersect + 1
                    intersection(n_intersect) = idx_list_1(i)
                end if
            end do
        end do

    end subroutine find_list_intersection



    subroutine pack_idx_lists(idx_list_1, idx_list_2, package_len, package)
        integer, dimension(:), intent(in) :: idx_list_1
        integer, dimension(:), intent(in) :: idx_list_2
        integer, intent(out) :: package_len
        integer, dimension(:), allocatable, intent(out) :: package

        ! internal variables
        integer :: n_dim_1, n_dim_2

        n_dim_1 = size(idx_list_1)
        n_dim_2 = size(idx_list_2)
        package_len = n_dim_1 + n_dim_2 + 2

        ! allocate and initialize the package
        allocate (package(package_len))

        ! write the dimensions in the head
        package(1) = n_dim_1
        package(2) = n_dim_2
        ! copy the idx_list_1 and idx_list_2 into the package
        package(3:n_dim_1 + 2) = idx_list_1(:)
        package(n_dim_1 + 3:package_len) = idx_list_2(:)

    end subroutine pack_idx_lists



    subroutine unpack_idx_lists(package_in, pkg)
        integer, dimension(:), target, intent(in) :: package_in
        type(idx_list_package), intent(out) :: pkg

        ! internal variables
        integer :: n1, n2

        ! Read sizes
        n1 = package_in(1)
        n2 = package_in(2)

        ! Assign metadata
        pkg%n_dim_1 = n1
        pkg%n_dim_2 = n2

        ! Point to full package
        pkg%package => package_in

        ! Point idx_list_1 and idx_list_2 to slices of the package
        pkg%idx_list_1 => pkg%package(3:n1 + 2)
        pkg%idx_list_2 => pkg%package(n1 + 3:n1 + n2 + 2)

    end subroutine unpack_idx_lists

end module distributed_sum

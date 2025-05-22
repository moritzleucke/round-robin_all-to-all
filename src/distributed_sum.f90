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
        integer, dimension(:), pointer :: package    => null()
    end type idx_list_package

contains

    subroutine sum_and_redistribute(idx_list_1, idx_list_2, loc_mat, end_idx_list_1, end_idx_list_2, end_mat)
        integer, dimension(:), intent(in) :: idx_list_1
        integer, dimension(:), intent(in) :: idx_list_2
        real(kind=8), dimension(:,:), intent(in) :: loc_mat
        integer, dimension(:), intent(in) :: end_idx_list_1
        integer, dimension(:), intent(in) :: end_idx_list_2
        real(kind=8), dimension(:,:), allocatable, intent(out) :: end_mat

        
        ! internal variables
        integer :: loc_n_dim_1, loc_n_dim_2
        integer :: end_n_dim_1, end_n_dim_2
        integer :: n_comm, i_comm
        integer :: i_error
        integer :: max_end_package_len, my_max
        integer :: end_package_len
        integer :: oponent 
        integer, dimension(:), allocatable, target :: end_package
        integer, dimension(:), allocatable, target :: recv_end_package
        integer, dimension(:), allocatable :: comm_schedule
        integer, dimension(MPI_STATUS_SIZE) :: status
        type(idx_list_package) :: oponent_end


        loc_n_dim_1 = size(idx_list_1)
        loc_n_dim_2 = size(idx_list_2)
        end_n_dim_1 = size(end_idx_list_1)
        end_n_dim_2 = size(end_idx_list_2)

        ! consistency check

        ! allocate and initialize the end_mat
        if (allocated(end_mat)) deallocate(end_mat)
        allocate(end_mat(end_n_dim_1, end_n_dim_2))
        end_mat(:, :) = 0.0d0

        ! get communication schedule
        call get_round_robin_turnament_schedule(mpi_size, mpi_rank, n_comm, comm_schedule)

        ! prepare the package for communicating the end matrix layout
        call pack_idx_lists(end_idx_list_1, end_idx_list_2, end_package_len, end_package)

        ! get the maximum idx_list size of all processors
        call MPI_Allreduce(end_package_len, max_end_package_len, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, i_error)
        allocate(recv_end_package(max_end_package_len))

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

            ! compute the sendbuffer

            ! comminicate the send layout 

            ! send buffer in chunks 

            

        end do 



        deallocate(comm_schedule)
        deallocate(end_package)
    end subroutine sum_and_redistribute

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
        if (allocated(package)) deallocate(package)
        allocate(package(package_len))

        ! write the dimensions in the head
        package(1) = n_dim_1
        package(2) = n_dim_2
        ! copy the idx_list_1 and idx_list_2 into the package
        package(3:n_dim_1+2) = idx_list_1(:)
        package(n_dim_1+3:package_len) = idx_list_2(:)

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
        pkg%idx_list_1 => pkg%package(3:n1+2)
        pkg%idx_list_2 => pkg%package(n1+3:n1+n2+2)

    end subroutine unpack_idx_lists

end module distributed_sum

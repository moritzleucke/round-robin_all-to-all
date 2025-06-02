module distributed_sum

    use mpi, only: MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_ALLREDUCE, &
        MPI_MAX, MPI_STATUS_SIZE, MPI_Sendrecv, MPI_STATUS_IGNORE, MPI_Get_count, &
        MPI_Send, MPI_Recv
    use mpi_env, only: mpi_rank, mpi_size
    use round_robin, only: get_round_robin_turnament_schedule

    implicit none

    private
    public :: sum_and_redistribute

    !> wrapper for a package containing two index lists
    type :: idx_list_package
        !> number of indices in the first index list
        integer :: n_dim_1
        !> number of indices in the second index list
        integer :: n_dim_2
        !> pointer to the list of global indices in first dimension
        integer, dimension(:), pointer :: idx_list_1 => null()
        !> pointer to the list of global indices in second dimension
        integer, dimension(:), pointer :: idx_list_2 => null()
        !> pointer to the package containing both index lists
        integer, dimension(:), pointer :: package => null()
    end type idx_list_package

    !> max package size in real(kind=8) bit units
    integer, parameter :: MAX_PACKAGE_LEN = 1024000

    !> flag for enabling buffered exchange of matrix elements (memory efficient)
    logical, parameter :: flag_exchange_buffered = .true.

contains

    !> @brief Sums and redistributes the local matrices elements to the end data
    !!        layout among all processes
    !!
    !!        Take the local matrices from all processes and redistribute them
    !!        in the end data layout. Some local matrices of processes may overlap,
    !!        in these cases sum up the matrix elements.
    !!
    !! @param[in] idx_list_1 -- global index list 1 of the local matrix
    !! @param[in] idx_list_2 -- global index list 2 of the local matrix
    !! @param[in] loc_mat -- local matrix
    !! @param[in] end_idx_list_1 -- global index list 1 of the end data layout
    !! @param[in] end_idx_list_2 -- global index list 2 of the end data layout
    !! @param[out] end_mat -- end matrix, where the results are accumulated
    subroutine sum_and_redistribute(idx_list_1, idx_list_2, loc_mat, &
                                    end_idx_list_1, end_idx_list_2, end_mat)
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
        integer :: max_end_package_len
        integer :: end_package_len
        integer :: len_offering, max_oponent_offering
        integer :: oponent
        integer, dimension(:), allocatable, target :: end_package
        integer, dimension(:), allocatable, target :: my_offering
        integer, dimension(:), allocatable, target :: recv_end_package
        integer, dimension(:), allocatable, target :: oponent_offering
        integer, dimension(:), allocatable :: comm_schedule
        integer, dimension(:), allocatable :: inv_idx_list_1, inv_idx_list_2
        integer, dimension(:), allocatable :: inv_end_idx_list_1, inv_end_idx_list_2
        type(idx_list_package) :: oponent_end
        type(idx_list_package) :: oponent_offer
        type(idx_list_package) :: my_offer

        loc_n_dim_1 = size(idx_list_1)
        loc_n_dim_2 = size(idx_list_2)
        end_n_dim_1 = size(end_idx_list_1)
        end_n_dim_2 = size(end_idx_list_2)

        ! get inverse index lists
        call get_inv_idx_list(idx_list_1, inv_idx_list_1)
        call get_inv_idx_list(idx_list_2, inv_idx_list_2)
        call get_inv_idx_list(end_idx_list_1, inv_end_idx_list_1)
        call get_inv_idx_list(end_idx_list_2, inv_end_idx_list_2)

        ! allocate and initialize the end_mat
        call initialize_end_mat(idx_list_1, idx_list_2, loc_mat, &
                                end_n_dim_1, end_n_dim_2, &
                                inv_end_idx_list_1, inv_end_idx_list_2, end_mat)

        ! nothing to do if serial execution
        if (mpi_size == 1) return 

        ! get communication schedule
        call get_round_robin_turnament_schedule(mpi_size, mpi_rank, n_comm, comm_schedule)

        ! prepare the package for communicating the end matrix layout
        call pack_idx_lists(end_idx_list_1, end_idx_list_2, end_package_len, end_package)

        ! get the maximum idx_list size of end data layout of all processors
        call MPI_Allreduce(end_package_len, max_end_package_len, 1, MPI_INTEGER, &
                           MPI_MAX, MPI_COMM_WORLD, i_error)
        allocate(recv_end_package(max_end_package_len))

        ! offering of oponent idices can only be so large as my end data layout
        max_oponent_offering = end_n_dim_1 + end_n_dim_2 + 2
        allocate(oponent_offering(max_oponent_offering))

        ! communicate with every process
        do i_comm = 1, n_comm

            oponent = comm_schedule(i_comm)
            ! skip step if oponent is a ghost
            if (oponent == -1) cycle

            ! communicate the end data layout
            call MPI_Sendrecv(end_package, end_package_len, MPI_INTEGER, oponent, i_comm, &
                              recv_end_package, max_end_package_len, MPI_INTEGER, oponent, i_comm, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, i_error)

            call unpack_idx_lists(recv_end_package, oponent_end)

            ! compute the sendbuffer layout
            call get_sendbuffer_layout(oponent_end, idx_list_1, idx_list_2, &
                                       len_offering, my_offering)

            ! comminicate the send layout
            call MPI_Sendrecv(my_offering, len_offering, MPI_INTEGER, oponent, i_comm, &
                              oponent_offering, max_oponent_offering, MPI_INTEGER, oponent, i_comm, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, i_error)
            call unpack_idx_lists(oponent_offering, oponent_offer)
            call unpack_idx_lists(my_offering, my_offer)

            ! send and recv buffer
            if (flag_exchange_buffered) then
                call exchange_matrix_elements_buffered(oponent, my_offer, oponent_offer, &
                                                       idx_list_1, idx_list_2, &
                                                       inv_idx_list_1, inv_idx_list_2, &
                                                       end_idx_list_1, end_idx_list_2, &
                                                       inv_end_idx_list_1, inv_end_idx_list_2, &
                                                       loc_mat, end_mat)
            else
                call exchange_matrix_elements(oponent, my_offer, oponent_offer, &
                                              idx_list_1, idx_list_2, &
                                              end_idx_list_1, end_idx_list_2, &
                                              loc_mat, end_mat)
            end if
            if (mpi_rank == 0) print *, "step ", i_comm, "/", n_comm, " done"

        end do

        !call print_matrix(end_mat, end_idx_list_1, end_idx_list_2)
        !call print_matrix(loc_mat, idx_list_1, idx_list_2)

        deallocate(comm_schedule)
        deallocate(end_package)
        deallocate(my_offering)
        deallocate(recv_end_package)
        deallocate(oponent_offering)
        deallocate(inv_idx_list_1)
        deallocate(inv_idx_list_2)
        deallocate(inv_end_idx_list_1)
        deallocate(inv_end_idx_list_2)
    end subroutine sum_and_redistribute




    subroutine print_matrix(mat, idx_list_1, idx_list_2)
        real(kind=8), dimension(:, :), intent(in) :: mat
        integer, dimension(:), intent(in) :: idx_list_1
        integer, dimension(:), intent(in) :: idx_list_2
        integer :: i, j, idx_1, idx_2

        do i = 1, size(mat, 1)
            idx_1 = idx_list_1(i)
            do j = 1, size(mat, 2)
                idx_2 = idx_list_2(j)
                print *, idx_1, idx_2, mat(i, j)
            end do
        end do
    end subroutine print_matrix



    !> @brief allocate, initialize and fill the end matrix with the local matrix
    !!
    !! @param[in] idx_list_1 -- global index list 1
    !! @param[in] idx_list_2 -- global index list 2
    !! @param[in] loc_mat -- local matrix
    !! @param[in] end_idx_list_1 -- global index list 1 (end data layout)
    !! @param[in] end_idx_list_2 -- global index list 2 (end data layout)
    !! @param[out] end_mat -- end matrix
    subroutine initialize_end_mat(idx_list_1, idx_list_2, loc_mat, n_end_1, n_end_2, &
                                  inv_end_idx_list_1, inv_end_idx_list_2, end_mat)
        integer, dimension(:), intent(in) :: idx_list_1
        integer, dimension(:), intent(in) :: idx_list_2
        real(kind=8), dimension(:, :), intent(in) :: loc_mat
        integer, intent(in) :: n_end_1
        integer, intent(in) :: n_end_2
        integer, dimension(:), intent(in) :: inv_end_idx_list_1
        integer, dimension(:), intent(in) :: inv_end_idx_list_2
        real(kind=8), dimension(:, :), allocatable, intent(out) :: end_mat

        ! internal variables
        integer :: i, j, idx_1, idx_2, max_end_idx_1, max_end_idx_2

        allocate(end_mat(n_end_1, n_end_2))
        end_mat(:, :) = 0.0d0

        max_end_idx_1 = size(inv_end_idx_list_1)
        max_end_idx_2 = size(inv_end_idx_list_2)

        ! fill the end_mat with the data from the local matrix
        do j = 1, size(idx_list_2)

            if (idx_list_2(j) > max_end_idx_2) cycle
            idx_2 = inv_end_idx_list_2(idx_list_2(j))
            if (idx_2 == -1) cycle

            do i = 1, size(idx_list_1)

                if (idx_list_1(i) > max_end_idx_1) cycle
                idx_1 = inv_end_idx_list_1(idx_list_1(i))
                if (idx_1 == -1) cycle
                end_mat(idx_1, idx_2) = loc_mat(i, j)

            end do
        end do

    end subroutine initialize_end_mat



    !> @brief Exchange matrix elements with the oponent
    !!
    !! @param[in] oponent_rank -- rank of the oponent
    !! @param[in] my_offer -- my offering of indices
    !! @param[in] oponent_offer -- oponent's offering of indices
    !! @param[in] idx_list_1 -- global index list 1
    !! @param[in] idx_list_2 -- global index list 2
    !! @param[in] end_idx_list_1 -- global index list 1 (end data layout)
    !! @param[in] end_idx_list_2 -- global index list 2 (end data layout)
    !! @param[in] loc_mat -- local matrix
    !! @param[inout] end_mat -- end matrix, where the results are accumulated
    subroutine exchange_matrix_elements(oponent_rank, my_offer, oponent_offer, &
                                        idx_list_1, idx_list_2, end_idx_list_1, end_idx_list_2, &
                                        loc_mat, end_mat)
        integer, intent(in) :: oponent_rank
        type(idx_list_package), intent(in) :: my_offer
        type(idx_list_package), intent(in) :: oponent_offer
        integer, dimension(:), intent(in) :: idx_list_1
        integer, dimension(:), intent(in) :: idx_list_2
        integer, dimension(:), intent(in) :: end_idx_list_1
        integer, dimension(:), intent(in) :: end_idx_list_2
        real(kind=8), dimension(:, :), intent(in) :: loc_mat
        real(kind=8), dimension(:, :), intent(inout) :: end_mat

        ! internal variables
        integer :: send_size, recv_size, i, j, idx_1, idx_2, i_error
        real(kind=8), dimension(:,:), allocatable :: send_buffer
        real(kind=8), dimension(:,:), allocatable :: recv_buffer

        ! check if there is anything to exchange
        send_size = my_offer%n_dim_1 * my_offer%n_dim_2
        recv_size = oponent_offer%n_dim_1 * oponent_offer%n_dim_2
        if (send_size == 0 .and. recv_size == 0) return

        ! prepare sendbuffer
        allocate (send_buffer(my_offer%n_dim_1, my_offer%n_dim_2))
        do j = 1, my_offer%n_dim_2
            idx_2 = find_index(idx_list_2, my_offer%idx_list_2(j))
            do i = 1, my_offer%n_dim_1
                idx_1 = find_index(idx_list_1, my_offer%idx_list_1(i))
                send_buffer(i, j) = loc_mat(idx_1, idx_2)
            end do
        end do

        ! prepare recvbuffer
        allocate (recv_buffer(oponent_offer%n_dim_1, oponent_offer%n_dim_2))

        ! exchange the matrices
        call MPI_Sendrecv(send_buffer, send_size, MPI_DOUBLE_PRECISION, oponent_rank, 0, &
                          recv_buffer, recv_size, MPI_DOUBLE_PRECISION, oponent_rank, 0, &
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, i_error)
        deallocate(send_buffer)

        ! unpack the received matrix into the end_mat
        do j = 1, oponent_offer%n_dim_2
            idx_2 = find_index(end_idx_list_2, oponent_offer%idx_list_2(j))
            do i = 1, oponent_offer%n_dim_1
                idx_1 = find_index(end_idx_list_1, oponent_offer%idx_list_1(i))
                end_mat(idx_1, idx_2) = end_mat(idx_1, idx_2) + recv_buffer(i, j)
            end do
        end do

    end subroutine exchange_matrix_elements



    !> @brief Exchange matrix elements with the oponent using buffered communication
    !!
    !!        This is a memory efficient version of exchange_matrix_elements()
    !!
    !! @param[in] oponent_rank -- rank of the oponent
    !! @param[in] my_offer -- my offering of indices
    !! @param[in] oponent_offer -- oponent's offering of indices
    !! @param[in] idx_list_1 -- global index list 1
    !! @param[in] idx_list_2 -- global index list 2
    !! @param[in] end_idx_list_1 -- global index list 1 (end data layout)
    !! @param[in] end_idx_list_2 -- global index list 2 (end data layout)
    !! @param[in] loc_mat -- local matrix
    !! @param[inout] end_mat -- end matrix, where the results are accumulated
    subroutine exchange_matrix_elements_buffered(oponent_rank, my_offer, oponent_offer, &
                                                 idx_list_1, idx_list_2, &
                                                 inv_idx_list_1, inv_idx_list_2, &
                                                 end_idx_list_1, end_idx_list_2, &
                                                 inv_end_idx_list_1, inv_end_idx_list_2, &
                                                 loc_mat, end_mat)
        integer, intent(in) :: oponent_rank
        type(idx_list_package), intent(in) :: my_offer
        type(idx_list_package), intent(in) :: oponent_offer
        integer, dimension(:), intent(in) :: idx_list_1
        integer, dimension(:), intent(in) :: idx_list_2
        integer, dimension(:), intent(in) :: inv_idx_list_1
        integer, dimension(:), intent(in) :: inv_idx_list_2
        integer, dimension(:), intent(in) :: end_idx_list_1
        integer, dimension(:), intent(in) :: end_idx_list_2
        integer, dimension(:), intent(in) :: inv_end_idx_list_1
        integer, dimension(:), intent(in) :: inv_end_idx_list_2
        real(kind=8), dimension(:, :), intent(in) :: loc_mat
        real(kind=8), dimension(:, :), intent(inout) :: end_mat

        ! internal variables
        real(kind=8), dimension(:), allocatable :: send_buffer
        real(kind=8), dimension(:), allocatable :: recv_buffer
        integer :: send_size, recv_size, steps, i_step, i_send, i_recv
        integer :: send_msg_size, recv_msg_size, i_error
        integer, dimension(MPI_STATUS_SIZE) :: status

        send_size = my_offer%n_dim_1 * my_offer%n_dim_2
        recv_size = oponent_offer%n_dim_1 * oponent_offer%n_dim_2

        if (send_size == 0 .and. recv_size == 0) return

        ! round up the division by MAX_PACKAGE_LEN
        steps = (max(send_size, recv_size) + MAX_PACKAGE_LEN - 1)/ MAX_PACKAGE_LEN



        ! prepare the buffers
        allocate(send_buffer(MAX_PACKAGE_LEN))
        allocate(recv_buffer(MAX_PACKAGE_LEN))

        i_send = 0
        i_recv = 0
        do i_step = 1, steps
            if ((i_send < send_size) .and. (i_recv < recv_size)) then
                ! send AND receive

                ! fill the send buffer
                call prep_sendbuffer_low_scaling(i_send, send_msg_size)

                ! send the buffer
                call MPI_Sendrecv(send_buffer, send_msg_size, MPI_DOUBLE_PRECISION, oponent_rank, i_step, &
                                  recv_buffer, MAX_PACKAGE_LEN, MPI_DOUBLE_PRECISION, oponent_rank, i_step, &
                                  MPI_COMM_WORLD, status, i_error)
                call MPI_Get_count(status, MPI_DOUBLE_PRECISION, recv_msg_size, i_error)

                ! unpack the received matrix into the end_mat
                call unpack_recvbuffer_low_scaling(i_recv, recv_msg_size)

                i_send = i_send + send_msg_size
                i_recv = i_recv + recv_msg_size

            else if (i_send < send_size .and. (i_recv == recv_size)) then
                ! send only

                call prep_sendbuffer_low_scaling(i_send, send_msg_size)
                call MPI_Send(send_buffer, send_msg_size, MPI_DOUBLE_PRECISION, oponent_rank, i_step, &
                              MPI_COMM_WORLD, i_error)
                i_send = i_send + send_msg_size

            else if ((i_send == send_size) .and. (i_recv < recv_size)) then
                ! receive only

                call MPI_Recv(recv_buffer, MAX_PACKAGE_LEN, MPI_DOUBLE_PRECISION, oponent_rank, i_step, &
                              MPI_COMM_WORLD, status, i_error)
                call MPI_Get_count(status, MPI_DOUBLE_PRECISION, recv_msg_size, i_error)

                call unpack_recvbuffer_low_scaling(i_recv, recv_msg_size)
                i_recv = i_recv + recv_msg_size

            end if
        end do

        deallocate(send_buffer)
        deallocate(recv_buffer)

    contains

        subroutine prep_sendbuffer(start_idx_cmb, this_msg_size)
            integer , intent(in) :: start_idx_cmb
            integer, intent(out) :: this_msg_size

            ! internal variables
            integer :: i, j,  tmp_cmb_idx
            integer :: idx_1, idx_2

            tmp_cmb_idx = 0
            this_msg_size = 0
            do i = 1, my_offer%n_dim_1

                idx_1 = inv_idx_list_1(my_offer%idx_list_1(i))

                do j = 1, my_offer%n_dim_2
                    tmp_cmb_idx = tmp_cmb_idx + 1
                    if (tmp_cmb_idx > start_idx_cmb) then

                        idx_2 = inv_idx_list_2(my_offer%idx_list_2(j))

                        this_msg_size = this_msg_size + 1

                        send_buffer(this_msg_size) = loc_mat(idx_1, idx_2)

                        if (this_msg_size == MAX_PACKAGE_LEN) return
                    end if
                end do
            end do

        end subroutine prep_sendbuffer

        subroutine prep_sendbuffer_low_scaling(start_idx_cmb, this_msg_size)
            integer, intent(in)  :: start_idx_cmb
            integer, intent(out) :: this_msg_size
        
            ! internal variables
            integer :: total_combinations, offset, linear_idx
            integer :: i, j, idx_1, idx_2
            integer :: dim1, dim2
        
            dim1 = my_offer%n_dim_1
            dim2 = my_offer%n_dim_2
        
            total_combinations = dim1 * dim2
            offset = start_idx_cmb
            this_msg_size = 0
        
            do linear_idx = offset + 1, min(offset + MAX_PACKAGE_LEN, total_combinations)
                i = mod(linear_idx - 1, dim1) + 1
                j = (linear_idx - 1) / dim1 + 1
        
                idx_1 = inv_idx_list_1(my_offer%idx_list_1(i))
                idx_2 = inv_idx_list_2(my_offer%idx_list_2(j))
        
                this_msg_size = this_msg_size + 1
                send_buffer(this_msg_size) = loc_mat(idx_1, idx_2)
            end do
        end subroutine prep_sendbuffer_low_scaling

        subroutine unpack_recvbuffer(start_idx_cmb, this_msg_size)
            integer, intent(in) :: start_idx_cmb
            integer, intent(in) :: this_msg_size

            ! internal variables
            integer :: i, j, tmp_cmb_idx
            integer :: idx_1, idx_2, counter

            tmp_cmb_idx = 0
            counter = 0
            do i = 1, oponent_offer%n_dim_1
                idx_1 = inv_end_idx_list_1(oponent_offer%idx_list_1(i))
                do j = 1, oponent_offer%n_dim_2
                    tmp_cmb_idx = tmp_cmb_idx + 1
                    if ((tmp_cmb_idx > start_idx_cmb) .and. (counter < this_msg_size)) then

                        idx_2 = inv_end_idx_list_2(oponent_offer%idx_list_2(j))

                        counter = counter + 1
                        end_mat(idx_1, idx_2) = end_mat(idx_1, idx_2) + recv_buffer(counter)

                        if (counter == this_msg_size) return
                    end if
                end do
            end do

        end subroutine unpack_recvbuffer

        subroutine unpack_recvbuffer_low_scaling(start_idx_cmb, this_msg_size)
            integer, intent(in) :: start_idx_cmb
            integer, intent(in) :: this_msg_size
        
            ! internal variables
            integer :: total_combinations, offset, linear_idx
            integer :: i, j, idx_1, idx_2
            integer :: dim1, dim2
            integer :: counter
        
            dim1 = oponent_offer%n_dim_1
            dim2 = oponent_offer%n_dim_2
            total_combinations = dim1 * dim2
            offset = start_idx_cmb
            counter = 0
        
            do linear_idx = offset + 1, min(offset + this_msg_size, total_combinations)
                i = mod(linear_idx - 1, dim1) + 1
                j = (linear_idx - 1) / dim1 + 1
        
                idx_1 = inv_end_idx_list_1(oponent_offer%idx_list_1(i))
                idx_2 = inv_end_idx_list_2(oponent_offer%idx_list_2(j))
        
                counter = counter + 1
                end_mat(idx_1, idx_2) = end_mat(idx_1, idx_2) + recv_buffer(counter)
            end do
        end subroutine unpack_recvbuffer_low_scaling

    end subroutine exchange_matrix_elements_buffered


    subroutine get_inv_idx_list(idx_list, inv_idx_list)
        integer, dimension(:), intent(in) :: idx_list 
        integer, dimension(:), allocatable, intent(out) :: inv_idx_list

        ! internal variables
        integer :: i, n, max_idx

        n = size(idx_list)
        max_idx = maxval(idx_list)
        allocate (inv_idx_list(max_idx))

        inv_idx_list(:) = -1 
        do i = 1, n
            if (idx_list(i) < 1 .or. idx_list(i) > max_idx) then
                print *, "Error: index out of bounds in get_inv_idx_list"
                stop
            end if
            inv_idx_list(idx_list(i)) = i
        end do

    end subroutine get_inv_idx_list



    !> @brief compare the indices that my_process has with the indices that oponent
    !!        needs and return the layout of the resulting send buffer
    !!
    !! @param[in] oponent_end -- end data layout of the oponent
    !! @param[in] idx_list_1 -- my global index list 1
    !! @param[in] idx_list_2 -- my global index list 2
    !! @param[out] len_offering -- length of the resulting send buffer
    !! @param[out] my_offering -- resulting send buffer
    subroutine get_sendbuffer_layout(oponent_end, idx_list_1, idx_list_2, &
                                     len_offering, my_offering)
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
        call find_list_intersection_low_scaling(oponent_end%n_dim_1, oponent_end%idx_list_1, &
                                    size(idx_list_1), idx_list_1, &
                                    n_intersect_1, intersection_dim_1)
        call find_list_intersection_low_scaling(oponent_end%n_dim_2, oponent_end%idx_list_2, &
                                    size(idx_list_2), idx_list_2, &
                                    n_intersect_2, intersection_dim_2)

        ! pack the send buffer
        call pack_idx_lists(intersection_dim_1, intersection_dim_2, &
                            len_offering, my_offering)

    end subroutine get_sendbuffer_layout



    !> @brief find the intersection of two index lists
    !!
    !! @param[in] n_1 -- number of indices in the first list
    !! @param[in] idx_list_1 -- first list of indices
    !! @param[in] n_2 -- number of indices in the second list
    !! @param[in] idx_list_2 -- second list of indices
    !! @param[out] n_intersect -- number of indices in the intersection
    !! @param[out] intersection -- idices that are in both lists
    subroutine find_list_intersection(n_1, idx_list_1, n_2, idx_list_2, &
                                      n_intersect, intersection)
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

    

    !> @brief find the intersection of two index lists with low scaling
    !!
    !!        !! the two index lists must be sorted in ascending order !!
    !!
    !! @param[in] n_1 -- number of indices in the first list
    !! @param[in] idx_list_1 -- first list of indices
    !! @param[in] n_2 -- number of indices in the second list
    !! @param[in] idx_list_2 -- second list of indices
    !! @param[out] n_intersect -- number of indices in the intersection
    !! @param[out] intersection -- idices that are in both lists
    subroutine find_list_intersection_low_scaling(n_1, idx_list_1, n_2, idx_list_2, &
                                      n_intersect, intersection)
        integer, intent(in) :: n_1
        integer, dimension(:), intent(in) :: idx_list_1
        integer, intent(in) :: n_2
        integer, dimension(:), intent(in) :: idx_list_2
        integer, intent(out) :: n_intersect
        integer, dimension(:), allocatable, intent(out) :: intersection

        ! internal variables
        integer :: i, j
        integer, dimension(:), allocatable :: tmp_intersect 

        ! consistency check
        if (n_1 == 0 .or. n_2 == 0) then
            n_intersect = 0
            return
        end if
        
        allocate(tmp_intersect(min(n_1, n_2)))
        n_intersect = 0

        i = 1
        j = 1
    
        ! Two-pointer intersection loop
        do while (i <= n_1 .and. j <= n_2)
            if (idx_list_1(i) < idx_list_2(j)) then
                i = i + 1
            else if (idx_list_1(i) > idx_list_2(j)) then
                j = j + 1
            else
                ! Match found
                n_intersect = n_intersect + 1
                tmp_intersect(n_intersect) = idx_list_1(i)
                i = i + 1
                j = j + 1
            end if
        end do

        if (n_intersect == 0) return

        allocate (intersection(n_intersect))
        intersection(:) = tmp_intersect(1:n_intersect)

        deallocate(tmp_intersect)

    end subroutine find_list_intersection_low_scaling



    !> @brief Packs two index lists into a single array for communication
    !!
    !!        layout:
    !!        [
    !!            n_dim_1, n_dim_2, idx_list_1(1), ..., idx_list_1(n_dim_1),
    !!            idx_list_2(1), ..., idx_list_2(n_dim_2)
    !!        ]
    !!
    !! @param[in] idx_list_1 -- first index list
    !! @param[in] idx_list_2 -- second index list
    !! @param[out] package_len -- length of the resulting package
    !! @param[out] package -- resulting package containing both index lists
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



    !> @brief Unpacks the index lists from a index lists package
    !!
    !!        layout of package_in:
    !!        [
    !!            n_dim_1, n_dim_2, idx_list_1(1), ..., idx_list_1(n_dim_1),
    !!            idx_list_2(1), ..., idx_list_2(n_dim_2)
    !!        ]
    !!
    !! @param[in] package_in -- input package containing the index lists
    !! @param[out] pkg -- struct containing pointers to the index lists
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


    !> @brief Finds the index of a value in an array
    !!
    !!        fortran intrinsic findloc() not available in all compilers
    !!
    !! @param array -- array to search in
    !! @param value -- value to find
    integer function find_index(array, value) result(index)
        integer, dimension(:), intent(in) :: array
        integer, intent(in) :: value

        ! internal variables
        integer :: i

        index = -1  ! Default: not found
        do i = 1, size(array)
            if (array(i) == value) then
                index = i
                return
            end if
        end do
    end function find_index

end module distributed_sum

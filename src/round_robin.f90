module round_robin

    implicit none

    private
    public :: get_round_robin_turnament_schedule

contains

    !> @brief for a given rank get the communication schedule to communicate once
    !!        with all other n_ranks in a pair-wise fashion
    !!
    !!        Implemented as a round robin tournament (circle method) with
    !!        n_ranks = number of players. When n_ranks is odd, one ghost player is added
    !!        ghost player has rank index -1 in the schedule
    !!
    !! @param[in] n_ranks -- number of mpi ranks
    !! @param[in] my_rank -- rank of the current process (0-based)
    !! @param[out] n_ops -- number of operations to be performed (either n_ranks-1 or n_ranks)
    !! @param[out] my_schedule -- schedule of communication partners (0-based idxs)
    subroutine get_round_robin_turnament_schedule(n_ranks, my_rank, n_ops, my_schedule)
        integer, intent(in) :: n_ranks
        integer, intent(in) :: my_rank
        integer, intent(out) :: n_ops
        integer, dimension(:), allocatable, intent(out) :: my_schedule

        ! internal variables
        integer :: i_turnament, i_player, halfe_players
        integer, dimension(:, :), allocatable :: turnament_schedule, old_schedule

        ! consistency checks
        if (n_ranks == 0) then 
            n_ops = 0
            allocate(my_schedule(1)) 
            my_schedule(1) = -1
            return
        end if 
        if (my_rank < 0 .or. my_rank > (n_ranks-1)) then
            print *, "Error: my_rank out of bounds"
            stop
        end if

        ! number of turnaments shall be even shall be even
        n_ops = (n_ranks -1)+ mod(n_ranks, 2)
        halfe_players = (n_ops + 1) /2
        if (allocated(my_schedule)) deallocate(my_schedule)
        allocate(my_schedule(n_ops))
        allocate(turnament_schedule(halfe_players, 2))
        allocate(old_schedule(halfe_players, 2))

        ! obtaining the turnament table for the first match
        do i_player = 0, halfe_players-1
            turnament_schedule(i_player+1, 1) = i_player
            turnament_schedule(i_player+1, 2) = 2*halfe_players - i_player -1
        end do
        my_schedule(1) = find_oponent(my_rank, n_ops, turnament_schedule)

        ! play all turnaments
        do i_turnament = 2, n_ops
            old_schedule(:, :) = turnament_schedule(:, :)

            ! turnament_schedule(1,1) is fixed
            ! rotate other players
            turnament_schedule(2, 1) = old_schedule(1, 2)
            turnament_schedule(halfe_players, 2) = old_schedule(halfe_players, 1)
            turnament_schedule(1, 2) = old_schedule(2, 2)
            do i_player = 3, halfe_players
                turnament_schedule(i_player, 1) = old_schedule(i_player - 1, 1)
                turnament_schedule(i_player-1, 2) = old_schedule(i_player, 2)
            end do

            my_schedule(i_turnament) = find_oponent(my_rank, halfe_players, turnament_schedule)
        end do

        ! give ghost player an index of -1
        if (mod(n_ranks, 2) /= 0) then
            do i_turnament = 1, n_ops
                if (my_schedule(i_turnament) > (n_ranks-1)) then
                    my_schedule(i_turnament) = -1
                end if
            end do
        end if

        deallocate(turnament_schedule)
        deallocate(old_schedule)

    end subroutine get_round_robin_turnament_schedule



    !> @brief find the oponent of the given player in the turnament schedule
    !!
    !! @param[in] this_player -- rank of the current process (0-based)
    !! @param[in] halfe_players -- number of players in the turnament devided by two
    !! @param[in] schedule -- turnament table 
    !! @return oponent -- rank of the oponent (0-based)
    integer function find_oponent(this_player, halfe_players, schedule) result(oponent)
        integer, intent(in) :: this_player
        integer, intent(in) :: halfe_players
        integer, dimension(:, :), intent(in) :: schedule

        ! internal variables
        integer :: i_oponent

        do i_oponent = 1, halfe_players
            if (schedule(i_oponent, 1) == this_player) then
                oponent = schedule(i_oponent, 2)
                return
            else if (schedule(i_oponent, 2) == this_player) then
                oponent = schedule(i_oponent, 1)
                return
            end if
        end do
        print *, "Error: Oponent not found for player ", this_player
        stop
    end function find_oponent

end module round_robin


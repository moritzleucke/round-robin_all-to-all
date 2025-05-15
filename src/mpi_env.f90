module mpi_env
    use mpi
    implicit none
    private
    public :: init_mpi, finalize_mpi, mpi_rank, mpi_size

    integer :: mpi_rank = -1
    integer :: mpi_size = -1

contains

    subroutine init_mpi()
        integer :: ierr

        call MPI_Init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)
    end subroutine init_mpi

    subroutine finalize_mpi()
        integer :: ierr
        call MPI_Finalize(ierr)
    end subroutine finalize_mpi

end module mpi_env
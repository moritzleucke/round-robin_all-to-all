cmake_minimum_required(VERSION 3.15)
project(FortranMPIProject LANGUAGES Fortran)

# Enable MPI
find_package(MPI REQUIRED)

# Set source files
set(SOURCES
    src/main.f90
    src/mpi_env.f90
    src/round_robin.f90
    src/distributed_sum.f90
)

# Create the executable
add_executable(round_rob ${SOURCES})

# Link MPI Fortran
target_link_libraries(round_rob PRIVATE MPI::MPI_Fortran)

# Include MPI module path if needed
target_include_directories(round_rob PRIVATE ${MPI_Fortran_INCLUDE_DIRS})

# Optional: Add compile flags (e.g., for debugging)
target_compile_options(round_rob PRIVATE -O2 -g -fcheck=bounds,do -fsanitize=address,undefined -Wall)
target_link_options(round_rob PRIVATE -fsanitize=address,undefined)
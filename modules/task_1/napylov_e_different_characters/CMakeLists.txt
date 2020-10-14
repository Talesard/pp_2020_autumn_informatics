cmake_minimum_required(VERSION 3.9)
project(ITK_MPI)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

add_executable(main main.cpp count_diff_char_mpi.cpp)
target_link_libraries(main ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(main PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(main PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
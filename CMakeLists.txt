cmake_minimum_required( VERSION 3.14 )
project( mobbrmsd LANGUAGES C Fortran )

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR
    "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
    )
endif()

cmake_policy( SET CMP0076 NEW )
set( CMAKE_POSITION_INDEPENDENT_CODE ON )
add_subdirectory( src )


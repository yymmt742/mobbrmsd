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

find_package(OpenMP QUIET)

#set(BLA_STATIC ON)

if( USE_EXTERNAL_LAPACK )
  #foreach(var Intel10_64ilp_seq OpenBLAS)
  foreach(var OpenBLAS)
    set(BLA_VENDOR ${var})
    find_package(BLAS QUIET)
    find_package(LAPACK QUIET)
    if( BLAS_FOUND AND LAPACK_FOUND )
      message(STATUS "Found BLAS and LAPACK with BLA_VENDOR=${BLA_VENDOR}")
      break()
    else()
      message(STATUS "Not found BLAS and LAPACK with BLA_VENDOR=${BLA_VENDOR}")
    endif( BLAS_FOUND AND LAPACK_FOUND )
    set(BLA_VENDOR)
    set(BLAS_FOUND)
    set(LAPACK_FOUND)
  endforeach()
endif()

set( mobbrmsd_PRECISIONS sp dp )
set( mobbrmsd_DIMENTIONS 2d 3d xd )

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  enable_testing()
  include( CTest )
  find_package(unittest QUIET)
  if(NOT unittest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      unittest
      GIT_REPOSITORY https://github.com/yymmt742/unittest
      GIT_TAG origin/main
      GIT_PROGRESS TRUE
      )
    set(CMAKE_MAXIMUM_RANK 2)
    set(CMAKE_BUILD_TYPE Rerease)
    FetchContent_MakeAvailable(unittest)
    set(CMAKE_BUILD_TYPE Debug)
    set(unittest_INCLUDE_DIRS ${unittest_BINARY_DIR}/src/mod_files)
  else()
    set(unittest_INCLUDE_DIRS ${unittest_DIR}/mod_files)
  endif()
endif()

add_subdirectory( src )


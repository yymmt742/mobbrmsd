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

find_package(OpenMP REQUIRED QUIET)
set(BLA_STATIC ON)
find_package(BLAS REQUIRED)

set( mobbrmsd_PRECISIONS sp dp )
set( mobbrmsd_DIMENTIONS 2d 3d )
if( BLAS_FOUND )
  list( APPEND mobbrmsd_DIMENTIONS xd )
endif( BLAS_FOUND )

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
    FetchContent_MakeAvailable(unittest)
    set(unittest_INCLUDE_DIRS ${unittest_BINARY_DIR}/src/mod_files)
  else()
    set(unittest_INCLUDE_DIRS ${unittest_DIR}/mod_files)
  endif()
endif()

add_subdirectory( src )


cmake_minimum_required( VERSION 3.14 )
project( symRMSD )

cmake_policy( SET CMP0076 NEW )

enable_testing()
include( CTest )

include( CMakePrintHelpers )
enable_language( C )
enable_language( Fortran )

set( CMAKE_POSITION_INDEPENDENT_CODE ON )

add_subdirectory( src )


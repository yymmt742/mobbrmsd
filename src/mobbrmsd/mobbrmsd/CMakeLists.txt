set(LIB_MOD_GD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
if(NOT EXISTS "${LIB_MOD_GD_DIR}")
  make_directory("${LIB_MOD_GD_DIR}")
endif()

set(TARGET_GD ${PROJECT_NAME})
add_library(${TARGET_GD} OBJECT)

if ( REAL )
  message ( STATUS "Configuring build for ${REAL}-bit reals" )
  target_compile_definitions(${TARGET_GD} PUBLIC REAL${REAL})
endif ()

if ( INT )
  message ( STATUS "Configuring build for ${INT}-bit integers" )
  target_compile_definitions(${TARGET_GD} PUBLIC INT${INT})
endif ()

target_include_directories(${TARGET_GD} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_sources(${TARGET_GD} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_f_matrix.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_block.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_list.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90 )
target_link_libraries( ${TARGET_GD} PRIVATE blas )

set_target_properties(${TARGET_GD} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_GD_DIR})
add_library(${TARGET_GD}_shared SHARED $<TARGET_OBJECTS:mobbrmsd>)
set_target_properties(${TARGET_GD}_shared PROPERTIES OUTPUT_NAME ${TARGET_GD})
add_library(${TARGET_GD}_static STATIC $<TARGET_OBJECTS:mobbrmsd>)
set_target_properties(${TARGET_GD}_static PROPERTIES OUTPUT_NAME ${TARGET_GD})
target_link_libraries( ${TARGET_GD}_static PRIVATE blas )

install(TARGETS ${TARGET_GD}_static
        EXPORT ${TARGET_GD}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_GD}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_GD}/lib")
install(TARGETS ${TARGET_GD}_shared
        EXPORT ${TARGET_GD}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_GD}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_GD}/lib")
install(EXPORT ${TARGET_GD}-export
        NAMESPACE ${TARGET_GD}::
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_GD}")
install(EXPORT ${TARGET_GD}-export
        FILE ${TARGET_GD}-config.cmake
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_GD}")

install(DIRECTORY ${LIB_MOD_DIR} DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_GD}")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_GD}")

if( CMAKE_BUILD_TYPE MATCHES "Debug" )

  find_package( purelapack REQUIRED )
  find_package( unittest REQUIRED )

  add_executable( test_rotation )
  target_sources( test_rotation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_rotation.f90 )
  target_include_directories( test_rotation PRIVATE ${unittest_DIR}/mod_files )
  target_include_directories( test_rotation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_link_libraries( test_rotation PRIVATE unittest )
  target_link_libraries( test_rotation PRIVATE purelapack )
  add_test( NAME test_rotation COMMAND test_rotation )

  add_executable( test_c_matrix )
  target_sources( test_c_matrix PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_c_matrix.f90 )
  target_include_directories( test_c_matrix PRIVATE ${unittest_DIR}/mod_files )
  target_include_directories( test_c_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_link_libraries( test_c_matrix PRIVATE unittest )
  target_link_libraries( test_c_matrix PRIVATE purelapack )
  add_test( NAME test_c_matrix COMMAND test_c_matrix )

  add_executable( test_f_matrix )
  target_sources( test_f_matrix PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_f_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_f_matrix.f90 )
  target_include_directories( test_f_matrix PRIVATE ${unittest_DIR}/mod_files )
  target_include_directories( test_f_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_link_libraries( test_f_matrix PRIVATE unittest )
  target_link_libraries( test_f_matrix PRIVATE purelapack )
  add_test( NAME test_f_matrix COMMAND test_f_matrix )

  add_executable( test_bb_block )
  target_sources( test_bb_block PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_f_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_bb_block.f90 )
  target_include_directories( test_bb_block PRIVATE ${unittest_DIR}/mod_files )
  target_include_directories( test_bb_block PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_link_libraries( test_bb_block PRIVATE unittest )
  target_link_libraries( test_bb_block PRIVATE purelapack )
  add_test( NAME test_bb_block COMMAND test_bb_block )

  add_executable( test_bb_list )
  target_sources( test_bb_list PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_f_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_list.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_bb_list.f90 )
  target_include_directories( test_bb_list PRIVATE ${unittest_DIR}/mod_files )
  target_include_directories( test_bb_list PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_link_libraries( test_bb_list PRIVATE unittest )
  target_link_libraries( test_bb_list PRIVATE purelapack )
  add_test( NAME test_bb_list COMMAND test_bb_list 4 )

  add_executable( test_mobbrmsd )
  target_sources( test_mobbrmsd PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_f_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_list.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_mobbrmsd.f90 )
  target_include_directories( test_mobbrmsd PRIVATE ${unittest_DIR}/mod_files )
  target_include_directories( test_mobbrmsd PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_link_libraries( test_mobbrmsd PRIVATE unittest )
  target_link_libraries( test_mobbrmsd PRIVATE purelapack )
  add_test( NAME test_mobbrmsd COMMAND test_mobbrmsd 4 )
endif()


set(LIB_MOD_3D_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
if(NOT EXISTS "${LIB_MOD_3D_DIR}")
  make_directory("${LIB_MOD_3D_DIR}")
endif()

set(TARGET_3D ${PROJECT_NAME}_3d)
add_library(${TARGET_3D} OBJECT)

if ( REAL )
  message ( STATUS "Configuring build for ${REAL}-bit reals" )
  target_compile_definitions(${TARGET_3D} PUBLIC REAL${REAL})
endif ()

if ( INT )
  message ( STATUS "Configuring build for ${INT}-bit integers" )
  target_compile_definitions(${TARGET_3D} PUBLIC INT${INT})
endif ()

target_include_directories(${TARGET_3D} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_sources(${TARGET_3D} PRIVATE
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

set_target_properties(${TARGET_3D} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_3D_DIR})
add_library(${TARGET_3D}_shared SHARED $<TARGET_OBJECTS:mobbrmsd_3d>)
set_target_properties(${TARGET_3D}_shared PROPERTIES OUTPUT_NAME ${TARGET_3D})
add_library(${TARGET_3D}_static STATIC $<TARGET_OBJECTS:mobbrmsd_3d>)
set_target_properties(${TARGET_3D}_static PROPERTIES OUTPUT_NAME ${TARGET_3D})

install(TARGETS ${TARGET_3D}_static
        EXPORT ${TARGET_3D}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}/lib")
install(TARGETS ${TARGET_3D}_shared
        EXPORT ${TARGET_3D}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}/lib")
install(EXPORT ${TARGET_3D}-export
        NAMESPACE ${TARGET_3D}::
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")
install(EXPORT ${TARGET_3D}-export
        FILE ${TARGET_3D}-config.cmake
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")

install(DIRECTORY ${LIB_MOD_DIR} DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  find_package( unittest REQUIRED )
  find_package( purelapack REQUIRED )
  add_executable( test_gemm_3d )
  target_sources( test_gemm_3d PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_gemm.f90 )
  target_include_directories( test_gemm_3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_include_directories( test_gemm_3d PRIVATE ${unittest_DIR}/mod_files )
  target_link_libraries( test_gemm_3d PRIVATE unittest )
  add_test( NAME test_gemm_3d COMMAND test_gemm_3d )

  add_executable( test_rotation_3d )
  target_sources( test_rotation_3d PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_rotation.f90 )
  target_include_directories( test_rotation_3d PRIVATE ${unittest_DIR}/mod_files )
  target_include_directories( test_rotation_3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_link_libraries( test_rotation_3d PRIVATE unittest )
  target_link_libraries( test_rotation_3d PRIVATE purelapack )
  add_test( NAME test_rotation_3d COMMAND test_rotation_3d )

  add_executable( test_c_matrix_3d )
  target_sources( test_c_matrix_3d PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_c_matrix.f90 )
  target_include_directories( test_c_matrix_3d PRIVATE ${unittest_DIR}/mod_files )
  target_include_directories( test_c_matrix_3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_link_libraries( test_c_matrix_3d PRIVATE unittest )
  add_test( NAME test_c_matrix_3d COMMAND test_c_matrix_3d )

  add_executable( test_bb_block_3d )
  target_sources( test_bb_block_3d PRIVATE
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
  target_include_directories( test_bb_block_3d PRIVATE ${unittest_DIR}/mod_files )
  target_include_directories( test_bb_block_3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  target_link_libraries( test_bb_block_3d PRIVATE unittest )
  add_test( NAME test_bb_block_3d COMMAND test_bb_block_3d )

  add_executable( test_bb_list_3d )
  target_sources( test_bb_list_3d PRIVATE
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
  target_include_directories( test_bb_list_3d PRIVATE ${unittest_DIR}/mod_files )
  target_link_libraries( test_bb_list_3d PRIVATE unittest )
  target_include_directories( test_bb_list_3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
  add_test( NAME test_bb_list_3d COMMAND test_bb_list_3d 0 )

endif()

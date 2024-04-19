set(LIB_MOD_2D_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
if(NOT EXISTS "${LIB_MOD_2D_DIR}")
  make_directory("${LIB_MOD_2D_DIR}")
endif()

set(TARGET_2D ${PROJECT_NAME}_2d)
add_library(${TARGET_2D} OBJECT)

if(REAL)
  message(STATUS "Configuring build for ${REAL}-bit reals")
  target_compile_definitions(${TARGET_2D} PUBLIC REAL${REAL})
endif()

if(INT)
  message(STATUS "Configuring build for ${INT}-bit integers")
  target_compile_definitions(${TARGET_2D} PUBLIC INT${INT})
endif()

target_sources(
  ${TARGET_2D}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_f_matrix.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_block.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_list.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_header.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_state.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90)

target_include_directories(${TARGET_2D}
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_include_directories(${TARGET_2D}
                           PRIVATE ${forbar_BINARY_DIR}/src/include)
target_link_libraries(${TARGET_2D} PRIVATE ${OpenMP_Fortran_LIB_NAMES})
target_link_libraries(${TARGET_2D} PRIVATE forbar)
set_target_properties(${TARGET_2D} PROPERTIES Fortran_MODULE_DIRECTORY
                                              ${LIB_MOD_2D_DIR})

add_library(${TARGET_2D}_shared SHARED $<TARGET_OBJECTS:mobbrmsd_2d>)
set_target_properties(${TARGET_2D}_shared PROPERTIES OUTPUT_NAME ${TARGET_2D})
target_link_libraries(${TARGET_2D}_shared PRIVATE ${OpenMP_Fortran_LIB_NAMES})
target_link_libraries(${TARGET_2D}_shared PRIVATE forbar)

add_library(${TARGET_2D}_static STATIC $<TARGET_OBJECTS:mobbrmsd_2d>)
set_target_properties(${TARGET_2D}_static PROPERTIES OUTPUT_NAME ${TARGET_2D})
target_link_libraries(${TARGET_2D}_static PRIVATE ${OpenMP_Fortran_LIB_NAMES})
target_link_libraries(${TARGET_2D}_static PRIVATE forbar)

install(
  TARGETS ${TARGET_2D}_static
  EXPORT ${TARGET_2D}-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}/lib")
install(
  TARGETS ${TARGET_2D}_shared
  EXPORT ${TARGET_2D}-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}/lib")
install(
  EXPORT ${TARGET_2D}-export
  NAMESPACE ${TARGET_2D}::
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}")
install(
  EXPORT ${TARGET_2D}-export
  FILE ${TARGET_2D}-config.cmake
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}")

install(DIRECTORY ${LIB_MOD_DIR}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}")

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  find_package(unittest REQUIRED)
  find_package(purelapack REQUIRED)

  add_executable(test_gemm_2d)
  target_sources(
    test_gemm_2d
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_gemm.f90)
  target_include_directories(test_gemm_2d
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  target_include_directories(test_gemm_2d PRIVATE ${unittest_DIR}/mod_files)
  target_link_libraries(test_gemm_2d PRIVATE unittest)
  add_test(NAME test_gemm_2d COMMAND test_gemm_2d)

  add_executable(test_rotation_2d)
  target_sources(
    test_rotation_2d
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_rotation.f90)
  target_include_directories(test_rotation_2d PRIVATE ${unittest_DIR}/mod_files)
  target_include_directories(test_rotation_2d
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  target_link_libraries(test_rotation_2d PRIVATE unittest)
  target_link_libraries(test_rotation_2d PRIVATE purelapack)
  add_test(NAME test_rotation_2d COMMAND test_rotation_2d)

  add_executable(test_c_matrix_2d)
  target_sources(
    test_c_matrix_2d
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_c_matrix.f90)
  target_include_directories(test_c_matrix_2d PRIVATE ${unittest_DIR}/mod_files)
  target_include_directories(test_c_matrix_2d
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  target_link_libraries(test_c_matrix_2d PRIVATE unittest)
  add_test(NAME test_c_matrix_2d COMMAND test_c_matrix_2d)

  add_executable(test_f_matrix_2d)
  target_sources(
    test_f_matrix_2d
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_f_matrix.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_f_matrix.f90)
  target_include_directories(test_f_matrix_2d PRIVATE ${unittest_DIR}/mod_files)
  target_include_directories(test_f_matrix_2d
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  target_link_libraries(test_f_matrix_2d PRIVATE unittest)
  add_test(NAME test_f_matrix_2d COMMAND test_f_matrix_2d)

  add_executable(test_bb_block_2d)
  target_sources(
    test_bb_block_2d
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
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
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_bb_block.f90)
  target_include_directories(test_bb_block_2d PRIVATE ${unittest_DIR}/mod_files)
  target_include_directories(test_bb_block_2d
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  target_link_libraries(test_bb_block_2d PRIVATE unittest)
  add_test(NAME test_bb_block_2d COMMAND test_bb_block_2d)

  add_executable(test_bb_list_2d)
  target_sources(
    test_bb_list_2d
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
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
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_bb_list.f90)
  target_include_directories(test_bb_list_2d PRIVATE ${unittest_DIR}/mod_files)
  target_link_libraries(test_bb_list_2d PRIVATE unittest)
  target_include_directories(test_bb_list_2d
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  add_test(NAME test_bb_list_2d COMMAND test_bb_list_2d 0)

  add_executable(test_mobbrmsd_2d)
  target_sources(
    test_mobbrmsd_2d
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/blas_lapack_interface.f90
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
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_header.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_state.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/mod_testutil.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_mobbrmsd.f90)
  target_include_directories(test_mobbrmsd_2d PRIVATE ${unittest_DIR}/mod_files)
  target_include_directories(test_mobbrmsd_2d
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  target_include_directories(test_mobbrmsd_2d
                             PRIVATE ${forbar_BINARY_DIR}/src/include)
  target_link_libraries(test_mobbrmsd_2d PRIVATE unittest)
  target_link_libraries(test_mobbrmsd_2d PRIVATE purelapack)
  target_link_libraries(test_mobbrmsd_2d PRIVATE forbar)
  target_link_libraries(test_mobbrmsd_2d PRIVATE ${OpenMP_Fortran_LIB_NAMES})
  add_test(NAME test_mobbrmsd_2d COMMAND test_mobbrmsd_2d 0)
endif()

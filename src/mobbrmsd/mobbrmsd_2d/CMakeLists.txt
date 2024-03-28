set(LIB_MOD_2D_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
if(NOT EXISTS "${LIB_MOD_2D_DIR}")
  make_directory("${LIB_MOD_2D_DIR}")
endif()

set(TARGET_2D ${PROJECT_NAME}_2d)
add_library(${TARGET_2D} OBJECT)
target_include_directories(${TARGET_2D} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_sources(${TARGET_2D} PRIVATE
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
  gemm.f90
  blas_lapack_interface.f90
  mod_rotation.f90 )

set_target_properties(${TARGET_2D} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_2D_DIR})
add_library(${TARGET_2D}_shared SHARED $<TARGET_OBJECTS:mobbrmsd_2d>)
set_target_properties(${TARGET_2D}_shared PROPERTIES OUTPUT_NAME ${TARGET_2D})
add_library(${TARGET_2D}_static STATIC $<TARGET_OBJECTS:mobbrmsd_2d>)
set_target_properties(${TARGET_2D}_static PROPERTIES OUTPUT_NAME ${TARGET_2D})

install(TARGETS ${TARGET_2D}_static
        EXPORT ${TARGET_2D}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}/lib")
install(TARGETS ${TARGET_2D}_shared
        EXPORT ${TARGET_2D}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install(EXPORT ${TARGET_2D}-export
        NAMESPACE ${TARGET_2D}::
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}")
install(EXPORT ${TARGET_2D}-export
        FILE ${TARGET_2D}-config.cmake
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}")

install(DIRECTORY ${LIB_MOD_DIR} DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_2D}")

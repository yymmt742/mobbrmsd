set(LIB_MOD_3D_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
if(NOT EXISTS "${LIB_MOD_3D_DIR}")
  make_directory("${LIB_MOD_3D_DIR}")
endif()

set(TARGET_3D ${PROJECT_NAME}_3d)
add_library(${TARGET_3D} OBJECT)
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
  gemm.f90
  blas_lapack_interface.f90
  mod_rotation.f90 )

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
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install(EXPORT ${TARGET_3D}-export
        NAMESPACE ${TARGET_3D}::
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")
install(EXPORT ${TARGET_3D}-export
        FILE ${TARGET_3D}-config.cmake
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")

install(DIRECTORY ${LIB_MOD_DIR} DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")

set(LIB_MOD_XDDP_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_xddp/mod_files)
if(NOT EXISTS "${LIB_MOD_XDDP_DIR}")
  make_directory("${LIB_MOD_XDDP_DIR}")
endif()
set(LIB_MOD_XDSP_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_xdsp/mod_files)
if(NOT EXISTS "${LIB_MOD_XDSP_DIR}")
  make_directory("${LIB_MOD_XDSP_DIR}")
endif()

list(APPEND sources
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
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_state.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_batch_run.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_mst.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_dimspec_functions.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
  )

add_library(${PROJECT_NAME}_xdsp OBJECT ${sources})
add_library(${PROJECT_NAME}_xddp OBJECT ${sources})
set_target_properties(${PROJECT_NAME}_xdsp PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_XDSP_DIR})
set_target_properties(${PROJECT_NAME}_xddp PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_XDDP_DIR})
ADDOPTS(${PROJECT_NAME}_xdsp 2 32 32 "" 1)
ADDOPTS(${PROJECT_NAME}_xddp 2 64 32 "" 1)
add_library(${PROJECT_NAME}_xdsp_static STATIC $<TARGET_OBJECTS:${PROJECT_NAME}_xdsp>)
add_library(${PROJECT_NAME}_xdsp_shared SHARED $<TARGET_OBJECTS:${PROJECT_NAME}_xdsp>)
add_library(${PROJECT_NAME}_xddp_static STATIC $<TARGET_OBJECTS:${PROJECT_NAME}_xddp>)
add_library(${PROJECT_NAME}_xddp_shared SHARED $<TARGET_OBJECTS:${PROJECT_NAME}_xddp>)
set_target_properties(${PROJECT_NAME}_xddp_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_xddp)
set_target_properties(${PROJECT_NAME}_xddp_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_xddp)
set_target_properties(${PROJECT_NAME}_xdsp_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_xdsp)
set_target_properties(${PROJECT_NAME}_xdsp_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_xdsp)
install(
  TARGETS ${PROJECT_NAME}_xddp_shared
  EXPORT ${PROJECT_NAME}_xd-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_xddp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_xddp/lib"
  )
install(
  TARGETS ${PROJECT_NAME}_xddp_static
  EXPORT ${PROJECT_NAME}_xd-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_xddp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_xddp/lib"
  )
install(
  TARGETS ${PROJECT_NAME}_xdsp_shared
  EXPORT ${PROJECT_NAME}_xd-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_xdsp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_xdsp/lib"
  )
install(
  TARGETS ${PROJECT_NAME}_xdsp_static
  EXPORT ${PROJECT_NAME}_xd-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_xdsp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_xdsp/lib"
  )
install(
  EXPORT ${PROJECT_NAME}_xd-export
  NAMESPACE ${PROJECT_NAME}_xd::
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")
install(DIRECTORY ${LIB_MOD_XDDP_DIR}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")
install(DIRECTORY ${LIB_MOD_XDSP_DIR}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")

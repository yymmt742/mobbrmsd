set(LIB_MOD_GDDP_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_gddp/mod_files)
if(NOT EXISTS "${LIB_MOD_GDDP_DIR}")
  make_directory("${LIB_MOD_GDDP_DIR}")
endif()
set(LIB_MOD_GDSP_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_gdsp/mod_files)
if(NOT EXISTS "${LIB_MOD_GDSP_DIR}")
  make_directory("${LIB_MOD_GDSP_DIR}")
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
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_header.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_state.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_batch_run.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_mst.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_dimspec_functions.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_rotation.f90
  )

add_library(${PROJECT_NAME}_gddp_shared SHARED ${sources})
add_library(${PROJECT_NAME}_gddp_static STATIC ${sources})
add_library(${PROJECT_NAME}_gdsp_shared SHARED ${sources})
add_library(${PROJECT_NAME}_gdsp_static STATIC ${sources})
ADDOPTS(${PROJECT_NAME}_gddp_shared 64 32 1 1)
ADDOPTS(${PROJECT_NAME}_gddp_static 64 32 1 1)
ADDOPTS(${PROJECT_NAME}_gdsp_shared 32 32 1 1)
ADDOPTS(${PROJECT_NAME}_gdsp_static 32 32 1 1)
set_target_properties(${PROJECT_NAME}_gddp_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_gddp)
set_target_properties(${PROJECT_NAME}_gddp_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_gddp)
set_target_properties(${PROJECT_NAME}_gdsp_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_gdsp)
set_target_properties(${PROJECT_NAME}_gdsp_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_gdsp)
set_target_properties(${PROJECT_NAME}_gddp_shared PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_GDDP_DIR})
set_target_properties(${PROJECT_NAME}_gddp_static PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_GDDP_DIR})
set_target_properties(${PROJECT_NAME}_gdsp_shared PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_GDSP_DIR})
set_target_properties(${PROJECT_NAME}_gdsp_static PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_GDSP_DIR})
install(
  TARGETS ${PROJECT_NAME}_gddp_shared
  EXPORT ${PROJECT_NAME}_gd-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_gddp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_gddp/lib"
  )
install(
  TARGETS ${PROJECT_NAME}_gddp_static
  EXPORT ${PROJECT_NAME}_gd-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_gddp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_gddp/lib"
  )
install(
  TARGETS ${PROJECT_NAME}_gdsp_shared
  EXPORT ${PROJECT_NAME}_gd-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_gdsp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_gdsp/lib"
  )
install(
  TARGETS ${PROJECT_NAME}_gdsp_static
  EXPORT ${PROJECT_NAME}_gd-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_gdsp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_gdsp/lib"
  )
install(
  EXPORT ${PROJECT_NAME}_gd-export
  NAMESPACE ${PROJECT_NAME}_gd::
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")
install(DIRECTORY ${LIB_MOD_GDDP_DIR}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")
install(DIRECTORY ${LIB_MOD_GDSP_DIR}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")

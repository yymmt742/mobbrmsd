set(LIB_MOD_2DDP_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_2ddp/mod_files)
if(NOT EXISTS "${LIB_MOD_2DDP_DIR}")
  make_directory("${LIB_MOD_2DDP_DIR}")
endif()
set(LIB_MOD_2DSP_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_2dsp/mod_files)
if(NOT EXISTS "${LIB_MOD_2DSP_DIR}")
  make_directory("${LIB_MOD_2DSP_DIR}")
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

add_library(${PROJECT_NAME}_2dsp OBJECT ${sources})
add_library(${PROJECT_NAME}_2ddp OBJECT ${sources})
set_target_properties(${PROJECT_NAME}_2dsp PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_2DSP_DIR})
set_target_properties(${PROJECT_NAME}_2ddp PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_2DDP_DIR})
ADDOPTS(${PROJECT_NAME}_2dsp 2 32 32 "" 1)
ADDOPTS(${PROJECT_NAME}_2ddp 2 64 32 "" 1)
add_library(${PROJECT_NAME}_2dsp_static STATIC $<TARGET_OBJECTS:${PROJECT_NAME}_2dsp>)
add_library(${PROJECT_NAME}_2dsp_shared SHARED $<TARGET_OBJECTS:${PROJECT_NAME}_2dsp>)
add_library(${PROJECT_NAME}_2ddp_static STATIC $<TARGET_OBJECTS:${PROJECT_NAME}_2ddp>)
add_library(${PROJECT_NAME}_2ddp_shared SHARED $<TARGET_OBJECTS:${PROJECT_NAME}_2ddp>)
set_target_properties(${PROJECT_NAME}_2ddp_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_2ddp)
set_target_properties(${PROJECT_NAME}_2ddp_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_2ddp)
set_target_properties(${PROJECT_NAME}_2dsp_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_2dsp)
set_target_properties(${PROJECT_NAME}_2dsp_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_2dsp)
install(
  TARGETS ${PROJECT_NAME}_2ddp_shared
  EXPORT ${PROJECT_NAME}-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2ddp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2ddp/lib"
  )
install(
  TARGETS ${PROJECT_NAME}_2ddp_static
  EXPORT ${PROJECT_NAME}-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2ddp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2ddp/lib"
  )
install(
  TARGETS ${PROJECT_NAME}_2dsp_shared
  EXPORT ${PROJECT_NAME}-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2dsp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2dsp/lib"
  )
install(
  TARGETS ${PROJECT_NAME}_2dsp_static
  EXPORT ${PROJECT_NAME}-export
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2dsp/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2dsp/lib"
  )
install(
  EXPORT ${PROJECT_NAME}-export
  NAMESPACE ${PROJECT_NAME}_2ddp::
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2ddp")
install(DIRECTORY ${LIB_MOD_2DDP_DIR}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2ddp")

install(
  EXPORT ${PROJECT_NAME}-export
  NAMESPACE ${PROJECT_NAME}_2dsp::
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2dsp")
install(DIRECTORY ${LIB_MOD_2DSP_DIR}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2dsp")


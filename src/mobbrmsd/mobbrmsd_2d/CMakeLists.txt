function(BUILD_2D prec)
  set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_2d${prec}/mod_files)
  if(NOT EXISTS "${LIB_MOD_DIR}")
    make_directory("${LIB_MOD_DIR}")
  endif()
  add_library(${PROJECT_NAME}_2d${prec}
    OBJECT
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
  set_target_properties(${PROJECT_NAME}_2d${prec} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR})
  if(prec STREQUAL "sp")
    ADDOPTS(${PROJECT_NAME}_2d${prec} 2 32 32 "" 1)
  elseif(prec STREQUAL "dp")
    ADDOPTS(${PROJECT_NAME}_2d${prec} 2 64 32 "" 1)
  endif()
  add_library(${PROJECT_NAME}_2d${prec}_static STATIC $<TARGET_OBJECTS:${PROJECT_NAME}_2d${prec}>)
  add_library(${PROJECT_NAME}_2d${prec}_shared SHARED $<TARGET_OBJECTS:${PROJECT_NAME}_2d${prec}>)
  set_target_properties(${PROJECT_NAME}_2d${prec}_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_2d${prec})
  set_target_properties(${PROJECT_NAME}_2d${prec}_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_2d${prec})
  install(
      TARGETS ${PROJECT_NAME}_2d${prec}_shared
    EXPORT ${PROJECT_NAME}-export
    ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2d${prec}/lib"
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2d${prec}/lib"
    )
  install(
    TARGETS ${PROJECT_NAME}_2d${prec}_static
    EXPORT ${PROJECT_NAME}-export
    ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2d${prec}/lib"
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2d${prec}/lib"
    )
  install(
    EXPORT ${PROJECT_NAME}-export
    NAMESPACE ${PROJECT_NAME}_2d${prec}::
    DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2d${prec}"
    )
  install(
    DIRECTORY ${LIB_MOD_DIR}
    DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_2d${prec}"
    )
endfunction(BUILD_2D)

foreach(prec IN ITEMS ${mobbrmsd_PRECISIONS})
  BUILD_2D(${prec})
endforeach()

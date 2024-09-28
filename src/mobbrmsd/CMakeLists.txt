function(BUILD_LIB dims prec)
  set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd_${dims}/${PROJECT_NAME}_${dims}${prec}/mod_files)
  if(NOT EXISTS "${LIB_MOD_DIR}")
    make_directory("${LIB_MOD_DIR}")
  endif()
  add_library(${PROJECT_NAME}_${dims}${prec}
    OBJECT
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_group_permutation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_mol_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_Hungarian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_c_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_f_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_tree.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_bb_block.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_bb_list.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_mobbrmsd_state.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_mobbrmsd.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_mobbrmsd_batch_run.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_mobbrmsd_mst.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mobbrmsd_${dims}/mod_dimspec_functions.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mobbrmsd_${dims}/mod_rotation.f90
  )
  set_target_properties(${PROJECT_NAME}_${dims}${prec}
                        PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR})
  if(prec STREQUAL "sp")
    if(dims STREQUAL "xd")
      ADDOPTS(${PROJECT_NAME}_${dims}${prec} 2 32 32 1 1)
    else()
      ADDOPTS(${PROJECT_NAME}_${dims}${prec} 2 32 32 "" 1)
    endif()
  elseif(prec STREQUAL "dp")
    if(dims STREQUAL "xd")
      ADDOPTS(${PROJECT_NAME}_${dims}${prec} 2 64 32 1 1)
    else()
      ADDOPTS(${PROJECT_NAME}_${dims}${prec} 2 64 32 "" 1)
    endif()
  endif()
  add_library(${PROJECT_NAME}_${dims}${prec}_static STATIC $<TARGET_OBJECTS:${PROJECT_NAME}_${dims}${prec}>)
  add_library(${PROJECT_NAME}_${dims}${prec}_shared SHARED $<TARGET_OBJECTS:${PROJECT_NAME}_${dims}${prec}>)
  set_target_properties(${PROJECT_NAME}_${dims}${prec}_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_${dims}${prec})
  set_target_properties(${PROJECT_NAME}_${dims}${prec}_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_${dims}${prec})
  install(
      TARGETS ${PROJECT_NAME}_${dims}${prec}_shared
    EXPORT ${PROJECT_NAME}-export
    ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_${dims}${prec}/lib"
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_${dims}${prec}/lib"
    )
  install(
    TARGETS ${PROJECT_NAME}_${dims}${prec}_static
    EXPORT ${PROJECT_NAME}-export
    ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_${dims}${prec}/lib"
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_${dims}${prec}/lib"
    )
  install(
    EXPORT ${PROJECT_NAME}-export
    NAMESPACE ${PROJECT_NAME}_${dims}${prec}::
    DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_${dims}${prec}"
    )
  install(
    DIRECTORY ${LIB_MOD_DIR}
    DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}_${dims}${prec}"
    )
endfunction(BUILD_LIB)


foreach(prec IN ITEMS ${mobbrmsd_PRECISIONS})
  foreach(dims IN ITEMS ${mobbrmsd_DIMENTIONS})
    BUILD_LIB( ${dims} ${prec} )
    #add_subdirectory(mobbrmsd_${dims})
  endforeach()
endforeach()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_subdirectory(test)
endif()


function( BUILD_LIB dims prec )
  set( LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${dims}${prec}/mod_files )
  set( target ${PROJECT_NAME}_${dims}${prec} )
  if( NOT EXISTS "${LIB_MOD_DIR}" )
    make_directory( "${LIB_MOD_DIR}" )
  endif()
  if( dims STREQUAL "xd" )
    set( blas 1 )
  else()
    set( blas 0 )
  endif()
  add_library(${target}
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
    ${CMAKE_CURRENT_SOURCE_DIR}/lapack/mod_mobbrmsd_lapack_routines.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/lapack/mod_mobbrmsd_lapack.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mobbrmsd_${dims}/mod_dimspec_functions.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mobbrmsd_${dims}/mod_rotation.f90
  )
  set_target_properties(${target}
                        PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR})

  ADDOPTS( ${target} ${dims} ${prec} 32 ${blas} 1 )

  install(
    EXPORT ${target}-export
    NAMESPACE ${target}::
    DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${target}"
    )
  install(
    DIRECTORY ${LIB_MOD_DIR}
    DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${target}"
    )
endfunction(BUILD_LIB)

foreach(prec IN ITEMS ${mobbrmsd_PRECISIONS})
  foreach(dims IN ITEMS ${mobbrmsd_DIMENTIONS})
    BUILD_LIB( ${dims} ${prec} )
  endforeach()
endforeach()

option( mobbrmsd_BUILD_SHARED_LIBS "build mylib as a shared library" ${BUILD_SHARED_LIBS})

if( mobbrmsd_BUILD_SHARED_LIBS )
  add_subdirectory(shared)
elseif( NOT mobbrmsd_BUILD_SHARED_LIBS OR SKBUILD )
  add_subdirectory(static)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_subdirectory(test)
endif()


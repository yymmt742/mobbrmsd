list( APPEND f90_src mod_params.f90 )
list( APPEND f90_src mod_group_permutation.f90 )
list( APPEND f90_src mod_mol_symmetry.f90 )
list( APPEND f90_src mod_det.f90 )
list( APPEND f90_src mod_estimate_rotation_matrix.f90 )
list( APPEND f90_src mod_mol_block.f90 )
list( APPEND f90_src mod_Hungarian.f90 )
list( APPEND f90_src mod_d_matrix.f90 )
list( APPEND f90_src mod_tree.f90 )
list( APPEND f90_src mod_branch_and_bound.f90 )
list( APPEND f90_src mod_symRMSD.f90 )

set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files)

if(NOT EXISTS "${LIB_MOD_DIR}")
  make_directory("${LIB_MOD_DIR}")
endif()

add_library(${PROJECT_NAME} OBJECT)
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${LIB_MOD_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)

target_sources(${PROJECT_NAME} PRIVATE ${f90_src})
target_include_directories(${PROJECT_NAME} PRIVATE include)
target_link_libraries(${PROJECT_NAME} PRIVATE ${BLAS_LIBRARIES})

add_library(${PROJECT_NAME}_shared SHARED $<TARGET_OBJECTS:symRMSD>)
set_target_properties(${PROJECT_NAME}_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_shared PRIVATE ${BLAS_LIBRARIES})
add_library(${PROJECT_NAME}_static STATIC $<TARGET_OBJECTS:symRMSD>)
set_target_properties(${PROJECT_NAME}_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME}_static
        EXPORT ${PROJECT_NAME}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}/lib")

install(TARGETS ${PROJECT_NAME}_shared
        EXPORT ${PROJECT_NAME}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

install(EXPORT ${PROJECT_NAME}-export
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")

install(DIRECTORY ${LIB_MOD_DIR}
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")

install(EXPORT ${PROJECT_NAME}-export
  FILE ${PROJECT_NAME}-config.cmake
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_subdirectory(test)
endif()

list( APPEND f90_src mod_kinds.f90 )
list( APPEND f90_src mod_params.f90 )
list( APPEND f90_src mod_group_permutation.f90 )
list( APPEND f90_src mod_mol_block.f90 )
list( APPEND f90_src mod_Hungarian.f90 )
list( APPEND f90_src mod_c_matrix.f90 )
list( APPEND f90_src mod_f_matrix.f90 )
list( APPEND f90_src mod_tree.f90 )
list( APPEND f90_src mod_bb_block.f90 )
list( APPEND f90_src mod_bb_list.f90 )

set(LIB_MOD_2D_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files_2d)
if(NOT EXISTS "${LIB_MOD_2D_DIR}")
  make_directory("${LIB_MOD_2D_DIR}")
endif()
set(LIB_MOD_3D_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files_3d)
if(NOT EXISTS "${LIB_MOD_3D_DIR}")
  make_directory("${LIB_MOD_3D_DIR}")
endif()
set(LIB_MOD_GD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files_g)
if(NOT EXISTS "${LIB_MOD_GD_DIR}")
  make_directory("${LIB_MOD_GD_DIR}")
endif()
set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
if(NOT EXISTS "${LIB_MOD_DIR}")
  make_directory("${LIB_MOD_DIR}")
endif()

set(TARGET_2D ${PROJECT_NAME}2d)
set(TARGET_3D ${PROJECT_NAME}3d)
set(TARGET_GD ${PROJECT_NAME}gd)
set(TARGET_INTERFACE ${PROJECT_NAME})

add_library(${TARGET_2D} OBJECT)
set_target_properties(${TARGET_2D} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_2D_DIR})
target_sources(${TARGET_2D} PRIVATE ${f90_src})
add_subdirectory(d2routines)
target_include_directories(${TARGET_2D} PRIVATE include)
target_include_directories(${TARGET_2D} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} )
target_include_directories(${TARGET_2D} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/d2routines )

add_library(${TARGET_3D} OBJECT)
set_target_properties(${TARGET_3D} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_3D_DIR})
target_sources(${TARGET_3D} PRIVATE ${f90_src})
add_subdirectory(d3routines)
target_include_directories(${TARGET_3D} PRIVATE include)

add_library(${TARGET_GD} OBJECT)
set_target_properties(${TARGET_GD} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_GD_DIR})
target_sources(${TARGET_GD} PRIVATE ${f90_src})
add_subdirectory(dgroutines)
target_include_directories(${TARGET_GD} PRIVATE include)
target_link_libraries(${TARGET_GD} PRIVATE ${BLAS_LIBRARIES})

add_library(${TARGET_INTERFACE} OBJECT)
set_target_properties(${TARGET_INTERFACE} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR})
target_include_directories(${TARGET_INTERFACE} PUBLIC 
  $<BUILD_INTERFACE:${LIB_MOD_DIR}>
  $<BUILD_INTERFACE:${LIB_MOD_2D_DIR}>
  $<BUILD_INTERFACE:${LIB_MOD_3D_DIR}>
  $<BUILD_INTERFACE:${LIB_MOD_GD_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)
target_sources(${TARGET_INTERFACE} PRIVATE mod_mobbrmsd.f90 )
target_link_libraries(${TARGET_INTERFACE} PRIVATE ${TARGET_2D})
target_link_libraries(${TARGET_INTERFACE} PRIVATE ${TARGET_3D})
target_link_libraries(${TARGET_INTERFACE} PRIVATE ${TARGET_GD})
target_link_libraries(${TARGET_INTERFACE} PRIVATE ${BLAS_LIBRARIES})

add_library(${TARGET_INTERFACE}_shared SHARED $<TARGET_OBJECTS:mobbrmsd>)
set_target_properties(${TARGET_INTERFACE}_shared PROPERTIES OUTPUT_NAME ${TARGET_INTERFACE})
target_link_libraries(${TARGET_INTERFACE}_shared PRIVATE ${BLAS_LIBRARIES})
add_library(${TARGET_INTERFACE}_static STATIC $<TARGET_OBJECTS:mobbrmsd>)
set_target_properties(${TARGET_INTERFACE}_static PROPERTIES OUTPUT_NAME ${TARGET_INTERFACE})

install(TARGETS ${TARGET_INTERFACE}_static
        EXPORT ${TARGET_INTERFACE}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}/lib")
install(TARGETS ${TARGET_INTERFACE}_shared
        EXPORT ${TARGET_INTERFACE}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install(EXPORT ${TARGET_INTERFACE}-export
        NAMESPACE ${TARGET_INTERFACE}::
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")
install(EXPORT ${TARGET_INTERFACE}-export
        FILE ${TARGET_INTERFACE}-config.cmake
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")

install(DIRECTORY ${LIB_MOD_DIR}
  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_subdirectory(test)
endif()

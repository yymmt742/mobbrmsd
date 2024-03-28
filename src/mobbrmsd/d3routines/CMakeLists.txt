list( APPEND f90_3d_src gemm.f90 )
list( APPEND f90_3d_src blas_lapack_interface.f90 )
list( APPEND f90_3d_src mod_rotation.f90 )
list( APPEND f90_3d_src mod_bb_interface_3d.f90 )

target_sources(${TARGET_3D} PRIVATE ${f90_3d_src})

set(LIB_MOD_3D_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
if(NOT EXISTS "${LIB_MOD_3D_DIR}")
  make_directory("${LIB_MOD_3D_DIR}")
endif()

set_target_properties(${TARGET_3D} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_3D_DIR})
target_include_directories(${TARGET_3D} PRIVATE include)

set_target_properties(${TARGET_3D_INTERFACE} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_3D_DIR})
target_include_directories(${TARGET_3D_INTERFACE} PUBLIC $<BUILD_INTERFACE:${LIB_MOD_3D_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)
target_link_libraries(${TARGET_3D_INTERFACE} PRIVATE ${TARGET_3D})

add_library(${TARGET_3D_INTERFACE}_shared SHARED $<TARGET_OBJECTS:mobbrmsd_2d>)
set_target_properties(${TARGET_3D_INTERFACE}_shared PROPERTIES OUTPUT_NAME ${TARGET_3D_INTERFACE})

add_library(${TARGET_3D_INTERFACE}_static STATIC $<TARGET_OBJECTS:mobbrmsd_2d>)
set_target_properties(${TARGET_3D_INTERFACE}_static PROPERTIES OUTPUT_NAME ${TARGET_3D_INTERFACE})

install(TARGETS ${TARGET_3D_INTERFACE}_static
        EXPORT ${TARGET_3D_INTERFACE}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}/lib")
install(TARGETS ${TARGET_3D_INTERFACE}_shared
        EXPORT ${TARGET_3D_INTERFACE}-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install(EXPORT ${TARGET_3D_INTERFACE}-export
        NAMESPACE ${TARGET_3D_INTERFACE}::
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")
install(EXPORT ${TARGET_3D_INTERFACE}-export
        FILE ${TARGET_3D_INTERFACE}-config.cmake
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")

install(DIRECTORY ${LIB_MOD_3D_DIR} DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${TARGET_3D}")

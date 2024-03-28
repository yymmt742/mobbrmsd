#set(TARGET_3D ${PROJECT_NAME}_3d)
#add_library(${TARGET_3D} OBJECT)
#target_sources(${TARGET_3D} PRIVATE ${f90_src})
#target_include_directories(${TARGET_3D} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

#set(TARGET_G ${PROJECT_NAME})
#add_library(${TARGET_G}  OBJECT)
#target_sources(${TARGET_G}  PRIVATE ${f90_src})
#target_include_directories(${TARGET_GD} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(mobbrmsd_2d)
#add_subdirectory(d3routines)
#add_subdirectory(dgroutines)

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_subdirectory(test)
endif()

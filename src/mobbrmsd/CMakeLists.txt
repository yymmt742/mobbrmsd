add_subdirectory(mobbrmsd_2d)
add_subdirectory(mobbrmsd_3d)
add_subdirectory(mobbrmsd_gd)
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_subdirectory(test)
endif()


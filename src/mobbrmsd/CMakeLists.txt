add_subdirectory(mobbrmsd_2d)
add_subdirectory(mobbrmsd_3d)
add_subdirectory(mobbrmsd_xd)
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_subdirectory(test)
endif()


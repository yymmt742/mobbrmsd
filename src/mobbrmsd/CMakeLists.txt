add_subdirectory(mobbrmsd_2d)
add_subdirectory(mobbrmsd_3d)
add_subdirectory(mobbrmsd_gd)
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_subdirectory(test)
endif()
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}")


foreach(dims IN ITEMS ${mobbrmsd_DIMENTIONS})
 add_subdirectory(mobbrmsd_${dims})
endforeach()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_subdirectory(test)
endif()


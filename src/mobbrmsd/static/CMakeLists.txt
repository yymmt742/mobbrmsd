foreach(prec IN ITEMS ${mobbrmsd_PRECISIONS})
  foreach(dims IN ITEMS ${mobbrmsd_DIMENTIONS})
    set( target ${PROJECT_NAME}_${dims}${prec} )
    set( target_static ${target}_static )
    add_library(${target_static} STATIC $<TARGET_OBJECTS:${target}>)
    set_target_properties(${target_static} PROPERTIES OUTPUT_NAME ${target})
    install(
      TARGETS ${target_static}
      EXPORT ${target}-export
      ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${target}/lib"
      LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/${target}/lib"
      )
  endforeach()
endforeach()


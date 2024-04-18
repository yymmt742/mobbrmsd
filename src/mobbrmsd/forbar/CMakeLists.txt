set(FORBAR_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
if(NOT EXISTS "${FORBAR_MOD_DIR}")
  make_directory("${FORBAR_MOD_DIR}")
endif()

add_library(forbar STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/isatty.c
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_optarg.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_iolib.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_word_iterator.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_repeat_iterator.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_cycle_iterator.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_progbar_iterator.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/mod_forbar.f90
  )

target_include_directories(forbar PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(forbar PROPERTIES Fortran_MODULE_DIRECTORY ${FORBAR_MOD_DIR})

install(TARGETS forbar
        EXPORT forbar-export
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/forbar/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/forbar/lib")
install(EXPORT forbar-export
        NAMESPACE forbar::
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/forbar")

#install(DIRECTORY ${FORBAR_MOD_DIR} DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/forbar")

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  find_package( unittest REQUIRED )

  add_executable( test_progress_bar )
  target_sources( test_progress_bar PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/isatty.c
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_optarg.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_iolib.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_word_iterator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_repeat_iterator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_cycle_iterator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_progbar_iterator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_forbar.f90
    test_progress_bar.f90 )
  target_include_directories( test_progress_bar PRIVATE ${unittest_DIR}/mod_files )
  target_link_libraries( test_progress_bar PRIVATE unittest )
  add_test( NAME test_progress_bar COMMAND test_progress_bar )
endif()

find_package( unittest REQUIRED )

add_executable( test_cov )
target_include_directories( test_cov PRIVATE ${unittest_DIR}/mod_files )

target_sources( test_cov PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90 )
target_sources( test_cov PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_cov.f90 )
target_sources( test_cov PRIVATE test_cov.f90 )
target_link_libraries( test_cov PRIVATE unittest )

add_test( NAME test_cov
          COMMAND test_cov )


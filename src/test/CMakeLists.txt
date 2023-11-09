find_package( unittest REQUIRED )

add_executable( test_cov_time_1 )
target_include_directories( test_cov_time_1 PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_cov_time_1 PRIVATE cov_1.f90 test_cov_time.f90 )
target_link_libraries( test_cov_time_1 PRIVATE unittest )
target_link_libraries( test_cov_time_1 PRIVATE blas )
add_test( NAME test_cov_time_1 COMMAND test_cov_time_1 )

add_executable( test_cov_time_2 )
target_include_directories( test_cov_time_2 PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_cov_time_2 PRIVATE cov_2.f90 test_cov_time.f90 )
target_link_libraries( test_cov_time_2 PRIVATE unittest )
add_test( NAME test_cov_time_2 COMMAND test_cov_time_2 )

add_executable( test_cov_time_3 )
target_include_directories( test_cov_time_3 PRIVATE ${unittest_DIR}/mod_files )
target_include_directories( test_cov_time_3 PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_cov_time_3 PRIVATE ${purelapack_DIR}/include )
target_sources( test_cov_time_3 PRIVATE cov_3.f90 test_cov_time.f90 )
target_link_libraries( test_cov_time_3 PRIVATE purelapack )
target_link_libraries( test_cov_time_3 PRIVATE unittest )
add_test( NAME test_cov_time_3 COMMAND test_cov_time_3 )

add_executable( test_cov_time_4 )
target_include_directories( test_cov_time_4 PRIVATE ${unittest_DIR}/mod_files )
target_include_directories( test_cov_time_4 PRIVATE ${purelapack_DIR}/include )
target_include_directories( test_cov_time_4 PRIVATE ${purelapack_DIR}/mod_files )
target_sources( test_cov_time_4 PRIVATE cov_3.f90 test_cov_time.f90 )
target_link_libraries( test_cov_time_4 PRIVATE blas )
target_link_libraries( test_cov_time_4 PRIVATE unittest )
add_test( NAME test_cov_time_4 COMMAND test_cov_time_4 )

add_executable( test_cov )
target_include_directories( test_cov PRIVATE ${unittest_DIR}/mod_files )
target_include_directories( test_cov PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_cov PRIVATE ${purelapack_DIR}/include )
target_sources( test_cov PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../mod_cov.f90
                                 test_cov.f90 )
target_link_libraries( test_cov PRIVATE unittest )
target_link_libraries( test_cov PRIVATE blas )
add_test( NAME test_cov COMMAND test_cov )

add_executable( test_det )
target_include_directories( test_det PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_det PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                 test_det.f90 )
target_link_libraries( test_det PRIVATE unittest )
target_include_directories( test_det PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_det PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_det PRIVATE purelapack )
add_test( NAME test_det COMMAND test_det )

add_executable( test_svd )
target_include_directories( test_svd PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_svd PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                 test_svd.f90 )
target_link_libraries( test_svd PRIVATE unittest )
target_include_directories( test_svd PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_svd PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_svd PRIVATE purelapack )
add_test( NAME test_svd COMMAND test_svd )

add_executable( test_rot )
target_include_directories( test_rot PRIVATE ${unittest_DIR}/mod_files )
target_include_directories( test_rot PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_rot PRIVATE ${purelapack_DIR}/include )
target_sources( test_rot PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rot.f90
                                 test_rot.f90 )
target_link_libraries( test_rot PRIVATE unittest )
target_link_libraries( test_rot PRIVATE blas )
add_test( NAME test_rot COMMAND test_rot )

add_executable( test_Kabsch )
target_include_directories( test_Kabsch PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_Kabsch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_cov.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rot.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                    test_Kabsch.f90 )
target_link_libraries( test_Kabsch PRIVATE unittest )
target_include_directories( test_Kabsch PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_Kabsch PRIVATE ${purelapack_DIR}/include )
#target_link_libraries( test_Kabsch PRIVATE purelapack )
target_link_libraries( test_Kabsch PRIVATE blas )
add_test( NAME test_Kabsch COMMAND test_Kabsch )

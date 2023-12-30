find_package( unittest REQUIRED )
find_package( purelapack REQUIRED )

add_executable( test_optarg )
target_include_directories( test_optarg PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_optarg PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                    test_optarg.f90 )
target_link_libraries( test_optarg PRIVATE unittest )
add_test( NAME test_optarg COMMAND test_optarg )

add_executable( test_group_permutation )
target_include_directories( test_group_permutation PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_group_permutation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                    test_group_permutation.f90 )
target_link_libraries( test_group_permutation PRIVATE unittest )
add_test( NAME test_group_permutation COMMAND test_group_permutation )

add_executable( test_Hungarian )
target_include_directories( test_Hungarian PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_Hungarian PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                       ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                       test_Hungarian.f90 )
target_link_libraries( test_Hungarian PRIVATE unittest )
add_test( NAME test_Hungarian COMMAND test_Hungarian )

add_executable( test_det )
target_include_directories( test_det PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_det PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                 test_det.f90 )
target_include_directories( test_det PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_det PRIVATE unittest )
target_link_libraries( test_det PRIVATE purelapack )
add_test( NAME test_det COMMAND test_det )

add_executable( test_Kabsch )
target_include_directories( test_Kabsch PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_Kabsch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                    test_Kabsch.f90 )
target_include_directories( test_Kabsch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_Kabsch PRIVATE unittest )
target_link_libraries( test_Kabsch PRIVATE purelapack )
add_test( NAME test_Kabsch COMMAND test_Kabsch )

add_executable( test_d_matrix )
target_include_directories( test_d_matrix PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_d_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_molecular_rotation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_d_matrix.f90
                                      test_d_matrix.f90 )
target_include_directories( test_d_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_d_matrix PRIVATE unittest )
target_link_libraries( test_d_matrix PRIVATE purelapack )
add_test( NAME test_d_matrix COMMAND test_d_matrix )

add_executable( test_tree )
target_include_directories( test_tree PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_tree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
                                  test_tree.f90 )
target_link_libraries( test_tree PRIVATE unittest )
add_test( NAME test_tree COMMAND test_tree )

add_executable( test_branch_and_prune )
target_include_directories( test_branch_and_prune PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_branch_and_prune PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_molecular_rotation.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_d_matrix.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_branch_and_prune.f90
                                  test_branch_and_prune.f90 )
target_include_directories( test_branch_and_prune PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_branch_and_prune PRIVATE unittest )
target_link_libraries( test_branch_and_prune PRIVATE purelapack )
add_test( NAME test_branch_and_prune COMMAND test_branch_and_prune )

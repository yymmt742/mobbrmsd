find_package( unittest REQUIRED )
find_package( purelapack REQUIRED )

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

add_executable( test_mol_block )
target_include_directories( test_mol_block PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_mol_block PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                       ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                       test_mol_block.f90 )
target_link_libraries( test_mol_block PRIVATE unittest )
add_test( NAME test_mol_block COMMAND test_mol_block )

add_executable( test_det )
target_include_directories( test_det PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_det PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                 test_det.f90 )
target_include_directories( test_det PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_det PRIVATE unittest )
target_link_libraries( test_det PRIVATE purelapack )
add_test( NAME test_det COMMAND test_det )

add_executable( test_tree )
target_include_directories( test_tree PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_tree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
                                  test_tree.f90 )
target_include_directories( test_tree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_tree PRIVATE unittest )
add_test( NAME test_tree COMMAND test_tree )

add_executable( test_rotation_matrix )
target_include_directories( test_rotation_matrix PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_rotation_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rotation_matrix.f90
                                    mod_testutil.f90
                                    test_rotation_matrix.f90 )
target_include_directories( test_rotation_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_rotation_matrix PRIVATE unittest )
target_link_libraries( test_rotation_matrix PRIVATE purelapack )
add_test( NAME test_rotation_matrix COMMAND test_rotation_matrix )

add_executable( test_lowerbound )
target_include_directories( test_lowerbound PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_lowerbound PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rotation_matrix.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_lowerbound.f90
                                        mod_testutil.f90
                                        test_lowerbound.f90 )
target_include_directories( test_lowerbound PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_lowerbound PRIVATE unittest )
target_link_libraries( test_lowerbound PRIVATE purelapack )
add_test( NAME test_lowerbound COMMAND test_lowerbound )

add_executable( test_c_matrix )
target_include_directories( test_c_matrix PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_c_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rotation_matrix.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_symmetry.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
                                      mod_testutil.f90
                                      test_c_matrix.f90 )
target_include_directories( test_c_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_c_matrix PRIVATE unittest )
target_link_libraries( test_c_matrix PRIVATE purelapack )
add_test( NAME test_c_matrix COMMAND test_c_matrix )

add_executable( test_s_matrix )
target_include_directories( test_s_matrix PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_s_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rotation_matrix.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_symmetry.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_s_matrix.f90
                                      mod_testutil.f90
                                      test_s_matrix.f90 )
target_include_directories( test_s_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_s_matrix PRIVATE unittest )
target_link_libraries( test_s_matrix PRIVATE purelapack )
add_test( NAME test_s_matrix COMMAND test_s_matrix )

add_executable( test_bb_manager )
target_include_directories( test_bb_manager PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_bb_manager PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rotation_matrix.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_symmetry.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_s_matrix.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_lowerbound.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_manager.f90
                                        mod_testutil.f90
                                        test_bb_manager.f90 )
target_include_directories( test_bb_manager PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_bb_manager PRIVATE unittest )
target_link_libraries( test_bb_manager PRIVATE purelapack )
add_test( NAME test_bb_manager COMMAND test_bb_manager )

add_executable( test_branch_and_bound )
target_include_directories( test_branch_and_bound PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_branch_and_bound PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rotation_matrix.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_symmetry.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_s_matrix.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_branch_and_bound.f90
                                  test_branch_and_bound.f90 )
target_include_directories( test_branch_and_bound PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_branch_and_bound PRIVATE unittest )
target_link_libraries( test_branch_and_bound PRIVATE purelapack )
add_test( NAME test_branch_and_bound COMMAND test_branch_and_bound )

add_executable( test_symRMSD )
target_include_directories( test_symRMSD PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_symRMSD PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rotation_matrix.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_symmetry.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_s_matrix.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_branch_and_bound.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/../mod_symRMSD.f90
                                     test_symRMSD.f90 )
target_include_directories( test_symRMSD PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include )
target_link_libraries( test_symRMSD PRIVATE unittest )
target_link_libraries( test_symRMSD PRIVATE purelapack )
add_test( NAME test_symRMSD COMMAND test_symRMSD )

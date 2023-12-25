find_package( unittest REQUIRED )

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

add_executable( test_rmsd )
target_include_directories( test_rmsd PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_rmsd PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd.f90
                                 test_rmsd.f90 )
target_link_libraries( test_rmsd PRIVATE unittest )
add_test( NAME test_rmsd COMMAND test_rmsd )

add_executable( test_rmsd_brute )
target_include_directories( test_rmsd_brute PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_rmsd_brute PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_permutation.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd_brute.f90
                                        test_rmsd_brute.f90 )
target_link_libraries( test_rmsd_brute PRIVATE unittest )
add_test( NAME test_rmsd_brute COMMAND test_rmsd_brute )

#add_executable( test_cov_time_1 )
#target_include_directories( test_cov_time_1 PRIVATE ${unittest_DIR}/mod_files )
#target_sources( test_cov_time_1 PRIVATE cov_1.f90 test_cov_time.f90 )
#target_link_libraries( test_cov_time_1 PRIVATE unittest )
#add_test( NAME test_cov_time_1 COMMAND test_cov_time_1 )

#add_executable( test_cov_time_2 )
#target_include_directories( test_cov_time_2 PRIVATE ${unittest_DIR}/mod_files )
#target_sources( test_cov_time_2 PRIVATE cov_2.f90 test_cov_time.f90 )
#target_link_libraries( test_cov_time_2 PRIVATE unittest )
#add_test( NAME test_cov_time_2 COMMAND test_cov_time_2 )

#add_executable( test_cov_time_3 )
#target_include_directories( test_cov_time_3 PRIVATE ${unittest_DIR}/mod_files )
#target_include_directories( test_cov_time_3 PRIVATE ${purelapack_DIR}/mod_files )
#target_include_directories( test_cov_time_3 PRIVATE ${purelapack_DIR}/include )
#target_sources( test_cov_time_3 PRIVATE cov_3.f90 test_cov_time.f90 )
#target_link_libraries( test_cov_time_3 PRIVATE purelapack )
#target_link_libraries( test_cov_time_3 PRIVATE unittest )
#add_test( NAME test_cov_time_3 COMMAND test_cov_time_3 )

#add_executable( test_cov_time_4 )
#target_include_directories( test_cov_time_4 PRIVATE ${unittest_DIR}/mod_files )
#target_include_directories( test_cov_time_4 PRIVATE ${purelapack_DIR}/include )
#target_include_directories( test_cov_time_4 PRIVATE ${purelapack_DIR}/mod_files )
#target_sources( test_cov_time_4 PRIVATE cov_3.f90 test_cov_time.f90 )
#target_link_libraries( test_cov_time_4 PRIVATE unittest )
#add_test( NAME test_cov_time_4 COMMAND test_cov_time_4 )

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
target_link_libraries( test_det PRIVATE unittest )
target_include_directories( test_det PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_det PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_det PRIVATE purelapack )
add_test( NAME test_det COMMAND test_det )

add_executable( test_svd )
target_include_directories( test_svd PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_svd PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                 test_svd.f90 )
target_link_libraries( test_svd PRIVATE unittest )
target_include_directories( test_svd PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_svd PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_svd PRIVATE purelapack )
add_test( NAME test_svd COMMAND test_svd )

add_executable( test_Kabsch )
target_include_directories( test_Kabsch PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_Kabsch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                    test_Kabsch.f90 )
target_link_libraries( test_Kabsch PRIVATE unittest )
target_include_directories( test_Kabsch PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_Kabsch PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_Kabsch PRIVATE purelapack )
add_test( NAME test_Kabsch COMMAND test_Kabsch )

add_executable( test_procrustes )
target_include_directories( test_procrustes PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_procrustes PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_procrustes.f90
                                    test_procrustes.f90 )
target_link_libraries( test_procrustes PRIVATE unittest )
target_include_directories( test_procrustes PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_procrustes PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_procrustes PRIVATE purelapack )
add_test( NAME test_procrustes COMMAND test_procrustes )

add_executable( test_lower_bound )
target_include_directories( test_lower_bound PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_lower_bound PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_permutation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd_brute.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_procrustes.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_lower_bound.f90
                                      test_lower_bound.f90 )
target_link_libraries( test_lower_bound PRIVATE unittest )
target_include_directories( test_lower_bound PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_lower_bound PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_lower_bound PRIVATE purelapack )
add_test( NAME test_lower_bound COMMAND test_lower_bound )

add_executable( test_d_matrix )
target_include_directories( test_d_matrix PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_d_matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_molecular_rotation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_d_matrix.f90
                                      test_d_matrix.f90 )
target_link_libraries( test_d_matrix PRIVATE unittest )
target_include_directories( test_d_matrix PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_d_matrix PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_d_matrix PRIVATE purelapack )
add_test( NAME test_d_matrix COMMAND test_d_matrix )

add_executable( test_block_lower_bound )
target_include_directories( test_block_lower_bound PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_block_lower_bound PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_permutation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd_brute.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_procrustes.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_molecular_rotation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_block_lower_bound.f90
                                      test_block_lower_bound.f90 )
target_link_libraries( test_block_lower_bound PRIVATE unittest )
target_include_directories( test_block_lower_bound PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_block_lower_bound PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_block_lower_bound PRIVATE purelapack )
add_test( NAME test_block_lower_bound COMMAND test_block_lower_bound )

add_executable( test_partial_rmsd )
target_include_directories( test_partial_rmsd PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_partial_rmsd PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_permutation.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd_brute.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_procrustes.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_partial_rmsd.f90
                                      test_partial_rmsd.f90 )
target_link_libraries( test_partial_rmsd PRIVATE unittest )
target_include_directories( test_partial_rmsd PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_partial_rmsd PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_partial_rmsd PRIVATE purelapack )
add_test( NAME test_partial_rmsd COMMAND test_partial_rmsd )

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
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_permutation.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_rmsd_brute.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_det.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_svd.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Kabsch.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_molecular_rotation.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_d_matrix.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_branch_and_prune.f90
                                  test_branch_and_prune.f90 )
target_link_libraries( test_branch_and_prune PRIVATE unittest )
target_include_directories( test_branch_and_prune PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( test_branch_and_prune PRIVATE ${purelapack_DIR}/include )
target_link_libraries( test_branch_and_prune PRIVATE purelapack )
add_test( NAME test_branch_and_prune COMMAND test_branch_and_prune )

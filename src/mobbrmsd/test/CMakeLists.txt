find_package( unittest REQUIRED )

set(TARGET_2D ${PROJECT_NAME}_2d)
set(TARGET_3D ${PROJECT_NAME}_3d)

set(LIB_MOD_2D_DIR ${CMAKE_CURRENT_BINARY_DIR}/../${TARGET_2D}/mod_files)
set(LIB_MOD_3D_DIR ${CMAKE_CURRENT_BINARY_DIR}/../${TARGET_3D}/mod_files)

add_executable( test_group_permutation )
target_sources( test_group_permutation PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
  test_group_permutation.f90 )
target_include_directories( test_group_permutation PRIVATE ${unittest_DIR}/mod_files )
target_link_libraries( test_group_permutation PRIVATE unittest )
add_test( NAME test_group_permutation COMMAND test_group_permutation )

add_executable( test_Hungarian )
target_sources( test_Hungarian PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
  test_Hungarian.f90 )
target_include_directories( test_Hungarian PRIVATE ${unittest_DIR}/mod_files )
target_link_libraries( test_Hungarian PRIVATE unittest )
add_test( NAME test_Hungarian COMMAND test_Hungarian )

add_executable( test_tree )
target_sources( test_tree PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
  test_tree.f90 )
target_include_directories( test_tree PRIVATE ${unittest_DIR}/mod_files )
target_link_libraries( test_tree PRIVATE unittest )
add_test( NAME test_tree COMMAND test_tree )


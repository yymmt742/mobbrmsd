ADDTEST(test_group_permutation "" "" "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/test_group_permutation.f90
  )
ADDTEST(test_Hungarian "" "" "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/test_Hungarian.f90
  )
ADDTEST(test_tree "" "" "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/test_tree.f90
  )

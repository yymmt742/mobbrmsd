ADDTEST(test_invsqrt_sp X 32 32 "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../test_invsqrt.f90
  )
ADDTEST(test_invcbrt_sp X 32 32 "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../test_invcbrt.f90
  )
ADDTEST(test_quartic_sp X 32 32 "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../test_quartic.f90
  )

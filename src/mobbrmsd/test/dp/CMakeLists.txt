ADDTEST(test_invsqrt_dp 64 32 "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../test_invsqrt.f90
  )
ADDTEST(test_invcbrt_dp 64 32 "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../test_invcbrt.f90
  )
ADDTEST(test_quartic_dp 64 32 "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../test_quartic.f90
  )

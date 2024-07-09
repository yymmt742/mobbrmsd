DIMSPECTESTS(3ddp ${CMAKE_CURRENT_SOURCE_DIR}/../../mobbrmsd_3d 64 32)
ADDTEST(test_quartic_dp 64 32 "" "" ""
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_kinds.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../../mod_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../test_quartic.f90
  )

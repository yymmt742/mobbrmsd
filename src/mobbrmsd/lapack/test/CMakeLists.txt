foreach(prec IN ITEMS ${mobbrmsd_PRECISIONS})
  set(moddir ${CMAKE_CURRENT_BINARY_DIR}/${prec})
  ADDTEST(test_mobbrmsd_lapack_${prec} xd ${prec} 32 "" "" ${moddir} ""
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_lapack_routines.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_lapack.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/test_mobbrmsd_lapack.f90
    )
  target_include_directories(test_mobbrmsd_lapack_${prec} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lapack_routine)
endforeach()


foreach(prec IN ITEMS ${mobbrmsd_PRECISIONS})
  set(moddir ${CMAKE_CURRENT_BINARY_DIR}/${prec})
  set(sources "")
  list(APPEND sources
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
    )
  ADDTEST(test_group_permutation_${prec} xd ${prec} 32 "" "" ${moddir} ""
    ${sources}
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/test_group_permutation.f90
    )
  ADDTEST(test_Hungarian_${prec} xd ${prec} 32 "" "" ${moddir} ""
    ${sources}
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/test_Hungarian.f90
    )
  ADDTEST(test_tree_${prec} xd ${prec} 32 "" "" ${moddir} ""
    ${sources}
    ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/test_tree.f90
    )
endforeach()

foreach(dims IN ITEMS ${mobbrmsd_DIMENTIONS})
  foreach(prec IN ITEMS ${mobbrmsd_PRECISIONS})
    set(sources "")
    set(moddir ${CMAKE_CURRENT_BINARY_DIR}/${dims}${prec})
    list(APPEND sources
      ${CMAKE_CURRENT_SOURCE_DIR}/../lapack/mod_mobbrmsd_lapack_routines_${prec}.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../lapack/mod_mobbrmsd_lapack.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mobbrmsd_${dims}/mod_dimspec_functions.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_kinds.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
      )
    ADDTEST(test_covdot_${dims}${prec} ${dims} ${prec} 32 1 "" ${moddir} ""
      ${sources}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_covdot.f90
      )
    ADDTEST(test_gemm_${dims}${prec} ${dims} ${prec} 32 1 "" ${moddir} ""
      ${sources}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_gemm.f90
      )
    list(APPEND sources
      ${CMAKE_CURRENT_SOURCE_DIR}/mod_testutil.f90
      )
    ADDTEST(test_testutil_${dims}${prec} ${dims} ${prec} 32 1 "" ${moddir} ""
      ${sources}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_testutil.f90
      )
    list(APPEND sources
      ${CMAKE_CURRENT_SOURCE_DIR}/../mobbrmsd_${dims}/mod_rotation.f90
      )
    ADDTEST(test_rotation_${dims}${prec} ${dims} ${prec} 32 1 "" ${moddir} ""
      ${sources}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_rotation.f90
      )
    list(APPEND sources
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_group_permutation.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mol_block.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_c_matrix.f90
      )
    ADDTEST(test_c_matrix_${dims}${prec} ${dims} ${prec} 32 1 "" ${moddir} ""
      ${sources}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_c_matrix.f90
      )
    list(APPEND sources
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_Hungarian.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_f_matrix.f90
      )
    ADDTEST(test_f_matrix_${dims}${prec} ${dims} ${prec} 32 1 "" ${moddir} ""
      ${sources}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_f_matrix.f90
      )
    list(APPEND sources
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_tree.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_block.f90
      )
    ADDTEST(test_bb_block_${dims}${prec} ${dims} ${prec} 32 1 "" ${moddir} ""
      ${sources}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_bb_block.f90
      )
    list(APPEND sources
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bb_list.f90
      )
    ADDTEST(test_bb_list_${dims}${prec} ${dims} ${prec} 32 1 "" ${moddir} 4
      ${sources}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_bb_list.f90
      )
    list(APPEND sources
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_header.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_state.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_mst.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_mobbrmsd_batch_run.f90
      )
    ADDTEST(test_mobbrmsd_${dims}${prec} ${dims} ${prec} 32 1 1 ${moddir} 4
      ${sources}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_mobbrmsd.f90
      )
  endforeach()
endforeach()


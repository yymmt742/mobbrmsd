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

set(curr_dir ${CMAKE_CURRENT_SOURCE_DIR})

function(DIMSPECTESTS post dimspec_dir real int)
  set(sources "")
  list(APPEND sources
    ${dimspec_dir}/mod_dimspec_functions.f90
    ${curr_dir}/../mod_kinds.f90
    ${curr_dir}/../mod_params.f90
    )
  ADDTEST(test_covdot_${post} "${real}" "${int}" 1 "" ""
    ${sources}
    ${curr_dir}/test_covdot.f90
    )
  list(APPEND sources
    ${dimspec_dir}/mod_rotation.f90
    ${curr_dir}/mod_testutil.f90
    )
  ADDTEST(test_rotation_${post} "${real}" "${int}" 1 "" ""
    ${sources}
    ${curr_dir}/test_rotation.f90
    )
  list(APPEND sources
    ${curr_dir}/../mod_group_permutation.f90
    ${curr_dir}/../mod_mol_block.f90
    ${curr_dir}/../mod_c_matrix.f90
    )
  ADDTEST(test_c_matrix_${post} "${real}" "${int}" 1 "" ""
    ${sources}
    ${curr_dir}/test_c_matrix.f90
    )
  list(APPEND sources
    ${curr_dir}/../mod_Hungarian.f90
    ${curr_dir}/../mod_f_matrix.f90
    )
  ADDTEST(test_f_matrix_${post} "${real}" "${int}" 1 "" ""
    ${sources}
    ${curr_dir}/test_f_matrix.f90
    )
  list(APPEND sources
    ${curr_dir}/../mod_tree.f90
    ${curr_dir}/../mod_bb_block.f90
    )
  ADDTEST(test_bb_block_${post} "${real}" "${int}" 1 "" ""
    ${sources}
    ${curr_dir}/test_bb_block.f90
    )
  list(APPEND sources
    ${curr_dir}/../mod_bb_list.f90
    )
  ADDTEST(test_bb_list_${post} "${real}" "${int}" 1 "" 4
    ${sources}
    ${curr_dir}/test_bb_list.f90
    )
  list(APPEND sources
    ${curr_dir}/../mod_mobbrmsd_header.f90
    ${curr_dir}/../mod_mobbrmsd_state.f90
    ${curr_dir}/../mod_mobbrmsd.f90
    ${curr_dir}/../mod_mobbrmsd_mst.f90
    ${curr_dir}/../mod_mobbrmsd_batch_run.f90
    )
  ADDTEST(test_mobbrmsd_${post} "${real}" "${int}" 1 1 4
    ${sources}
    ${curr_dir}/test_mobbrmsd.f90
    )
endfunction(DIMSPECTESTS)

add_subdirectory(sp)
add_subdirectory(dp)
add_subdirectory(2dsp)
add_subdirectory(2ddp)
add_subdirectory(3dsp)
add_subdirectory(3ddp)
add_subdirectory(gdsp)
add_subdirectory(gddp)

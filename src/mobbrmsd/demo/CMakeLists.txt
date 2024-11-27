set(python_files)
PREPROC_PYTHON("${python_files}" python_files _demo)
PREPROC_PYTHON("${python_files}" python_files _coord_generator)
PREPROC_PYTHON("${python_files}" python_files cogen)
PREPROC_PYTHON("${python_files}" python_files bb)
PREPROC_PYTHON("${python_files}" python_files bb_2d)
PREPROC_PYTHON("${python_files}" python_files bb_multi)
PREPROC_PYTHON("${python_files}" python_files mst)
PREPROC_PYTHON("${python_files}" python_files batch)
PREPROC_PYTHON("${python_files}" python_files batch_tri)
install(FILES
  ${python_files}
  DESTINATION ${SKBUILD_PLATLIB_DIR}/mobbrmsd/demo
)

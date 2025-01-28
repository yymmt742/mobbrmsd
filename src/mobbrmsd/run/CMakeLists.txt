set(python_files)
PREPROC_PYTHON("${python_files}" python_files parse)
PREPROC_PYTHON("${python_files}" python_files verbose_run)
PREPROC_PYTHON("${python_files}" python_files batch_run)
install(FILES
  ${python_files}
  DESTINATION ${SKBUILD_PLATLIB_DIR}/mobbrmsd/run
)

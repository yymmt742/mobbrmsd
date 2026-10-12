if(SKBUILD)
  find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

  # Grab the variables from a local Python installation F2PY headers
  execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # Python >= 3.12 uses the meson f2py backend which determines the extension
  # suffix via sysconfig rather than CMAKE_SHARED_MODULE_SUFFIX.
  execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'), end='')"
    OUTPUT_VARIABLE Python_EXT_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  function(F2PY_BUILD dims prec)

    set( target ${PROJECT_NAME}_${dims}${prec}_f2py )
    set( sofile "${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd_${dims}${prec}${Python_EXT_SUFFIX}" )

    set(lomp)
    set(fflags)

    if( dims STREQUAL "xd" )
      set( blas 1 )
      if(BLAS_FOUND AND LAPACK_FOUND)
        list( APPEND lomp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
        list( APPEND fflags ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS} )
      endif()
    endif()

    if(OpenMP_FOUND)
      foreach(var IN ITEMS ${OpenMP_Fortran_LIB_NAMES})
        list( APPEND lomp "-l${var}" )
      endforeach()
      list( APPEND fflags ${OpenMP_Fortran_FLAGS} )
    endif()

    # Collapse fflags list to a single space-separated string for FFLAGS env var.
    string(JOIN " " fflags_str ${fflags}
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-O3>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-cpp>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-std=f2008>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fimplicit-none>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fall-intrinsics>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","NVHPC">:-fast>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","NVHPC">:-cpp>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","NVHPC">:-Mnobackslash>
      $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:-fast>
      $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:-free>
      $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:-fpp>
      $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:-fast>
      $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:-free>
      $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:-fpp>
      $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:/fast>
      $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:/free>
      $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:/fpp>
      $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:/fast>
      $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:/free>
      $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:/fpp>
    )

    message(STATUS ${dims} ${prec} ${sofile})
    add_custom_command(
      OUTPUT ${sofile}
      # f2py meson backend (Python>=3.12) reads .f2py_f2cmap from CWD; copy it there first.
      COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap_${prec}"
        "${CMAKE_CURRENT_BINARY_DIR}/.f2py_f2cmap"
      # Pass Fortran flags via FFLAGS env var (meson backend ignores --f90flags).
      COMMAND ${CMAKE_COMMAND} -E env "FFLAGS=${fflags_str}"
        ${Python_EXECUTABLE} -m "numpy.f2py"
          -c
          -m mobbrmsd_${dims}${prec}
          $<$<CONFIG:Debug>:--debug-capi>
          "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
          -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd_fortran/${dims}${prec}/mod_files
          -L${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd_fortran/static
          -lmobbrmsd_${dims}${prec}
          ${lomp}
          $<$<CONFIG:Debug>:--debug>
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      DEPENDS "mobbrmsd_${dims}${prec}_static"
              "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap_${prec}"
              "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
      )

    add_custom_target(${target} ALL DEPENDS ${sofile})

    install(FILES
      ${sofile}
      DESTINATION ${SKBUILD_PLATLIB_DIR}/mobbrmsd
    )

  endfunction(F2PY_BUILD)

  foreach(dims IN ITEMS ${mobbrmsd_DIMENTIONS})
    foreach(prec IN ITEMS ${mobbrmsd_PRECISIONS})
      F2PY_BUILD( ${dims} ${prec} )
    endforeach()
  endforeach()

  set(preprocessor ${CMAKE_CURRENT_SOURCE_DIR}/../preprocessor.py)

  function(PREPROC_PYTHON list newlist py)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${py}.py
        COMMAND ${Python_EXECUTABLE}
        ${preprocessor}
        ${CMAKE_CURRENT_SOURCE_DIR}/${py}.py
        ${CMAKE_CURRENT_BINARY_DIR}/${py}.py
      )
    add_custom_target(${py} ALL
      DEPENDS
      ${preprocessor}
      ${CMAKE_CURRENT_SOURCE_DIR}/${py}.py
      ${CMAKE_CURRENT_BINARY_DIR}/${py}.py
      )
    set(list_)
    foreach(t IN LISTS list)
      list(APPEND list_ ${t})
    endforeach()
    list(APPEND list_ ${CMAKE_CURRENT_BINARY_DIR}/${py}.py)
    set(${newlist} ${list_} PARENT_SCOPE)
  endfunction(PREPROC_PYTHON)

  PREPROC_PYTHON("${python_files}" python_files mobbrmsd)
  PREPROC_PYTHON("${python_files}" python_files dataclass)
  PREPROC_PYTHON("${python_files}" python_files __init__)
  PREPROC_PYTHON("${python_files}" python_files __main__)
  install(FILES
    ${python_files}
    DESTINATION ${SKBUILD_PLATLIB_DIR}/mobbrmsd
  )
  add_subdirectory(run)
  add_subdirectory(demo)
endif()

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  ADDTEST(test_driver 3d dp 32 1 1 ${CMAKE_CURRENT_BINARY_DIR} ""
    ${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/test_driver.f90
    )
  target_include_directories(test_driver PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd_fortran/3ddp/mod_files
    )
  target_link_libraries(test_driver PRIVATE
    mobbrmsd_3ddp_static
    )
  if(OpenMP_FOUND)
    target_link_libraries(test_driver PRIVATE OpenMP::OpenMP_Fortran)
  endif()
endif()

if(SKBUILD)
  find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

  # Grab the variables from a local Python installation F2PY headers
  execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  function(F2PY_BUILD dims prec)

    set( target ${PROJECT_NAME}_${dims}${prec}_f2py )

    if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
      set (sofile "${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd_${dims}${prec}.${SKBUILD_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
    else()
      set (sofile "${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd_${dims}${prec}.${SKBUILD_SOABI}.so")
    endif()

    set(lomp)
    set(lflag)

    if( dims STREQUAL "xd" )
      set( blas 1 )
      if(BLAS_FOUND AND LAPACK_FOUND)
        list( APPEND lomp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
        list( APPEND lflag ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS} )
      endif()
    endif()

    if(OpenMP_FOUND)
      foreach(var IN ITEMS ${OpenMP_Fortran_LIB_NAMES})
        list( APPEND lomp "-l${var} " )
      endforeach()
      list( APPEND lflag ${OpenMP_Fortran_FLAGS} )
    endif()

    message(STATUS ${dims} ${prec} ${sofile})
    add_custom_command(
      OUTPUT ${sofile}
      COMMAND ${Python_EXECUTABLE} -m "numpy.f2py"
        -c
        -m mobbrmsd_${dims}${prec}
        --no-wrap-functions
        $<$<CONFIG:Debug>:--debug-capi>
        --f2cmap "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap_${prec}"
        "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
        -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd_fortran/${dims}${prec}/mod_files
        -L${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd_fortran/static
        -lmobbrmsd_${dims}${prec}
        ${lomp}
        --f90flags='
        "${lflag}"
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-O3">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-cpp">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-std=f2008">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fimplicit-none">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fall-intrinsics">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","NVHPC">:"-fast">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","NVHPC">:"-cpp">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","NVHPC">:"-Mnobackslash">
        $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:"-fast">
        $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:"-free">
        $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:"-fpp">
        $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:"-fast">
        $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:"-free">
        $<$<AND:$<BOOL:${UNIX}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:"-fpp">
        $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:"/fast">
        $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:"/free">
        $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">>:"/fpp">
        $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:"/fast">
        $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:"/free">
        $<$<AND:$<BOOL:${WIN32}>,$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">>:"/fpp">
        '
        $<$<CONFIG:Debug>:--debug>
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

  function(PREPROC_PYTHON list newlist py)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${py}.py
        COMMAND ${Python_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/preprocessor.py
        ${CMAKE_CURRENT_SOURCE_DIR}/${py}.py
        > ${CMAKE_CURRENT_BINARY_DIR}/${py}.py 2>&1
      )
    add_custom_target(${py} ALL
      DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/preprocessor.py
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

  PREPROC_PYTHON("${python_files}" python_files _mobbrmsd)
  PREPROC_PYTHON("${python_files}" python_files _demo)
  PREPROC_PYTHON("${python_files}" python_files coord_generator)
  PREPROC_PYTHON("${python_files}" python_files demo_cogen)
  PREPROC_PYTHON("${python_files}" python_files demo_bb)
  PREPROC_PYTHON("${python_files}" python_files demo_bb_2d)
  PREPROC_PYTHON("${python_files}" python_files demo_bb_multi)
  PREPROC_PYTHON("${python_files}" python_files demo_mst)
  PREPROC_PYTHON("${python_files}" python_files demo_batch)
  PREPROC_PYTHON("${python_files}" python_files demo_batch_tri)
  PREPROC_PYTHON("${python_files}" python_files __init__)
  PREPROC_PYTHON("${python_files}" python_files __main__)
  install(FILES
    ${python_files}
    DESTINATION ${SKBUILD_PLATLIB_DIR}/mobbrmsd
  )
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
endif()

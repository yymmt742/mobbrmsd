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
        -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/${dims}${prec}/mod_files
        -L${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/static
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

  install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/__main__.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/_demo.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/_mobbrmsd.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/coord_generator.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_cogen.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_bb.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_bb_multi.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_mst.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_batch.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_batch_tri.py"
    DESTINATION ${SKBUILD_PLATLIB_DIR}/mobbrmsd
    )
endif()

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  ADDTEST(test_driver 3d dp 32 1 1 ${CMAKE_CURRENT_BINARY_DIR} ""
    ${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/test_driver.f90
    )
  target_include_directories(test_driver PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/3ddp/mod_files
    )
  target_link_libraries(test_driver PRIVATE
    mobbrmsd_3ddp_static
    )
endif()

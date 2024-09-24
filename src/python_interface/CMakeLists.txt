if(SKBUILD)
  find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

  # Grab the variables from a local Python installation
  # F2PY headers
  execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  foreach(var IN ITEMS ${OpenMP_Fortran_LIB_NAMES})
    list( APPEND lomp "-l${var} " )
  endforeach()

  function(F2PY_BUILD dim prec)
    if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
      set (sofile "${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd_${dim}${prec}.${SKBUILD_SOABI}.so")
    else()
      set (sofile "${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd_${dim}${prec}.${SKBUILD_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
    endif()
    add_custom_command(
      OUTPUT ${sofile}
      COMMAND ${Python_EXECUTABLE} -m "numpy.f2py"
        -c
        -m mobbrmsd_${dim}${prec}
        --no-wrap-functions
        $<$<CONFIG:Debug>:--debug-capi>
        --f2cmap "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap_${prec}"
        "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
        "${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_${dim}/libmobbrmsd_${dim}${prec}.a"
        -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_${dim}/mobbrmsd_${dim}${prec}/mod_files
        ${BLAS_LIBRARIES}
        ${lomp}
        --opt=-O3
        --f90flags='
        "${BLAS_LINKER_FLAGS}"
        "${OpenMP_Fortran_FLAGS}"
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
      DEPENDS "mobbrmsd_${dim}${prec}_static"
              "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap_${prec}"
              "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
      )
    add_custom_target(mobbrmsd_${dim}${prec}_so ALL DEPENDS ${sofile})
    install(FILES
      ${sofile}
      DESTINATION ${Python_SITELIB}/mobbrmsd
    )
  endfunction(F2PY_BUILD)

  F2PY_BUILD(xd dp)
  F2PY_BUILD(xd sp)
  F2PY_BUILD(2d dp)
  F2PY_BUILD(2d sp)
  F2PY_BUILD(3d dp)
  F2PY_BUILD(3d sp)

  install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/__main__.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/__version__.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/_demo.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/_mobbrmsd.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/coord_generator.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_cogen.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_bb.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_bb_multi.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_mst.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_batch.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_batch_tri.py"
    DESTINATION ${Python_SITELIB}/mobbrmsd
    )
endif()

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  ADDTEST(test_driver 3 64 32 1 1 ""
    ${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/test_driver.f90
    )
  target_include_directories(test_driver PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_3d/mobbrmsd_3ddp/mod_files
    )
  target_link_libraries(test_driver PRIVATE
    mobbrmsd_3ddp_static
    )
endif()

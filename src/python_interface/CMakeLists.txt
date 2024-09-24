find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
if(SKBUILD)
  find_package(Python 3.8 REQUIRED COMPONENTS NumPy)
  find_package(F2PY REQUIRED)

  if (NOT CMAKE_INSTALL_PYTHON_LIBDIR)
    execute_process(
     COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
               from distutils import sysconfig as sc
               print(sc.get_python_lib(prefix='', plat_specific=True))"
     OUTPUT_VARIABLE CMAKE_INSTALL_PYTHON_LIBDIR
     OUTPUT_STRIP_TRAILING_WHITESPACE)
   endif ()

# Grab the variables from a local Python installation
# NumPy headers
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}"
      -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
# F2PY headers
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}"
    -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  include_directories(${PYTHON_INCLUDE_DIR})
  include_directories(${NumPy_INCLUDE_DIRS})

  foreach(var IN ITEMS ${OpenMP_Fortran_LIB_NAMES})
    list( APPEND lomp "-l${var} " )
  endforeach()

  function(FYPP_BUILD dim prec)
    set (sofile "${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd_${dim}${prec}.${Python_SOABI}.so")
    add_custom_command(
      OUTPUT ${sofile}
      COMMAND ${F2PY_EXECUTABLE}
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
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">:"-fast">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">:"-free">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","IntelLLVM">:"-fpp">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-fast">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-free">
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-fpp">
        '
        $<$<CONFIG:Debug>:--debug>
      DEPENDS "mobbrmsd_${dim}${prec}_static"
              "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap_${prec}"
              "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
    )
    add_custom_target(mobbrmsd_${dim}${prec}_so ALL DEPENDS ${sofile})
    install(FILES
      ${sofile}
      DESTINATION
      "${CMAKE_INSTALL_PYTHON_LIBDIR}/mobbrmsd"
    )
  endfunction(FYPP_BUILD)

  FYPP_BUILD(xd dp)
  FYPP_BUILD(xd sp)
  FYPP_BUILD(2d dp)
  FYPP_BUILD(2d sp)
  FYPP_BUILD(3d dp)
  FYPP_BUILD(3d sp)

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
    DESTINATION
    "${CMAKE_INSTALL_PYTHON_LIBDIR}/mobbrmsd"
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

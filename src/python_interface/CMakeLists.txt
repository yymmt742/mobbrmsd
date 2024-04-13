find_package(F2PY REQUIRED)
find_package(Python 3.8 REQUIRED COMPONENTS NumPy)
find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)

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

set (initpy "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")
set (f2cmap "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap")
set (sofile "${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd.${Python_SOABI}.so")
set (sofile_2d "${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd_2d.${Python_SOABI}.so")
set (sofile_3d "${CMAKE_CURRENT_BINARY_DIR}/mobbrmsd_3d.${Python_SOABI}.so")

add_custom_target( mobbrmsd_ ALL DEPENDS ${sofile} ${sofile_2d} ${sofile_3d} )

include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${NumPy_INCLUDE_DIRS})

foreach(var IN ITEMS ${OpenMP_Fortran_LIB_NAMES})
  list( APPEND lomp "-l${var} " )
endforeach()

add_custom_command(
    OUTPUT ${sofile}
    COMMAND ${F2PY_EXECUTABLE}
      -c
      -m mobbrmsd
      --no-wrap-functions
      $<$<CONFIG:Debug>:--debug-capi>
      --f2cmap ${f2cmap}
      "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
      "${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd/libmobbrmsd.a"
      -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd/mod_files
      ${BLAS_LIBRARIES}
      ${lomp}
      --opt=-O2
      --f90flags='
      "${BLAS_LINKER_FLAGS}"
      "${OpenMP_Fortran_FLAGS}"
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-std=f2008">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fimplicit-none">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fall-intrinsics">
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-Wuninitialized">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fbounds-check">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-Wunused">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fbacktrace">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fcheck=array-temps,bounds,do,mem,pointer,recursion">>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-stand">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"f08">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-free">
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-warn">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:" all">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-check">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"bounds">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-traceback">>
      '
      $<$<CONFIG:Debug>:--debug>
    DEPENDS mobbrmsd
            "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap"
            "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
  )

add_custom_command(
    OUTPUT ${sofile_2d}
    COMMAND ${F2PY_EXECUTABLE}
      -c
      -m mobbrmsd_2d
      --no-wrap-functions
      $<$<CONFIG:Debug>:--debug-capi>
      --f2cmap ${f2cmap}
      "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
      "${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_2d/libmobbrmsd_2d.a"
      -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_2d/mod_files
      ${BLAS_LIBRARIES}
      ${lomp}
      --opt=-O2
      --f90flags='
      "${BLAS_LINKER_FLAGS}"
      "${OpenMP_Fortran_FLAGS}"
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-std=f2008">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fimplicit-none">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fall-intrinsics">
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-Wuninitialized">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fbounds-check">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-Wunused">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fbacktrace">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fcheck=array-temps,bounds,do,mem,pointer,recursion">>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-stand">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"f08">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-free">
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-warn">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:" all">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-check">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"bounds">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-traceback">>
      '
      $<$<CONFIG:Debug>:--debug>
    DEPENDS mobbrmsd_2d
            "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap"
            "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
  )

add_custom_command(
    OUTPUT ${sofile_3d}
    COMMAND ${F2PY_EXECUTABLE}
      -c
      -m mobbrmsd_3d
      --no-wrap-functions
      $<$<CONFIG:Debug>:--debug-capi>
      --f2cmap ${f2cmap}
      "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
      "${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_3d/libmobbrmsd_3d.a"
      -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_3d/mod_files
      ${BLAS_LIBRARIES}
      ${lomp}
      --opt=-O2
      --f90flags='
      "${BLAS_LINKER_FLAGS}"
      "${OpenMP_Fortran_FLAGS}"
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-std=f2008">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fimplicit-none">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fall-intrinsics">
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-Wuninitialized">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fbounds-check">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-Wunused">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fbacktrace">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:"-fcheck=array-temps,bounds,do,mem,pointer,recursion">>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-stand">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"f08">
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-free">
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-warn">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:" all">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-check">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"bounds">>
      $<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:"-traceback">>
      '
      $<$<CONFIG:Debug>:--debug>
    DEPENDS mobbrmsd_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap"
            "${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90"
  )

install(FILES
    ${sofile}
    ${sofile_2d}
    ${sofile_3d}
    "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/__main__.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/__mobbrmsd__.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/__coord_generator__.py"
  DESTINATION
    "${CMAKE_INSTALL_PYTHON_LIBDIR}/mobbrmsd")

if( CMAKE_BUILD_TYPE MATCHES "Debug" )

  add_compile_options(
    ${OpenMP_Fortran_FLAGS}
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-O3>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fimplicit-none>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fall-intrinsics>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fexternal-blas>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wuninitialized>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbounds-check>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wunused>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbacktrace>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fcheck=array-temps,bounds,do,mem,pointer>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-O3>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fpp>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-stand>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:f08>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-free>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-implicitnone>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-qopt-matmul>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-warn>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:all>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-check>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:bounds>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-traceback>
  )

  find_package(unittest REQUIRED)
  add_compile_options(${BLAS_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS})
  add_compile_options(${BLAS_LINKER_FLAGS})

  add_executable(test_driver)
  target_sources(test_driver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/test_driver.f90)
  target_include_directories(test_driver PRIVATE ${unittest_DIR}/mod_files)
  target_include_directories(test_driver PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_3d/mod_files)
  target_include_directories(test_driver PRIVATE ${OpenMP_Fortran_INCLUDE_DIRS} )
  target_link_libraries(test_driver PRIVATE blas )
  target_link_libraries(test_driver PRIVATE ${OpenMP_Fortran_LIB_NAMES} )
  target_link_libraries(test_driver PRIVATE mobbrmsd_3d_shared)
  target_link_libraries(test_driver PRIVATE unittest)
  add_test( NAME test_driver COMMAND test_driver )
endif()

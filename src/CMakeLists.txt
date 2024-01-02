find_package( BLAS REQUIRED )

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

include_directories( ${CMAKE_CURRENT_BINARY_DIR} )

add_compile_options(
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-O2>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fimplicit-none>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fall-intrinsics>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fexternal-blas>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-O2>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-stand>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:f08>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-free>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-implicitnone>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-qopt-matmul>
  )

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_compile_options(
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wuninitialized>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbounds-check>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wunused>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbacktrace>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fcheck=array-temps,bounds,do,mem,pointer,recursion>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-warn>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:all>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-check>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:bounds>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-traceback>
  )
endif()

add_subdirectory( symRMSD )


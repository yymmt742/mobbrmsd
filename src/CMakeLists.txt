find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)

if(BLAS_FOUND)
  add_compile_options(${BLAS_LINKER_FLAGS})
endif()

add_compile_options(
  ${OpenMP_Fortran_FLAGS}
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
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fcheck=array-temps,bounds,do,mem,pointer>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-warn>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:all>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-check>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:bounds>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-traceback>
  )
endif()

add_subdirectory( mobbrmsd )

if(SKBUILD)
  add_subdirectory( python_interface )
endif()

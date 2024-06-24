find_package(OpenMP REQUIRED)

add_compile_options(
  ${OpenMP_Fortran_FLAGS}
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-O3>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fimplicit-none>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fall-intrinsics>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-g>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wuninitialized>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbounds-check>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wunused>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbacktrace>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fcheck=array-temps,bounds,do,mem,pointer>>
  $<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-fast>
  $<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-cpp>
  $<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mnobackslash>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-C>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-g>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-traceback>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Wall>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mchkptr>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mchkstk>>
  )

if(UNIX)
  add_compile_options(
  $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-O3>
  $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-free>
  $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fpp>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-g>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-traceback>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-debug>>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-O3>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-free>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fpp>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-g>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-traceback>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-debug>>
  )
endif()

add_subdirectory( mobbrmsd )
add_subdirectory( python_interface )


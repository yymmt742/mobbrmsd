find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  find_package(unittest)
  if(NOT unittest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      unittest
      GIT_REPOSITORY https://github.com/yymmt742/unittest
      GIT_TAG origin/main
      GIT_PROGRESS TRUE
      )
    FetchContent_MakeAvailable(unittest)
    set(unittest_INCLUDE_DIRS ${unittest_BINARY_DIR}/src/mod_files)
  else()
    set(unittest_INCLUDE_DIRS ${unittest_DIR}/mod_files)
  endif()
endif()

function(ADDOPTS name real int blas openmp)
  target_compile_options(
    ${name} PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-O3>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fimplicit-none>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fall-intrinsics>
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
  target_compile_options(
    ${name} PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fast>
    $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-free>
    $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fpp>
    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-g>>
    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-traceback>>
    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-debug>>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fast>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-free>
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fpp>
    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-g>>
    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-traceback>>
    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-debug>>
  )
endif()
  target_compile_definitions(${name} PRIVATE -DUSE_REAL${real})
  target_compile_definitions(${name} PRIVATE -DUSE_INT${int})
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  if(blas)
    target_compile_options(${name} PRIVATE ${BLAS_LINKER_FLAGS})
    target_link_libraries(${name} PRIVATE blas)
    target_link_libraries(${name} PRIVATE lapack)
  endif()
  if(openmp)
    target_compile_options(${name} PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_libraries(${name} PRIVATE ${OpenMP_Fortran_LIB_NAMES})
  endif()
endfunction(ADDOPTS)

function(ADDTEST name real int blas openmp testarg)
  add_executable(${name} ${ARGN})
  ADDOPTS("${name}" "${real}" "${int}" "${blas}" "${openmp}")
  target_link_libraries(${name} PRIVATE unittest)
  target_include_directories(${name} PRIVATE ${unittest_INCLUDE_DIRS})
  add_test(NAME ${name} COMMAND ${name} ${testarg})
endfunction(ADDTEST)

#get_cmake_property(_variableNames CACHE_VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()


#add_compile_options(
#  #  ${OpenMP_Fortran_FLAGS}
#  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-O3>
#  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>
#  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>
#  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fimplicit-none>
#  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fall-intrinsics>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wuninitialized>>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbounds-check>>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wunused>>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbacktrace>>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fcheck=array-temps,bounds,do,mem,pointer>>
#  $<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-fast>
#  $<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-cpp>
#  $<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mnobackslash>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-C>>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-g>>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-traceback>>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Wall>>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mchkptr>>
#  $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mchkstk>>
#  )
#if(UNIX)
#  add_compile_options(
#    $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fast>
#    $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-free>
#    $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fpp>
#    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-g>>
#    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-traceback>>
#    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-debug>>
#    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fast>
#    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-free>
#    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fpp>
#    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-g>>
#    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-traceback>>
#    $<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-debug>>
#  )
#endif()

add_subdirectory( mobbrmsd )
add_subdirectory( python_interface )


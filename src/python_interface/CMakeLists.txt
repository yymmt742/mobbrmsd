find_package(F2PY REQUIRED)
find_package(Python 3.8 REQUIRED COMPONENTS NumPy)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)

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

add_custom_target( mobbrmsd_ ALL DEPENDS ${sofile} )

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
      "${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/libmobbrmsd.a"
      -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mod_files
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

install( FILES ${sofile} DESTINATION mobbrmsd )
install( FILES ${initpy} DESTINATION mobbrmsd )

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  find_package(unittest REQUIRED)
  add_compile_options(${BLAS_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS})

  add_executable(test_driver)
  target_sources(test_driver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/test_driver.f90)
  target_include_directories(test_driver PRIVATE ${unittest_DIR}/mod_files)
  target_include_directories(test_driver PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mod_files)
  target_include_directories(test_driver PRIVATE ${OpenMP_Fortran_INCLUDE_DIRS} )
  target_link_libraries(test_driver PRIVATE ${OpenMP_Fortran_LIB_NAMES} )
  target_link_libraries(test_driver PRIVATE mobbrmsd_shared)
  target_link_libraries(test_driver PRIVATE unittest)
  add_test( NAME test_driver COMMAND test_driver )
endif()

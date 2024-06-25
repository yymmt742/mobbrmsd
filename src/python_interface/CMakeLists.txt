list(APPEND f90_flags
    "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-O3>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fimplicit-none>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fall-intrinsics>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wuninitialized>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbounds-check>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wunused>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fbacktrace>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fcheck=array-temps,bounds,do,mem,pointer>>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-fast>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-cpp>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mnobackslash>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-C>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-g>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-traceback>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Wall>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mchkptr>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mchkstk>>"
)
if(UNIX)
  list(APPEND f90_flags
    "$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fast>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-free>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fpp>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-g>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-traceback>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-debug>>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fast>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-free>"
    "$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fpp>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-g>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-traceback>>"
    "$<$<CONFIG:Debug>:$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-debug>>"
  )
endif()
if(SKBUILD)
  find_package(F2PY REQUIRED)
  find_package(Python 3.8 REQUIRED COMPONENTS NumPy)

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
      "${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd/libmobbrmsd_gddp.a"
      -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_gd/mobbrmsd_gddp/mod_files
      ${BLAS_LIBRARIES}
      ${lomp}
      --opt=-O3
      --f90flags='
      "${BLAS_LINKER_FLAGS}"
      "${OpenMP_Fortran_FLAGS}"
      "${f90_flags}"
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
      "${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_2d/libmobbrmsd_2ddp.a"
      -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_2d/mobbrmsd_2ddp/mod_files
      ${lomp}
      --opt=-O3
      --f90flags='
      "${BLAS_LINKER_FLAGS}"
      "${OpenMP_Fortran_FLAGS}"
      "${f90_flags}"
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
      "${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_3d/libmobbrmsd_3ddp.a"
      -I${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_3d/mobbrmsd_3ddp/mod_files
      ${lomp}
      --opt=-O3
      --f90flags='
      "${BLAS_LINKER_FLAGS}"
      "${OpenMP_Fortran_FLAGS}"
      "${f90_flags}"
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
    "${CMAKE_CURRENT_SOURCE_DIR}/__version__.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/_mobbrmsd.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/coord_generator.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_cogen.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_bb.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_mst.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_batch.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/demo_batch_tri.py"
    DESTINATION
    "${CMAKE_INSTALL_PYTHON_LIBDIR}/mobbrmsd")

endif()

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_compile_options(${BLAS_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS})
  add_compile_options(${BLAS_LINKER_FLAGS})

  add_executable(test_driver)
  target_sources(test_driver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/python_driver.f90
                                     ${CMAKE_CURRENT_SOURCE_DIR}/test_driver.f90)
  target_include_directories(test_driver PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../mobbrmsd/mobbrmsd_3d/mobbrmsd_3ddp/mod_files)
  target_include_directories(test_driver PRIVATE ${OpenMP_Fortran_INCLUDE_DIRS} )
  target_include_directories(test_driver PRIVATE ${unittest_INCLUDE_DIRS})
  target_link_libraries(test_driver PRIVATE blas )
  target_link_libraries(test_driver PRIVATE ${OpenMP_Fortran_LIB_NAMES} )
  target_link_libraries(test_driver PRIVATE mobbrmsd_3ddp_shared)
  target_link_libraries(test_driver PRIVATE unittest)
  add_test( NAME test_driver COMMAND test_driver )
endif()

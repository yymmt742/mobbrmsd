find_package( F2PY REQUIRED )
find_package( Python 3.8 REQUIRED COMPONENTS NumPy )

set (initpy "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")
set (f2cmap "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap")
set (sofile "${CMAKE_CURRENT_BINARY_DIR}/symrmsd.${Python_SOABI}.so")

add_custom_target( symrmsd ALL
    DEPENDS ${sofile}
)

include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${F2PY_INCLUDE_DIR})
include_directories(${NumPy_INCLUDE_DIRS})

add_custom_command(
    OUTPUT ${sofile}
    COMMAND ${F2PY_EXECUTABLE}
      -c
      -m symrmsd
      --no-wrap-functions
      $<$<CONFIG:Debug>:--debug-capi>
      --f2cmap ${f2cmap}
      "${CMAKE_CURRENT_SOURCE_DIR}/symrmsd.pyf"
      -I${CMAKE_CURRENT_BINARY_DIR}/..
      -L${CMAKE_CURRENT_BINARY_DIR}/..
      -lsymRMSD
      -lblas
      --opt=-O2
      --f90flags='
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
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS symrmsd
            "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap"
            "${CMAKE_CURRENT_SOURCE_DIR}/symrmsd.pyf"
  )

install( FILES ${sofile} DESTINATION symrmsd )
install( FILES ${initpy} DESTINATION symrmsd )


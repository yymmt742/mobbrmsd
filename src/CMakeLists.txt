find_package( BLAS REQUIRED )
find_package( purelapack REQUIRED )

add_compile_options(
  $<$<COMPILE_LANG_AND_ID:C,GNU>:-O2>
  $<$<COMPILE_LANG_AND_ID:C,GNU>:-Wall>
  $<$<COMPILE_LANG_AND_ID:C,GNU>:-Wextra>
  )

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
    $<$<COMPILE_LANG_AND_ID:C,GNU>:-g>
    $<$<COMPILE_LANG_AND_ID:C,GNU>:-Wshadow>
    $<$<COMPILE_LANG_AND_ID:C,GNU>:-Wconversion>
    $<$<COMPILE_LANG_AND_ID:C,GNU>:-Wfloat-equal>
    $<$<COMPILE_LANG_AND_ID:C,GNU>:-ftrapv>
    $<$<COMPILE_LANG_AND_ID:C,GNU>:-fstack-protector-all>
    )
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

add_subdirectory( molutil )

add_library( symRMSD )
target_include_directories( symRMSD PRIVATE ${purelapack_DIR}/mod_files )
target_include_directories( symRMSD PRIVATE ${purelapack_DIR}/include )

list( APPEND f90_src mod_params.f90 )
list( APPEND f90_src mod_optarg.f90 )
list( APPEND f90_src mod_permutation.f90 )
list( APPEND f90_src mod_group_permutation.f90 )
list( APPEND f90_src mod_molecular_rotation.f90 )
list( APPEND f90_src mod_molecular_permutation.f90 )
list( APPEND f90_src mod_rmsd.f90 )
list( APPEND f90_src mod_rmsd_brute.f90 )
list( APPEND f90_src mod_det.f90 )
list( APPEND f90_src mod_svd.f90 )
list( APPEND f90_src mod_Kabsch.f90 )
list( APPEND f90_src mod_mol_block.f90 )
list( APPEND f90_src mod_Hungarian.f90 )
list( APPEND f90_src mod_d_matrix.f90 )
list( APPEND f90_src mod_tree.f90 )
#list( APPEND f90_src mod_symRMSD.f90 )

target_sources( symRMSD PRIVATE ${f90_src} )
target_link_libraries( symRMSD PRIVATE molutil )
target_link_libraries( symRMSD PRIVATE blas )
target_link_libraries( symRMSD PRIVATE purelapack )
#target_link_libraries ( symRMSD jsonfortran::jsonfortran-static )

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_subdirectory( test )
endif()


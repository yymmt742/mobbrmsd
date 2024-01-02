add_library( symRMSD )
target_include_directories( symRMSD PRIVATE include )

list( APPEND f90_src mod_params.f90 )
list( APPEND f90_src mod_group_permutation.f90 )
list( APPEND f90_src mod_mol_symmetry.f90 )
list( APPEND f90_src mod_det.f90 )
list( APPEND f90_src mod_estimate_rotation_matrix.f90 )
list( APPEND f90_src mod_mol_block.f90 )
list( APPEND f90_src mod_Hungarian.f90 )
list( APPEND f90_src mod_d_matrix.f90 )
list( APPEND f90_src mod_tree.f90 )
list( APPEND f90_src mod_branch_and_prune.f90 )
list( APPEND f90_src mod_symRMSD.f90 )
list( APPEND f90_src python_driver.f90 )

target_sources( symRMSD PRIVATE ${f90_src} )
target_link_libraries( symRMSD PRIVATE blas )

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_subdirectory( test )
endif()

add_subdirectory( python_interface )

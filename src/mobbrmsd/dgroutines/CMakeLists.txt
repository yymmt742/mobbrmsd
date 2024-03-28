list( APPEND f90_gd_src blas_lapack_interface.f90 )
list( APPEND f90_gd_src mod_rotation.f90 )
list( APPEND f90_gd_src mod_bb_interface_g.f90 )
target_sources(${TARGET_GD} PRIVATE ${f90_gd_src})

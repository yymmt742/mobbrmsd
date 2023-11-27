find_package( unittest REQUIRED )

add_executable( test_molecule )
target_include_directories( test_molecule PRIVATE ${unittest_DIR}/mod_files )
target_sources( test_molecule PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../mod_params.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_optarg.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_base.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_element.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_element_set.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_bonds.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_molecule.f90
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../mod_molecular_topology.f90
                                      test_molecule.f90 )
target_link_libraries( test_molecule PRIVATE unittest )
add_test( NAME test_molecule COMMAND test_molecule )

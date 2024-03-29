!| Module for handling C, F and tree.
module mod_bb_interface_d3
  use mod_bb_block
  use mod_bb_list
  implicit none
  private
!
  public :: bb_block
!
  public :: bb_list
  public :: bb_list_memsize
  public :: bb_list_n_block
  public :: bb_list_n_atoms
  public :: bb_list_log_n_nodes
  public :: bb_list_setup
  public :: bb_list_run
  public :: bb_list_swap_y
  public :: bb_list_rotation_matrix
!
end module mod_bb_interface_d3


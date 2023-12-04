module mod_branch_and_prune
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_molecular_permutation
  use mod_tree
  implicit none
  private
  public :: branch_and_prune
!
  type branch_and_prune
    type(mol_block_list) :: blk
  contains
  end type branch_and_prune
!
contains
!
!| generate node instance
  pure function branch_and_prune_new() result(res)
    class(molecule), intent(in)              :: mol
    !! molecular template
    class(molecular_permutation), intent(in) :: prm
    !! molecular permutation
    real(RK), intent(in)                     :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)                     :: y(*)
    !! target molecular coordinate, y(d,m,n)
    type(node)                               :: res
    integer(IK)                              :: d, m, n, mn, dmn
!
  end function branch_and_prune
!
end module mod_branch_and_prune

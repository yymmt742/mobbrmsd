module driver
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_symRMSD
  use mod_mol_block
  use mod_mol_symmetry
  implicit none
  private
  public add_molecule
  public clear
!
  type(symRMSD_input), save :: inp
  type(symRMSD), save       :: sym
!
contains
!
  subroutine add_molecule(m, n, s, sym)
    integer(IK), intent(in) :: m
    integer(IK), intent(in) :: n
    integer(IK), intent(in) :: s
    integer(IK), intent(in) :: sym(*)
    type(mol_block)         :: b_
    b_ = mol_block(1, s, m, n, m, n)
    call inp%add_molecule(b_, sym)
  end subroutine add_molecule
!
  subroutine clear()
    call inp%clear()
    call sym%clear()
  end subroutine clear
!
end module driver

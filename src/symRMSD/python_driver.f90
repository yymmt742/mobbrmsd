module driver
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_symRMSD
  use mod_mol_block
  use mod_mol_symmetry
  implicit none
  private
  public njob
  public swap_y
  public add_molecule
  public setup
  public run
  public clear
!
  integer(IK), save         :: njob = 1
  logical, save             :: swap_y = .TRUE.
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
print*,inp%blk%d
print*,inp%blk%mn
print*,inp%blk%mg
  end subroutine add_molecule
!
  subroutine setup()
print'(6i4)', inp%blk%b
print*,'njob=', njob
    sym = symRMSD(inp, njob)
  end subroutine setup
!
  subroutine run(x, y, res)
    real(RK), intent(in)    :: x(*)
    real(RK), intent(inout) :: y(*)
    real(RK), intent(out)   :: res
    call sym%run(1, swap_y, x, y, res)
  end subroutine run
!
  subroutine clear()
    call inp%clear()
    call sym%clear()
  end subroutine clear
!
end module driver

module driver
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_symRMSD
  use mod_mol_block
  use mod_mol_symmetry
  implicit none
  private
  public add_molecule
!
  type(symRMSD_input) :: inp
  type(symRMSD)       :: sym
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
! subroutine load(arg, ios)
!   integer(IK)              :: unit, nn(2)
!
!   open (NEWUNIT=unit, FILE=arg, STATUS='OLD', IOSTAT=ios)
!   if (ios > 0) stop
!   read (unit, *, IOSTAT=ios) d
!   if (ios > 0) stop
!   close (unit, IOSTAT=ios)
!
!   call ip%copy(d%ip)
!   n = d%nstate()
!   nn = [n, n]
!   ip = d%ip%index_pointer()
!
! end subroutine load
!
end module driver

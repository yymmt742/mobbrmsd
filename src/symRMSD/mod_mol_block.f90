!| Module for molecular coodinate block indicator.<br>
!  Coordinates must be stored in the following format.<br>
!    X(d,m,n)<br>
!    - d :: spatial dimension.<br>
!    - m :: number of atom in a molecule.<br>
!    - n :: number of molecule.
module mod_mol_block
  use mod_params, only: D, DD, IK, RK
  implicit none
  private
  public :: mol_block
  public :: mol_block_list_init
!
!| molecular block
  type mol_block
    sequence
    !| p :: pointer to X.
    integer(IK)       :: x = 1
    !| s :: number of molecular symmetry
    integer(IK)       :: s = 1
    !| m :: number of atom in a molecule
    integer(IK)       :: m = 1
    !| n :: number of molecule
    integer(IK)       :: n = 1
  end type mol_block
!
contains
!
! Initializer of mol_block array.
  pure subroutine mol_block_list_init(b)
    type(mol_block), intent(inout) :: b(:)
    integer(IK)                    :: i, x
    x = 1
    do i = 1, SIZE(b)
      b(i)%x = x
      x = x + D * b(i)%n * b(i)%m
    end do
  end subroutine mol_block_list_init
!
end module mod_mol_block


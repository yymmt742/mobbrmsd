module mod_symRMSD
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_branch_and_prune
  use mod_mol_block
  use mod_mol_symmetry
  implicit none
  public :: symRMSD
!
  type symRMSD_input
    type(mol_block_list)             :: blk
    type(mol_symmetry), allocatable  :: ms(:)
  contains
    procedure :: clear      => symRMSD_input_clear
    final     :: symRMSD_input_destroy
  end type symRMSD_input
!
  type symRMSD
    private
    type(branch_and_prune) :: bra
    real(RK), allocatable  :: w(:)
  contains
    procedure :: run        => symRMSD_run
    procedure :: clear      => symRMSD_clear
    final     :: symRMSD_destroy
  end type symRMSD
!
contains
!
  pure elemental subroutine symRMSD_input_clear(this)
    class(symRMSD_input), intent(inout) :: this
    call this%blk%clear()
    if (ALLOCATED(this%ms)) deallocate (this%ms)
  end subroutine symRMSD_input_clear
!
  pure elemental subroutine symRMSD_input_destroy(this)
    type(symRMSD_input), intent(inout) :: this
    call symRMSD_input_clear(this)
  end subroutine symRMSD_input_destroy
!
!!!
!
  pure function symRMSD_new(inp) result(res)
    type(symRMSD_input), intent(in) :: inp
    type(symRMSD)                   :: res
    integer(IK)                     :: nmem
!
    if (ALLOCATED(inp%ms)) then
      res%bra = branch_and_prune(inp%blk, 1, inp%ms)
    else
      res%bra = branch_and_prune(inp%blk, 1)
    end if
    nmem = res%bra%memsize()
    allocate (res%w(nmem))
!
  end function symRMSD_new
!
  pure subroutine symRMSD_run(this, x, y, res)
    class(symRMSD), intent(inout) :: this
    real(RK), intent(in)          :: x(*), y(*)
    real(RK), intent(inout)       :: res
!
    call this%bra%setup(x, y, this%w)
!
    call this%bra%run(this%w)
!
    res = this%bra%upperbound(this%W)
!
  end subroutine symRMSD_run
!
  pure elemental subroutine symRMSD_clear(this)
    class(symRMSD), intent(inout) :: this
    call this%bra%clear()
    if (ALLOCATED(this%w)) deallocate (this%w)
  end subroutine symRMSD_clear
!
  pure elemental subroutine symRMSD_destroy(this)
    type(symRMSD), intent(inout) :: this
    call symRMSD_clear(this)
  end subroutine symRMSD_destroy
!
end module mod_symRMSD

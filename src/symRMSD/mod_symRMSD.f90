module mod_symRMSD
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_branch_and_prune
  use mod_mol_block
  use mod_mol_symmetry
  implicit none
  public :: symRMSD
  public :: symRMSD_input
!
  type symRMSD_input
    type(mol_block_list)             :: blk
    type(mol_symmetry), allocatable  :: ms(:)
  contains
    procedure :: add_molecule => symRMSD_input_add_molecule
    procedure :: clear        => symRMSD_input_clear
    final     :: symRMSD_input_destroy
  end type symRMSD_input
!
  type symRMSD
    private
    integer(IK)            :: njob = 1
    type(branch_and_prune) :: bra
    real(RK), allocatable  :: w(:, :)
  contains
    procedure :: run        => symRMSD_run
    procedure :: lowerbound => symRMSD_lowerbound
    procedure :: upperbound => symRMSD_upperbound
    procedure :: clear      => symRMSD_clear
    final     :: symRMSD_destroy
  end type symRMSD
!
  interface symRMSD
    module procedure symRMSD_new
  end interface symRMSD
!
contains
!
  pure subroutine symRMSD_input_add_molecule(this, b, s)
    class(symRMSD_input), intent(inout) :: this
    type(mol_block), intent(in)         :: b
    integer(IK), intent(in)             :: s(*)
    type(mol_symmetry), allocatable     :: ms(:)
    integer(IK)                         :: l, i, n, h(2)
!
    call this%blk%add_molecule(b)
!
    l = this%blk%nspecies()
    h(1) = this%blk%b(l)%m
    h(2) = this%blk%b(l)%s
    n = h(1) * h(2)
!
    allocate (ms(l))
    do concurrent(i=1:l - 1)
      ms(i) = this%ms(i)
    end do
    ms(l) = mol_symmetry(RESHAPE(s(:n), h))
    call move_alloc(from=ms, to=this%ms)
!
  end subroutine symRMSD_input_add_molecule

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
  pure function symRMSD_new(inp, njob) result(res)
    type(symRMSD_input), intent(in) :: inp
    integer(IK), intent(in)         :: njob
    type(symRMSD)                   :: res
    integer(IK)                     :: nmem
!
    if (ALLOCATED(inp%ms)) then
      res%bra = branch_and_prune(inp%blk, 1, inp%ms)
    else
      res%bra = branch_and_prune(inp%blk, 1)
    end if
!
    nmem = res%bra%memsize()
    res%njob = MAX(njob, 1)
    allocate (res%w(nmem, res%njob))
!
  end function symRMSD_new
!
  pure subroutine symRMSD_run(this, ijob, yswap, x, y, res)
    class(symRMSD), intent(inout) :: this
    integer(IK), intent(in)       :: ijob
    logical, intent(in)           :: yswap
    real(RK), intent(in)          :: x(*)
    real(RK), intent(inout)       :: y(*), res
!
    if(ijob<1.or.this%njob<ijob) return
!
    call this%bra%setup(x, y, this%w)
!
    call this%bra%run(this%w)
!
    res = this%bra%upperbound(this%W)
!
    if(yswap) call this%bra%swap(y)
!
  end subroutine symRMSD_run
!
  pure elemental function symRMSD_lowerbound(this) result(res)
    class(symRMSD), intent(in) :: this
    real(RK)                   :: res
    res = this%bra%lowerbound(this%W)
  end function symRMSD_lowerbound
!
  pure elemental function symRMSD_upperbound(this) result(res)
    class(symRMSD), intent(in) :: this
    real(RK)                   :: res
    res = this%bra%upperbound(this%W)
  end function symRMSD_upperbound
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

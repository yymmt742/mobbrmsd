module mod_mol_symmetry
  use mod_params, only: IK, RK
  use mod_group_permutation
  use mod_group_permutation
  implicit none
  private
  public :: mol_symmetry
!
  type mol_symmetry
    private
    integer(IK)                          :: m = 0
    type(group_permutation), allocatable :: p(:)
  contains
    procedure :: n_atom       => mol_symmetry_n_atom
    procedure :: n_sym        => mol_symmetry_n_sym
    procedure :: swap         => mol_symmetry_swap
    procedure :: reverse      => mol_symmetry_reverse
    procedure :: clear        => mol_symmetry_clear
    final     :: mol_symmetry_destroy
  end type mol_symmetry
!
  interface mol_symmetry
    module procedure mol_symmetry_new
  end interface mol_symmetry
!
contains
!
!| Constructer
  pure function mol_symmetry_new(sym) result(res)
    integer(IK), intent(in), optional :: sym(:, :)
    !! symmetry indices
    type(mol_symmetry)                :: res
    integer(IK)                       :: i, n
!
    if (PRESENT(sym)) then
      res%m = SIZE(sym, 1)
      n = SIZE(sym, 2)
      ALLOCATE(res%p(n))
      do concurrent(i=1:n)
        res%p(i) = group_permutation(sym(:, i))
      enddo
    end if
!
  end function mol_symmetry_new
!
  pure elemental function mol_symmetry_n_atom(this) result(res)
    class(mol_symmetry), intent(in) :: this
    integer(IK)                           :: res
    res = this%m
  end function mol_symmetry_n_atom
!
  pure elemental function mol_symmetry_n_sym(this) result(res)
    class(mol_symmetry), intent(in) :: this
    integer(IK)                           :: res
    if (ALLOCATED(this%p)) then
      res = SIZE(this%p)
    else
      res = 0
    end if
  end function mol_symmetry_n_sym
!
  pure subroutine mol_symmetry_swap(this, d, X, isym)
    class(mol_symmetry), intent(in) :: this
    integer(IK), intent(in)               :: d, isym
    real(RK), intent(inout)               :: X(*)
    if (isym < 1 .or. this%n_sym() < isym) return
    call this%p(isym)%swap(d, X)
  end subroutine mol_symmetry_swap
!
  pure subroutine mol_symmetry_reverse(this, d, X, isym)
    class(mol_symmetry), intent(in) :: this
    integer(IK), intent(in)               :: d, isym
    real(RK), intent(inout)               :: X(*)
    if (isym < 1 .or. this%n_sym() < isym) return
    call this%p(isym)%reverse(d, X)
  end subroutine mol_symmetry_reverse
!
  pure elemental subroutine mol_symmetry_clear(this)
    class(mol_symmetry), intent(inout) :: this
    this%m = 0
    if (ALLOCATED(this%p)) deallocate (this%p)
  end subroutine mol_symmetry_clear
!
  pure elemental subroutine mol_symmetry_destroy(this)
    type(mol_symmetry), intent(inout) :: this
    call this%clear()
  end subroutine mol_symmetry_destroy
!
end module mod_mol_symmetry

module mod_molecular_rotation
  use mod_params, only: IK, RK
  use mod_group_permutation
  use mod_optarg
  use mod_group_permutation
  implicit none
  private
  public :: molecular_rotation
!
  type molecular_rotation
    private
    integer(IK)                          :: m = 0
    type(group_permutation), allocatable :: p(:)
  contains
    procedure :: n_atom       => molecular_rotation_n_atom
    procedure :: n_sym        => molecular_rotation_n_sym
    procedure :: swap         => molecular_rotation_swap
    procedure :: reverse      => molecular_rotation_reverse
    procedure :: clear        => molecular_rotation_clear
    final     :: molecular_rotation_destroy
  end type molecular_rotation
!
  interface molecular_rotation
    module procedure molecular_rotation_new
  end interface molecular_rotation
!
contains
!
!| Constructer
  pure function molecular_rotation_new(sym) result(res)
    integer(IK), intent(in), optional :: sym(:, :)
    !! symmetry indices
    type(molecular_rotation)          :: res
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
  end function molecular_rotation_new
!
  pure elemental function molecular_rotation_n_atom(this) result(res)
    class(molecular_rotation), intent(in) :: this
    integer(IK)                           :: res
    res = this%m
  end function molecular_rotation_n_atom
!
  pure elemental function molecular_rotation_n_sym(this) result(res)
    class(molecular_rotation), intent(in) :: this
    integer(IK)                           :: res
    if (ALLOCATED(this%p)) then
      res = SIZE(this%p)
    else
      res = 0
    end if
  end function molecular_rotation_n_sym
!
  pure subroutine molecular_rotation_swap(this, d, X, isym)
    class(molecular_rotation), intent(in) :: this
    integer(IK), intent(in)               :: d, isym
    real(RK), intent(inout)               :: X(*)
    if (isym < 1 .or. this%n_sym() < isym) return
    call this%p(isym)%swap(d, X)
  end subroutine molecular_rotation_swap
!
  pure subroutine molecular_rotation_reverse(this, d, X, isym)
    class(molecular_rotation), intent(in) :: this
    integer(IK), intent(in)               :: d, isym
    real(RK), intent(inout)               :: X(*)
    if (isym < 1 .or. this%n_sym() < isym) return
    call this%p(isym)%reverse(d, X)
  end subroutine molecular_rotation_reverse
!
  pure elemental subroutine molecular_rotation_clear(this)
    class(molecular_rotation), intent(inout) :: this
    this%m = 0
    if (ALLOCATED(this%p)) deallocate (this%p)
  end subroutine molecular_rotation_clear
!
  pure elemental subroutine molecular_rotation_destroy(this)
    type(molecular_rotation), intent(inout) :: this
    call this%clear()
  end subroutine molecular_rotation_destroy
!
end module mod_molecular_rotation

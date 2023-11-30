module mod_molecular_rotation
  use mod_params, only: IK, RK
  use mod_group_permutation
  use mod_optarg
  implicit none
  private
  public :: molecular_rotation
!
  type molecular_rotation
    private
    type(group_permutation), allocatable :: p(:)
  contains
    procedure :: nsym         => molecular_rotation_nsym
    procedure :: swap         => molecular_rotation_swap
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
    !! swap indices
    type(molecular_rotation)          :: res
    integer(IK)                       :: i, n, m
!
    if (PRESENT(sym)) then
      n = SIZE(sym, 1)
      m = SIZE(sym, 2)
      ALLOCATE(res%p(m))
      do concurrent(i=1:m)
        res%p(i) = group_permutation(sym(:,i))
      enddo
    end if
!
  end function molecular_rotation_new
!
  pure elemental function molecular_rotation_nsym(this) result(res)
    class(molecular_rotation), intent(in) :: this
    integer(IK)                           :: res
    if (ALLOCATED(this%p)) then
      res = SIZE(this%p)
    else
      res = 0
    end if
  end function molecular_rotation_nsym
!
  pure subroutine molecular_rotation_swap(this, d, X, isym)
    class(molecular_rotation), intent(in) :: this
    integer(IK), intent(in)               :: d, isym
    real(RK), intent(inout)               :: X(*)
    if (isym < 1 .or. this%nsym() < isym) return
    call this%p(isym)%swap(d, X)
  end subroutine molecular_rotation_swap
!
  pure elemental subroutine molecular_rotation_clear(this)
    class(molecular_rotation), intent(inout) :: this
    if (ALLOCATED(this%p)) deallocate (this%p)
  end subroutine molecular_rotation_clear
!
  pure elemental subroutine molecular_rotation_destroy(this)
    type(molecular_rotation), intent(inout) :: this
    call this%clear()
  end subroutine molecular_rotation_destroy
!
end module mod_molecular_rotation

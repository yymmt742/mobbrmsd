module mod_molecular_rotation
  use mod_params, only: IK, RK
  use mod_optarg
  use mod_group_permutation
  implicit none
  private
  public :: molecular_rotation
!
  integer(IK), parameter :: DEF_d = 3
  !! default spatial dimension.
  integer(IK), parameter :: DEF_m = 1
  !! default number of atom in molecule.
!
  type molecular_rotation
    private
    type(group_permutation), allocatable :: p(:)
    !! number of operation, must be f <= m.
    integer(IK), allocatable             :: f(:)
    !! free indices
  contains
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
    integer(IK)                       :: i, n, m
!
    if (PRESENT(sym)) then
      n = SIZE(sym, 1)
      m = SIZE(sym, 2)
      allocate (res%p(m))
      do concurrent(i=1:m)
        res%p(i) = group_permutation(sym(:, i))
      end do
      res%f = free_indices(n, m, res%p)
    else
      allocate (res%p(0))
      allocate (res%f(0))
    end if
!
  end function molecular_rotation_new
!
  pure function free_indices(n, m, p) result(res)
    integer(IK), intent(in)             :: n, m
    type(group_permutation), intent(in) :: p(*)
    logical                             :: w(n)
    integer(IK), allocatable            :: res(:)
    integer(IK)                         :: i
    do concurrent(i = 1:m)
      w(i) = .FALSE.
    enddo
    do i = 1, m
      w(p(i)%free_indices()) = .TRUE.
    end do
    res = PACK([(i, i=1, n)], w)
  end function free_indices
!
  pure elemental subroutine molecular_rotation_clear(this)
    class(molecular_rotation), intent(inout) :: this
    if (ALLOCATED(this%p)) deallocate (this%p)
    if (ALLOCATED(this%f)) deallocate (this%f)
  end subroutine molecular_rotation_clear
!
  pure elemental subroutine molecular_rotation_destroy(this)
    type(molecular_rotation), intent(inout) :: this
    call this%clear()
  end subroutine molecular_rotation_destroy
!
end module mod_molecular_rotation

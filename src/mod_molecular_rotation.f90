module mod_molecular_rotation
  use mod_params, only: IK, RK
  use mod_optarg
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
    integer(IK)              :: m = 0
    !! number of atom in molecule.
    integer(IK)              :: f = 0
    !! number of operation, must be f <= m.
    integer(IK), allocatable :: q(:, :)
    !! permuation matrix, q(2, f).
    !! dimension 1 is pointer to source index, 0 < q1 <= m.
    !! dimension 2 is pointer to destination index, 0 < q2 <= m.
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
  pure function molecular_rotation_new(m, sym) result(res)
    integer(IK), intent(in), optional :: m
    !! default number of molecule = 1.
    integer(IK), intent(in), optional :: sym(:, :)
    !! swap indices
    type(molecular_rotation)          :: res
!
    res%m = MAX(optarg(m, DEF_m), 1)
!
    if (PRESENT(sym)) then
      res%fix = fix
      res%free = free_indices(res%n, fix)
    end if
!
  end function molecular_rotation_new
!
  pure elemental subroutine molecular_rotation_clear(this)
    class(molecular_rotation), intent(inout) :: this
    this%m = 0
    this%f = 0
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine molecular_rotation_clear
!
  pure elemental subroutine molecular_rotation_destroy(this)
    type(molecular_rotation), intent(inout) :: this
    call this%clear()
  end subroutine molecular_rotation_destroy
!
end module mod_molecular_rotation

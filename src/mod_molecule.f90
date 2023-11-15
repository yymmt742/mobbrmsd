module mod_molecule
  use mod_params, only: IK, RK
  use mod_element
  implicit none
  private
  public :: molecule
!
  type :: molecule
    private
    type(element), allocatable :: e(:)
    integer(IK), allocatable   :: s(:, :)
  contains
    procedure :: natom      => molecule_natom
    procedure :: nsym       => molecule_nsym
    procedure :: sym_index  => molecule_sym_index
    procedure :: free_index => molecule_free_index
    final     :: molecule_destroy
  end type molecule
!
  interface molecule
    module procedure molecule_new
  end interface molecule
!
contains
!
  pure function molecule_new(n, sym) result(res)
    integer(IK), intent(in)           :: n
    integer(IK), intent(in), optional :: sym(:, :)
    type(molecule)                    :: res
    integer(IK)                       :: l
!
    l = MAX(n, 0)
    allocate (res%e(l))
!
    if (PRESENT(sym)) then
      block
        integer(IK) :: r, s
        r = MIN(l, SIZE(sym, 1))
        s = SIZE(sym, 2)
        allocate (res%s(l, s), source=0)
        res%s(:r, :s) = sym(:r, :s)
      end block
    else
      allocate (res%s(l, 0))
    end if
!
  end function molecule_new
!
  pure elemental function molecule_natom(this) result(res)
    class(molecule), intent(in) :: this
    integer(IK)                 :: res
    if (ALLOCATED(this%e)) then
      res = SIZE(this%e)
    else
      res = 0
    end if
  end function molecule_natom
!
  pure elemental function molecule_nsym(this) result(res)
    class(molecule), intent(in) :: this
    integer(IK)                 :: res
    if (ALLOCATED(this%s)) then
      res = SIZE(this%s, 2) + 1
    else
      res = 0
    end if
  end function molecule_nsym
!
  pure function molecule_sym_index(this, s) result(res)
    class(molecule), intent(in) :: this
    integer(IK), intent(in)     :: s
    integer(IK)                 :: res(this%natom())
    integer(IK)                 :: i
    if (1 < s .and. s <= this%nsym()) then
      do concurrent(i=1:this%natom())
        res(i) = this%s(i, s - 1)
      end do
    else
      do concurrent(i=1:this%natom())
        res(i) = i
      end do
    end if
  end function molecule_sym_index
!
  pure function molecule_free_index(this) result(res)
    class(molecule), intent(in) :: this
    integer(IK), allocatable    :: res(:)
    integer(IK)                 :: i
    allocate (res(0))
    if (ALLOCATED(this%s)) then
      do i = 1, this%natom()
        if (ALL(this%s(:, i)==i)) cycle
        res = [res, i]
      end do
    end if
  end function molecule_free_index
!
  pure elemental subroutine molecule_destroy(this)
    type(molecule), intent(inout) :: this
!
    if (ALLOCATED(this%e)) deallocate (this%e)
    if (ALLOCATED(this%s)) deallocate (this%s)
!
  end subroutine molecule_destroy
!
end module mod_molecule

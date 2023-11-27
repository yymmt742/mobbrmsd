module mod_molecule
  use mod_params, only: IK, RK
  use mod_element
  implicit none
  private
  public :: molecule, molecules
!
  type :: molecule
    private
    integer(IK)                :: f = 0
    type(element), allocatable :: e(:)
    integer(IK), allocatable   :: s(:, :)
  contains
    procedure :: natom      => molecule_natom
    procedure :: nsym       => molecule_nsym
    procedure :: nrot       => molecule_nrot
    procedure :: sym_index  => molecule_sym_index
    procedure :: free_index => molecule_free_index
    procedure :: clear      => molecule_clear
    final     :: molecule_destroy
  end type molecule
!
  type :: molecules
    private
    type(molecule), allocatable :: m(:)
    integer(IK), allocatable    :: l(:)
  contains
    procedure :: nmol       => molecules_nmol
    procedure :: nspecies   => molecules_nspecies
    procedure :: clear      => molecules_clear
    final     :: molecules_destroy
  end type molecules
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
        integer(IK) :: i, r, s
        r = MIN(l, SIZE(sym, 1))
        s = SIZE(sym, 2)
        allocate (res%s(l, s), source=0)
        res%s(:r, :s) = sym(:r, :s)
        do i = 1, r
          if (ALL(res%s(i, :) == i)) cycle
          res%f = res%f + 1
        end do
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
  pure elemental function molecule_nrot(this) result(res)
    class(molecule), intent(in) :: this
    integer(IK)                 :: res
    res = this%f
  end function molecule_nrot
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
    integer(IK)                 :: res(this%f)
    integer(IK)                 :: i, j
    if (ALLOCATED(this%s)) then
      j = 1
      do i = 1, this%natom()
        if (ALL(this%s(i, :)==i)) cycle
        res(j) = i
        j = j + 1
      end do
    end if
  end function molecule_free_index
!
  pure elemental subroutine molecule_clear(this)
    class(molecule), intent(inout) :: this
    if (ALLOCATED(this%e)) deallocate (this%e)
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine molecule_clear
!
  pure elemental subroutine molecule_destroy(this)
    type(molecule), intent(inout) :: this
    call this%clear()
  end subroutine molecule_destroy
!
  pure elemental function molecules_nmol(this) result(res)
    class(molecules), intent(in) :: this
    integer(IK)                 :: res
    if (ALLOCATED(this%l)) then
      res = SIZE(this%l)
    else
      res = 0
    end if
  end function molecules_nmol
!
  pure elemental function molecules_nspecies(this) result(res)
    class(molecules), intent(in) :: this
    integer(IK)                 :: res
    if (ALLOCATED(this%m)) then
      res = SIZE(this%m)
    else
      res = 0
    end if
  end function molecules_nspecies
!
  pure elemental subroutine molecules_clear(this)
    class(molecules), intent(inout) :: this
    if (ALLOCATED(this%m)) deallocate (this%m)
    if (ALLOCATED(this%l)) deallocate (this%l)
  end subroutine molecules_clear
!
  pure elemental subroutine molecules_destroy(this)
    type(molecules), intent(inout) :: this
    call this%clear()
  end subroutine molecules_destroy
!
end module mod_molecule

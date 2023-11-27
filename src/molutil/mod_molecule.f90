module mod_molecule
  use mod_params, only: IK, RK
  use mod_optarg
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
    procedure :: unfmtread  => molecule_unfmtread
    generic   :: read (unformatted) => unfmtread
    procedure :: unfmtwrite => molecule_unfmtwrite
    generic   :: write (unformatted) => unfmtwrite
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
  interface molecules
    module procedure molecules_new
  end interface molecules
!
contains
!
  pure function molecule_new(n, d, sym, ele) result(res)
    integer(IK), intent(in)             :: n
    integer(IK), intent(in), optional   :: d
    integer(IK), intent(in), optional   :: sym(*)
    type(element), intent(in), optional :: ele(*)
    type(molecule)                      :: res
    integer(IK)                         :: i, j, n_, d_
!
    n_ = MAX(n, 0)
    d_ = optarg(d, 0)
!
    allocate (res%e(n_))
!
    if (PRESENT(ele)) then
      do concurrent(i=1:n_)
        res%e(i) = ele(i)
      end do
    end if
!
    if (PRESENT(sym)) then
      allocate (res%s(n_, d_))
      do concurrent(i=1:n_, j=1:d_)
         block
           integer(IK) :: k
           k = i + (j - 1) * n_
           res%s(i, j) = sym(k)
         end block
      end do
    else
      allocate (res%s(n_, d_), source=0_IK)
    end if
!
    res%f = 0
    do i = 1, d_
      if (ALL(res%s(i, :) == i)) cycle
      res%f = res%f + 1
    end do
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
  subroutine molecule_unfmtread(this, unit, iostat, iomsg)
    class(molecule), intent(inout) :: this
    integer(IK), intent(in)        :: unit
    integer(IK), intent(out)       :: iostat
    character(*), intent(inout)    :: iomsg
    integer(IK)                    :: i, n, d
!
    call molecule_clear(this)
!
    read (unit, IOSTAT=iostat, IOMSG=iomsg) n, d
!
    if (n < 0) return
    allocate (this%e(n))
    read (unit, IOSTAT=iostat, IOMSG=iomsg) this%e
!
    if (d < 1) return
    allocate (this%s(n, d))
    read (unit, IOSTAT=iostat, IOMSG=iomsg) this%s
    do i = 1, d
      if (ALL(this%s(i, :) == i)) cycle
      this%f = this%f + 1
    end do
!
  end subroutine molecule_unfmtread
!
  subroutine molecule_unfmtwrite(this, unit, iostat, iomsg)
    class(molecule), intent(in) :: this
    integer(IK), intent(in)     :: unit
    integer(IK), intent(out)    :: iostat
    character(*), intent(inout) :: iomsg
    integer(IK)                 :: n, d
!
    n = molecule_natom(this)
    d = molecule_nsym(this) - 1
!
    write (unit, IOSTAT=iostat, IOMSG=iomsg) n, d
    if(n<1) return
    write (unit, IOSTAT=iostat, IOMSG=iomsg) this%e
    if(d<1) return
    write (unit, IOSTAT=iostat, IOMSG=iomsg) this%s
!
  end subroutine molecule_unfmtwrite
!
  pure elemental subroutine molecule_clear(this)
    class(molecule), intent(inout) :: this
    this%f = 0
    if (ALLOCATED(this%e)) deallocate (this%e)
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine molecule_clear
!
  pure elemental subroutine molecule_destroy(this)
    type(molecule), intent(inout) :: this
    call this%clear()
  end subroutine molecule_destroy
!
  pure function molecules_new(mol, list) result(res)
    class(molecule), intent(in) :: mol(:)
    integer(IK), intent(in)     :: list(:)
    type(molecules)             :: res
    integer(IK)                 :: i, s, n
!
    s = SIZE(mol)
    n = SIZE(list)
    allocate (res%m(s))
    do concurrent(i=1:s)
      res%m(i) = mol(i)
    end do
!
    allocate (res%l, source=list)
!
  end function molecules_new
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

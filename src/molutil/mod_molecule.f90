module mod_molecule
  use mod_params, only: IK, RK
  use mod_optarg
  use mod_base
  use mod_bonds
  implicit none
  private
  public :: molecule
!
  integer(IK), parameter  :: MOLECULE_CHARLEN = 64
  character(*), parameter :: DEF_molname      = 'default_name'
!
  type, extends(base) :: molecule
    private
    character(MOLECULE_CHARLEN)  :: c = DEF_molname
    integer(IK)                  :: f = 0
    integer(IK), allocatable     :: s(:, :)
    type(bonds)                  :: b
  contains
    procedure :: natom               => molecule_natom
    procedure :: nsym                => molecule_nsym
    procedure :: nrot                => molecule_nrot
    procedure :: sym_indices         => molecule_sym_indices
    procedure :: free_indices        => molecule_free_indices
    procedure :: molname             => molecule_molname
    procedure :: formatted_string    => molecule_formatted_string
    procedure :: clear               => molecule_clear
    procedure :: unfmtread           => molecule_unfmtread
    generic   :: read (unformatted)  => unfmtread
    procedure :: unfmtwrite          => molecule_unfmtwrite
    generic   :: write (unformatted) => unfmtwrite
    final     :: molecule_destroy
  end type molecule
!
  interface molecule
    module procedure molecule_new
  end interface molecule
!
contains
!
  pure function molecule_new(n, d, sym, atoms, molname) result(res)
    integer(IK), intent(in)             :: n
    integer(IK), intent(in), optional   :: d
    integer(IK), intent(in), optional   :: sym(*)
    integer(IK), intent(in), optional   :: atoms(*)
    character(*), intent(in), optional  :: molname
    type(molecule)                      :: res
    integer(IK)                         :: i, j, n_, d_
!
    n_ = MAX(n, 0)
    d_ = optarg(d, 1)
!
    allocate (res%s(n_, d_), source=0)
!
    if (PRESENT(atoms)) then
      do concurrent(i=1:n_)
        res%s(i, 1) = atoms(i)
      end do
    end if
!
    if (PRESENT(sym)) then
      do concurrent(i=1:n_, j=2:d_)
         block
           integer(IK) :: k
           k = i + (j - 2) * n_
           res%s(i, j) = sym(k)
         end block
      end do
    end if
!
    if (PRESENT(molname)) then
      res%c = molname
    else
      res%c = DEF_molname
    end if
!
    if(d_<2) return
!
    res%f = 0
    do i = 1, n_
      if (ALL(res%s(i, 2:) == i)) cycle
      res%f = res%f + 1
    end do
!
  end function molecule_new
!
  pure elemental function molecule_natom(this) result(res)
    class(molecule), intent(in) :: this
    integer(IK)                 :: res
    if (ALLOCATED(this%s)) then
      res = SIZE(this%s, 1)
    else
      res = 0
    end if
  end function molecule_natom
!
  pure elemental function molecule_nsym(this) result(res)
    class(molecule), intent(in) :: this
    integer(IK)                 :: res
    if (ALLOCATED(this%s)) then
      res = SIZE(this%s, 2)
    else
      res = 0
    end if
  end function molecule_nsym
!
  pure function molecule_molname(this) result(res)
    class(molecule), intent(in) :: this
    character(:), allocatable   :: res
    res = TRIM(ADJUSTL(this%c))
  end function molecule_molname
!
  pure elemental function molecule_nrot(this) result(res)
    class(molecule), intent(in) :: this
    integer(IK)                 :: res
    res = this%f
  end function molecule_nrot
!
  pure function molecule_sym_indices(this, s) result(res)
    class(molecule), intent(in) :: this
    integer(IK), intent(in)     :: s
    integer(IK)                 :: res(this%natom())
    integer(IK)                 :: i, n
    n = this%natom()
    if (1 < s .and. s <= this%nsym()) then
      do concurrent(i=1:n)
        res(i) = this%s(i, s)
      end do
    else
      do concurrent(i=1:n)
        res(i) = i
      end do
    end if
  end function molecule_sym_indices
!
  pure function molecule_free_indices(this) result(res)
    class(molecule), intent(in) :: this
    integer(IK)                 :: res(this%f)
    integer(IK)                 :: i, j
    if (this%nsym() < 1) return
    if (ALLOCATED(this%s)) then
      j = 1
      do i = 1, this%natom()
        if (ALL(this%s(i, 2:) == i)) cycle
        res(j) = i
        j = j + 1
      end do
    end if
  end function molecule_free_indices
!
  pure function molecule_formatted_string(this) result(res)
  use mod_params,  only: NL => NEWLINE
    class(molecule), intent(in) :: this
    character(:), allocatable   :: res
    integer(IK)                 :: i, j, n, s, m, w
!
    n = this%natom()
    s = this%nsym()
    m = LEN_TRIM(this%c) + LEN(NL)
    w = 8 + 8 * s + LEN(NL)
    i = n * w + m
    allocate (character(i) :: res)
!
    res(:m) = TRIM(this%c)//NL
!
    if (n < 1) return
!
    do concurrent(i=1:n)
      write (res(1 + (i - 1) * w + m:8 + (i - 1) * w + m), '(I8)') i
      write (res(9 + (i - 1) * w + m:16 + (i - 1) * w + m), '(I8)') i
      do concurrent(j=2:s)
        block
          integer(IK) :: l, u
          l = 8 * (j - 1) + (i - 1) * w + m + 8 + 1
          u = l + 8 - 1
          write (res(l:u), '(I8)') this%s(i, j)
        end block
      end do
      res(i * w + m - LEN(NL) + 1:i * w + m) = NL
    end do
!
  end function molecule_formatted_string
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
    read (unit, IOSTAT=iostat, IOMSG=iomsg) n, d, this%c
!
    if (n < 0) return
    allocate (this%s(n, d))
    read (unit, IOSTAT=iostat, IOMSG=iomsg) this%s
!
    if (d < 2) return
    do i = 1, d
      if (ALL(this%s(i, 2:) == i)) cycle
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
    write (unit, IOSTAT=iostat, IOMSG=iomsg) n, d, this%c
    if(.not.ALLOCATED(this%s)) return
    write (unit, IOSTAT=iostat, IOMSG=iomsg) this%s
!
  end subroutine molecule_unfmtwrite
!
  pure elemental subroutine molecule_clear(this)
    class(molecule), intent(inout) :: this
    this%f = 0
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine molecule_clear
!
  pure elemental subroutine molecule_destroy(this)
    type(molecule), intent(inout) :: this
    call this%clear()
  end subroutine molecule_destroy
!
end module mod_molecule

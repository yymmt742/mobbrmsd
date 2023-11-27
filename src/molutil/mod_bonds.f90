module mod_bonds
  use mod_params, only: IK, RK
  use mod_base
  implicit none
  private
  public :: bonds
!
  type bond
    sequence
    integer(IK) :: a, b
    real(RK)    :: c, e
  end type bond
!
  type, extends(base) :: bonds
    private
    type(bond), allocatable :: b(:)
  contains
    procedure :: nbond               => bonds_nbond
    procedure :: clear               => bonds_clear
    procedure :: unfmtread           => bonds_unfmtread
    generic   :: read (unformatted)  => unfmtread
    procedure :: unfmtwrite          => bonds_unfmtwrite
    generic   :: write (unformatted) => unfmtwrite
    final     :: bonds_destroy
  end type bonds
!
contains
!
  pure elemental function bonds_nbond(this) result(res)
    class(bonds), intent(in) :: this
    integer(IK)              :: res
    if (ALLOCATED(this%b)) then
      res = SIZE(this%b)
    else
      res = 0
    end if
  end function bonds_nbond
!
  subroutine bonds_unfmtread(this, unit, iostat, iomsg)
    class(bonds), intent(inout) :: this
    integer(IK), intent(in)     :: unit
    integer(IK), intent(out)    :: iostat
    character(*), intent(inout) :: iomsg
    integer(IK)                 :: n
!
    call bonds_clear(this)
    read (unit, IOSTAT=iostat, IOMSG=iomsg) n
    if (n < 0) return
    allocate (this%b(n))
    read (unit, IOSTAT=iostat, IOMSG=iomsg) this%b
!
  end subroutine bonds_unfmtread
!
  subroutine bonds_unfmtwrite(this, unit, iostat, iomsg)
    class(bonds), intent(in)    :: this
    integer(IK), intent(in)     :: unit
    integer(IK), intent(out)    :: iostat
    character(*), intent(inout) :: iomsg
    integer(IK)                 :: n
!
    n = bonds_nbond(this)
    write (unit, IOSTAT=iostat, IOMSG=iomsg) n
    if(.not.ALLOCATED(this%b)) return
    write (unit, IOSTAT=iostat, IOMSG=iomsg) this%b
!
  end subroutine bonds_unfmtwrite
!
  pure elemental subroutine bonds_clear(this)
    class(bonds), intent(inout) :: this
    if (ALLOCATED(this%b)) deallocate (this%b)
  end subroutine bonds_clear
!
  pure elemental subroutine bonds_destroy(this)
    type(bonds), intent(inout) :: this
    call this%clear()
  end subroutine bonds_destroy
!
end module mod_bonds

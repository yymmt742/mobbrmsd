module mod_base
  implicit none
  private
  public :: base
!
  type :: base
  contains
    procedure :: formatted_string  => base_formatted_string
    procedure :: fmtwrite          => base_fmtwrite
    generic   :: write (formatted) => fmtwrite
  end type base
!
contains
!
  pure function base_formatted_string(this) result(res)
    class(base), intent(in)   :: this
    character(:), allocatable :: res
    allocate (character(0) :: res)
  end function base_formatted_string
!
  subroutine base_fmtwrite(this, unit, iotype, vlist, iostat, iomsg)
    class(base), intent(in)     :: this
    integer, intent(in)         :: unit
    character(*), intent(in)    :: iotype
    integer, intent(in)         :: vlist(:)
    integer, intent(out)        :: iostat
    character(*), intent(inout) :: iomsg
    write (unit, '(A)', IOSTAT=iostat, IOMSG=iomsg) this%formatted_string()
  end subroutine base_fmtwrite
!
end module mod_base

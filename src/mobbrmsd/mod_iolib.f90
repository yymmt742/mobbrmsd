module mod_iolib
  use ISO_C_BINDING
  use ISO_FORTRAN_ENV, only: INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
  implicit none
  private
  public :: isatty
!
  interface
    function isatty_stdin() bind(c)
      use ISO_C_BINDING
      integer(C_INT) :: isatty_stdin
    end function isatty_stdin
!
    function isatty_stdout() bind(c)
      use ISO_C_BINDING
      integer(C_INT) :: isatty_stdout
    end function isatty_stdout
!
    function isatty_stderr() bind(c)
      use ISO_C_BINDING
      integer(C_INT) :: isatty_stderr
    end function isatty_stderr
  end interface
!
contains
!
  function isatty(unit) result(res)
    integer, intent(in), optional :: unit
    logical                       :: res
    integer                       :: cres
    if (.not. PRESENT(unit)) then
      res = isatty_stdout() > 0
      return
    end if
    select case (unit)
    case (INPUT_UNIT)
      cres = isatty_stdin()
    case (OUTPUT_UNIT)
      cres = isatty_stdout()
    case (ERROR_UNIT)
      cres = isatty_stderr()
    case default
      cres = 0
    end select
    res = cres > 0
  end function isatty
!
end module mod_iolib

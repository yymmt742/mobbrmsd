module mod_optarg
  use ISO_FORTRAN_ENV, only: &
 &   INT8, &
 &   INT16, &
 &   INT32, &
 &   INT64, &
 &   REAL32, &
 &   REAL64, &
 &   REAL128
  implicit none
  private
  public optarg
!
  interface optarg
    module procedure :: &
      optarg_l, optarg_c, &
   &  optarg_i8, optarg_i16, optarg_i32, optarg_i64, &
   &  optarg_r32, optarg_r64, optarg_r128
  end interface optarg
!
contains
!
  pure elemental function optarg_l(arg, def) result(res)
    logical, intent(in), optional :: arg
    logical, intent(in)           :: def
    logical                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_l
!
  pure elemental function optarg_c(arg, def) result(res)
    character, intent(in), optional :: arg
    character, intent(in)           :: def
    character                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_c
!
  pure elemental function optarg_i8(arg, def) result(res)
    integer(INT8), intent(in), optional :: arg
    integer(INT8), intent(in)           :: def
    integer(INT8)                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_i8
!
  pure elemental function optarg_i16(arg, def) result(res)
    integer(INT16), intent(in), optional :: arg
    integer(INT16), intent(in)           :: def
    integer(INT16)                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_i16
!
  pure elemental function optarg_i32(arg, def) result(res)
    integer(INT32), intent(in), optional :: arg
    integer(INT32), intent(in)           :: def
    integer(INT32)                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_i32
!
  pure elemental function optarg_i64(arg, def) result(res)
    integer(INT64), intent(in), optional :: arg
    integer(INT64), intent(in)           :: def
    integer(INT64)                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_i64
!
  pure elemental function optarg_r32(arg, def) result(res)
    real(REAL32), intent(in), optional :: arg
    real(REAL32), intent(in)           :: def
    real(REAL32)                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_r32
!
  pure elemental function optarg_r64(arg, def) result(res)
    real(REAL64), intent(in), optional :: arg
    real(REAL64), intent(in)           :: def
    real(REAL64)                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_r64
!
  pure elemental function optarg_r128(arg, def) result(res)
    real(REAL128), intent(in), optional :: arg
    real(REAL128), intent(in)           :: def
    real(REAL128)                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_r128
!
end module mod_optarg


!| Optarg
module mod_optarg
  use mod_params, only: IK, RK, ZERO => RZERO
  implicit none
  private
  public :: optarg
!
  interface optarg
    module procedure :: optarg_l0, optarg_l0_, &
                      & optarg_i0, optarg_r0
  end interface optarg
!
contains
!
  pure elemental function optarg_l0(arg) result(res)
    logical, intent(in), optional :: arg
    logical                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = .false.
    end if
  end function optarg_l0
!
  pure elemental function optarg_l0_(arg, def) result(res)
    logical, intent(in), optional :: arg
    logical, intent(in)           :: def
    logical                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_l0_
!
  pure elemental function optarg_i0(arg, def) result(res)
    integer(IK), intent(in), optional :: arg
    integer(IK), intent(in)           :: def
    integer(IK)                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_i0
!
  pure elemental function optarg_r0(arg, def) result(res)
    real(RK), intent(in), optional :: arg
    real(RK), intent(in)           :: def
    real(RK)                       :: res
    if (PRESENT(arg)) then
      res = arg
    else
      res = def
    end if
  end function optarg_r0
!
end module mod_optarg

pure elemental function invcbrt(x) result(res)
  use mod_kinds, only: I4, R4, I8
  implicit none
  real(RK), intent(in) :: x
  real(RK)             :: res, c
  real(R4)             :: y
  y = ABS(x)
  res = SIGN(1.0_RK, x) * TRANSFER(INT(z"548C2B4B", I4) - TRANSFER(y, 0_I4) / 3, y)
  c = res * res * res * x
  res = res * (1.752319676_RK - c * (1.2509524245_RK - 0.5093818292_RK * c))
  c = 1.0_RK - res * res * res * x
  res = res * (1.0_RK + 0.333333333333333_RK * c)
! res = res * ((4.0_RK / 3.0_RK) - res * res * res * fx)
! res = res * ((4.0_RK / 3.0_RK) - res * res * res * fx)
end function invcbrt

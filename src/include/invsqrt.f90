pure elemental function invsqrt(x) result(res)
  use mod_kinds, only: I4, R4, I8
  implicit none
  real(RK), intent(in) :: x
  real(RK)             :: res
  real(R4)             :: y
  y = x
  res = TRANSFER(INT(z"5F1FFFF9", I4) - ISHFT(TRANSFER(y, 0_I4), -1), y)
  !res = TRANSFER(INT(z"5F3759DF", I4) - ISHFT(TRANSFER(y, 0_I4), -1), y)
  res = 0.703952253_RK * res * (2.38924456_RK - x * res * res)
  res = res * (1.5_RK - 0.5_RK * x * res * res)
  res = res * (1.5_RK - 0.5_RK * x * res * res)
end function invsqrt

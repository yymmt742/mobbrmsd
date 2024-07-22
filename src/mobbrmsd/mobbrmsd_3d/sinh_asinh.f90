  pure elemental subroutine sinh_asinh(x, res, yy)
    real(RK), intent(in)    :: x
    real(RK), intent(inout) :: res, yy
    yy = x * x
    if (yy < 1.69_RK) then ! 1.3^2
      res = x * (1.0_RK / 3.0_RK) - x * yy * (4.0_RK / 81.0_RK)
    else
      ! yy = x**(-1/3)
      call invcbrt(x, yy, res)
      ! res = 2**(-2/3) * x**(1/3) - 2**(-4/3) * x**(-1/3)
      res = 0.6299605249474366_RK / yy - 0.3968502629920499_RK * yy
    end if
    yy = res * res
    res = res - ((4.0_RK * yy + 3.0_RK) * res - x) / (12.0_RK * yy + 3.0_RK)
    yy = res * res
    res = res - ((4.0_RK * yy + 3.0_RK) * res - x) / (12.0_RK * yy + 3.0_RK)
  end subroutine sinh_asinh


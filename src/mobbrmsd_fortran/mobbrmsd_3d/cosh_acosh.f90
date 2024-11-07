  pure elemental subroutine cosh_acosh(x, res, yy, df)
    real(RK), intent(in)    :: x
    real(RK), intent(inout) :: res, yy, df
    integer(IK)             :: i
    if (x < ZERO) return
    if (x > ONE) return
    if (x < 5.0000000000000000E-01_RK) then
      if (x < 1.4644660940672627E-01_RK) then
        res = (2.9546700898745084E+00_RK) * x + (-6.8732996465706475E-01_RK)
        res = res * x + (-1.3519484012953734E-02_RK)
        res = res * x + (-4.6741735498024902E-03_RK)
        res = res * x + (0.0000000000000000E+00_RK)
      else
        res = (-1.0623675833993933E-01_RK) * x + (1.5654255346132895E-01_RK)
        res = res * x + (-9.6794886272540959E-02_RK)
        res = res * x + (1.8177193622178180E-03_RK)
        res = res * x + (-4.0727563438808847E-04_RK)
      end if
    else
      if (x < 8.5355339059327373E-01_RK) then
        res = (2.3919197448076808E-02_RK) * x + (-5.9903084576736369E-02_RK)
        res = res * x + (5.0794386881566130E-02_RK)
        res = res * x + (-4.7880140922667278E-02_RK)
        res = res * x + (6.4652937375348487E-03_RK)
      else
        res = (-1.3435433998222815E-01_RK) * x + (5.0040317286635916E-01_RK)
        res = res * x + (-6.9989274807980062E-01_RK)
        res = res * x + (4.0248132481808807E-01_RK)
        res = res * x + (-9.5448197561904868E-02_RK)
      end if
    end if
    ! res = 2**(-2/3) * x**(1/3) - 2**(-4/3) * x**(-1/3)
    call invcbrt(x, yy, df)
    res = res + 0.6299605249474366_RK / yy + 0.3968502629920499_RK * yy
    do i = 1, 5
      yy = res * res
      df = 12.0_RK * yy - 3.0_RK
      if (ABS(df) < 1E-18_RK) return
      df = ((4.0_RK * yy - 3.0_RK) * res - ONE / x) / df
      res = res - df
      if (ABS(df) < 1E-18_RK) return
    end do
  end subroutine cosh_acosh


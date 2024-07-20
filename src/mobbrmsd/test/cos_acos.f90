  pure elemental function cos_acos(x) result(res)
    real(RK), intent(in) :: x
    real(RK)             :: res, yy, df
    if (x < -1.0000000000000000E+00_RK) then
      res = 1.0000000000000000E+00_RK / 2.0000000000000000E+00_RK
      return
    elseif (x < -0.3333333333333333E+00_RK) then
      block
        real(RK) :: a, b, c, d, e
        e = sqrt(x + ONE)
        if (x < -6.6666666666666663e-01_RK) then
          if (x < -9.0236892706218241e-01_RK) then
            a = 3.9847929621330180e-04 * e + (2.7150043732313152e-03_RK)
            b = 2.7756772944468894e-03 * e + (8.0584699091552827e-03_RK)
            c = 9.9787382004175255e-03 * e + (-3.0422013335204556e-04_RK)
            d = 3.2124745520119027e-02 * e + (-6.9480186781926440e-02_RK)
            e = 4.3277149578179819e-01 * e + (4.3616749888734951e-01_RK)
          else
            a = 3.9847929621330180e-04 * e + (-8.1880793225013100e-04_RK)
            b = 2.7756772944468894e-03 * e + (-5.3682021468007667e-03_RK)
            c = 9.9787382004175255e-03 * e + (-1.9427315377984092e-02_RK)
            d = 3.2124745520119027e-02 * e + (-8.1580601406840383e-02_RK)
            e = 4.3277149578179819e-01 * e + (4.3329732093336737e-01_RK)
          endif
        else
          if (x < -4.3096440627115074e-01_RK) then
            a = 3.9847929621330180e-04 * e + (-1.3962794333934664e-03_RK)
            b = 2.7756772944468894e-03 * e + (-6.7281227183138212e-03_RK)
            c = 9.9787382004175255e-03 * e + (-2.0640121298405912e-02_RK)
            d = 3.2124745520119027e-02 * e + (-8.2067600037457583e-02_RK)
            e = 4.3277149578179819e-01 * e + (4.3322280904921723e-01_RK)
          else
            a = 3.9847929621330180e-04 * e + (1.1713621213896104e-02_RK)
            b = 2.7756772944468894e-03 * e + (1.3560409301341475e-02_RK)
            c = 9.9787382004175255e-03 * e + (-8.8387594483930760e-03_RK)
            d = 3.2124745520119027e-02 * e + (-7.9004467646912643e-02_RK)
            e = 4.3277149578179819e-01 * e + (4.3352276164963782e-01_RK)
          endif
        endif
        call p4(a, b, c, d, e, x, res)
      end block
    elseif (x < 0.3333333333333333E+00) then
      if (x < 2.0410779985789219e-17_RK) then
        if (x < -2.3570226039551581e-01_RK) then
          call p4(6.3261501181342716e-01_RK, &
            &     7.7758615573901924e-01_RK, &
            &     2.7435256739631902e-01_RK, &
            &     2.2752170341132330e-01_RK, &
            &     8.7030281603695292e-01_RK, &
            &     x, res)
        else
          call p4(-5.1407513495624758e-02_RK, &
            &     1.0016161738288792e-02_RK, &
            &     -5.0149303892811525e-02_RK, &
            &     1.6657415217801974e-01_RK, &
            &     8.6602540378443871e-01_RK, &
            &     x, res)
        endif
      else
        if (x < 2.3570226039551587e-01_RK) then
          call p4(9.3008662532715076e-03_RK, &
            &     1.4372359674536180e-02_RK, &
            &     -4.6656484877630987e-02_RK, &
            &     1.6659959026260604e-01_RK, &
            &     8.6602540378443871e-01_RK, &
            &     x, res)
        else
          call p4(-3.4902175486090642e-01_RK, &
            &     4.1033381483828746e-01_RK, &
            &     -2.1260486589396638e-01_RK, &
            &     1.9761921280113542e-01_RK, &
            &     8.6385435215153961e-01_RK, &
            &     x, res)
        endif
      endif
    elseif (x <= 1.0000000000000000E+00_RK) then
      if (x < 6.6666666666666663e-01_RK) then
        if (x < 4.3096440627115085e-01_RK) then
          call p4(9.2989836420019442e-02_RK, &
            &     -1.3124380721245635e-01_RK, &
            &     3.9551533881730334e-02_RK, &
            &     1.4461452858482687e-01_RK, &
            &     8.6810669969142074e-01_RK, &
            &     x, res)
        else
          call p4(-6.7055964504391377e-03_RK, &
            &     2.2911999129408750e-02_RK, &
            &     -5.0039084385563551e-02_RK, &
            &     1.6784857275128981e-01_RK, &
            &     8.6583329929197761e-01_RK, &
            &     x, res)
        endif
      else
        if (x < 9.0236892706218252e-01_RK) then
          call p4(-1.0881827591406440e-03_RK, &
            &     9.1295179354582787e-03_RK, &
            &     -3.7133080284270058e-02_RK, &
            &     1.6236757221824985e-01_RK, &
            &     8.6672538337508409e-01_RK, &
            &     x, res)
        else
          call p4(4.6017500027684044e-03_RK, &
            &     -1.2962733889034270e-02_RK, &
            &     -5.1111183599694618e-03_RK, &
            &     1.4181596019335754e-01_RK, &
            &     8.7165614205287767e-01_RK, &
            &     x, res)
        endif
      endif
    else
      res = cosh_acosh(ONE / x)
      return
    end if
    yy = res * res
    df = 12.0_RK * yy - 3.0_RK
    if (ABS(df) < 1E-18_RK) return
    df = ((4.0_RK * yy - 3.0_RK) * res - x) / df
    res = res - df
  end function cos_acos
!
  pure elemental function cosh_acosh(x) result(res)
    real(RK), intent(in) :: x
    real(RK)             :: res, yy, df
    if (x < ZERO) then
      res = ONE
      return
    elseif (x > ONE) then
      res = cos_acos(ONE / x)
      return
    else
      if (x < 5.0000000000000000e-01_RK) then
        if (x < 1.4644660940672627e-01_RK) then
          call p4(2.9546700898745084e+00_RK, &
            &     -6.8732996465706475e-01_RK, &
            &     -1.3519484012953734e-02_RK, &
            &     -4.6741735498024902e-03_RK, &
            &     0.0000000000000000e+00_RK, &
            &     x, res)
        else
          call p4(-1.0623675833993933e-01_RK, &
            &     1.5654255346132895e-01_RK, &
            &     -9.6794886272540959e-02_RK, &
            &     1.8177193622178180e-03_RK, &
            &     -4.0727563438808847e-04_RK, &
            &     x, res)
        endif
      else
        if (x < 8.5355339059327373e-01_RK) then
          call p4(2.3919197448076808e-02_RK, &
            &     -5.9903084576736369e-02_RK, &
            &     5.0794386881566130e-02_RK, &
            &     -4.7880140922667278e-02_RK, &
            &     6.4652937375348487e-03_RK, &
            &     x, res)
        else
          call p4(-1.3435433998222815e-01_RK, &
            &     5.0040317286635916e-01_RK, &
            &     -6.9989274807980062e-01_RK, &
            &     4.0248132481808807e-01_RK, &
            &     -9.5448197561904868e-02_RK, &
            &     x, res)
        endif
      endif
    end if
    yy = 4.0_RK * x
    res = res + invcbrt(yy) + x * invcbrt(yy * yy)
    yy = res * res
    df = 12.0_RK * yy - 3.0_RK
    if (ABS(df) < 1E-18_RK) return
    df = ((4.0_RK * yy - 3.0_RK) * res - ONE / x) / df
    res = res - df
  end function cosh_acosh
!
  pure elemental function sinh_asinh(x) result(res)
    real(RK), intent(in) :: x
    real(RK)             :: res, yy
    yy = x * x
    if (yy < 1.69_RK) then ! 1.3^2
      res = x * (1.0_RK / 3.0_RK) - x * yy * (4.0_RK / 81.0_RK)
    else
      yy = invcbrt(2.0_RK * x)
      res = 0.5_RK * (1.0_RK / yy - yy)
    end if
    yy = res * res
    res = res - ((4.0_RK * yy + 3.0_RK) * res - x) / (12.0_RK * yy + 3.0_RK)
    yy = res * res
    res = res - ((4.0_RK * yy + 3.0_RK) * res - x) / (12.0_RK * yy + 3.0_RK)
  end function sinh_asinh
!
  pure elemental subroutine p4(a, b, c, d, e, x, res)
    real(RK), intent(in)    :: a, b, c, d, e, x
    real(RK), intent(inout) :: res
    res = a * x + b
    res = res * x + c
    res = res * x + d
    res = res * x + e
  end subroutine p4
!
#include "invcbrt.f90"

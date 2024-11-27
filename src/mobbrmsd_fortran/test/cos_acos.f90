  pure elemental function cos_acos(x) result(res)
    real(RK), intent(in) :: x
    real(RK)             :: res, yy, df
    if (x < -1.0000000000000000E+00_RK) then
      res = 1.0000000000000000E+00_RK / 2.0000000000000000E+00_RK
      return
    elseif (x < -0.3333333333333333E+00_RK) then
      yy = SQRT(x + ONE)
      if (x < -6.6666666666666663E-01_RK) then
        if (x < -9.0236892706218241E-01_RK) then
          res = (3.9847929621330180E-04 * yy + (2.7150043732313152E-03_RK)) * x + (2.7756772944468894E-03 * yy +&
              & (8.0584699091552827E-03_RK))
          res = res * x + (9.9787382004175255E-03 * yy + (-3.0422013335204556E-04_RK))
          res = res * x + (3.2124745520119027E-02 * yy + (-6.9480186781926440E-02_RK))
          res = res * x + (4.3277149578179819E-01 * yy + (4.3616749888734951E-01_RK))
        else
          res = (3.9847929621330180E-04 * yy + (-8.1880793225013100E-04_RK)) * x + (2.7756772944468894E-03 * yy +&
              & (-5.3682021468007667E-03_RK))
          res = res * x + (9.9787382004175255E-03 * yy + (-1.9427315377984092E-02_RK))
          res = res * x + (3.2124745520119027E-02 * yy + (-8.1580601406840383E-02_RK))
          res = res * x + (4.3277149578179819E-01 * yy + (4.3329732093336737E-01_RK))
        end if
      else
        if (x < -4.3096440627115074E-01_RK) then
          res = (3.9847929621330180E-04 * yy + (-1.3962794333934664E-03_RK)) * x + (2.7756772944468894E-03 * yy +&
              & (-6.7281227183138212E-03_RK))
          res = res * x + (9.9787382004175255E-03 * yy + (-2.0640121298405912E-02_RK))
          res = res * x + (3.2124745520119027E-02 * yy + (-8.2067600037457583E-02_RK))
          res = res * x + (4.3277149578179819E-01 * yy + (4.3322280904921723E-01_RK))
        else
          res = (3.9847929621330180E-04 * yy + (1.1713621213896104E-02_RK)) * x + (2.7756772944468894E-03 * yy +&
              & (1.3560409301341475E-02_RK))
          res = res * x + (9.9787382004175255E-03 * yy + (-8.8387594483930760E-03_RK))
          res = res * x + (3.2124745520119027E-02 * yy + (-7.9004467646912643E-02_RK))
          res = res * x + (4.3277149578179819E-01 * yy + (4.3352276164963782E-01_RK))
        end if
      end if
    elseif (x < 0.3333333333333333E+00) then
      if (x < 2.0410779985789219E-17_RK) then
        if (x < -2.3570226039551581E-01_RK) then
          res = (6.3261501181342716E-01_RK) * x + (7.7758615573901924E-01_RK)
          res = res * x + (2.7435256739631902E-01_RK)
          res = res * x + (2.2752170341132330E-01_RK)
          res = res * x + (8.7030281603695292E-01_RK)
        else
          res = (-5.1407513495624758E-02_RK) * x + (1.0016161738288792E-02_RK)
          res = res * x + (-5.0149303892811525E-02_RK)
          res = res * x + (1.6657415217801974E-01_RK)
          res = res * x + (8.6602540378443871E-01_RK)
        end if
      else
        if (x < 2.3570226039551587E-01_RK) then
          res = (9.3008662532715076E-03_RK) * x + (1.4372359674536180E-02_RK)
          res = res * x + (-4.6656484877630987E-02_RK)
          res = res * x + (1.6659959026260604E-01_RK)
          res = res * x + (8.6602540378443871E-01_RK)
        else
          res = (-3.4902175486090642E-01_RK) * x + (4.1033381483828746E-01_RK)
          res = res * x + (-2.1260486589396638E-01_RK)
          res = res * x + (1.9761921280113542E-01_RK)
          res = res * x + (8.6385435215153961E-01_RK)
        end if
      end if
    elseif (x <= 1.0000000000000000E+00_RK) then
      if (x < 6.6666666666666663E-01_RK) then
        if (x < 4.3096440627115085E-01_RK) then
          res = (9.2989836420019442E-02_RK) * x + (-1.3124380721245635E-01_RK)
          res = res * x + (3.9551533881730334E-02_RK)
          res = res * x + (1.4461452858482687E-01_RK)
          res = res * x + (8.6810669969142074E-01_RK)
        else
          res = (-6.7055964504391377E-03_RK) * x + (2.2911999129408750E-02_RK)
          res = res * x + (-5.0039084385563551E-02_RK)
          res = res * x + (1.6784857275128981E-01_RK)
          res = res * x + (8.6583329929197761E-01_RK)
        end if
      else
        if (x < 9.0236892706218252E-01_RK) then
          res = (-1.0881827591406440E-03_RK) * x + (9.1295179354582787E-03_RK)
          res = res * x + (-3.7133080284270058E-02_RK)
          res = res * x + (1.6236757221824985E-01_RK)
          res = res * x + (8.6672538337508409E-01_RK)
        else
          res = (4.6017500027684044E-03_RK) * x + (-1.2962733889034270E-02_RK)
          res = res * x + (-5.1111183599694618E-03_RK)
          res = res * x + (1.4181596019335754E-01_RK)
          res = res * x + (8.7165614205287767E-01_RK)
        end if
      end if
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
    end if
    call invcbrt(0.5_RK * x, yy, df)
    res = res + 0.5_RK * (1.0_RK / yy + yy)
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
      call invcbrt(2.0_RK * x, yy, res)
      res = 0.5_RK * (1.0_RK / yy - yy)
    end if
    yy = res * res
    res = res - ((4.0_RK * yy + 3.0_RK) * res - x) / (12.0_RK * yy + 3.0_RK)
    yy = res * res
    res = res - ((4.0_RK * yy + 3.0_RK) * res - x) / (12.0_RK * yy + 3.0_RK)
  end function sinh_asinh
!
!| https://www.mdpi.com/1996-1073/14/4/1058
  pure elemental subroutine invcbrt(x, res, c)
    use mod_kinds, only: I4, R4
    implicit none
    real(RK), intent(in)    :: x
    real(RK), intent(inout) :: res, c
    real(R4)                :: y
    y = ABS(x)
    res = SIGN(ONE, x) * TRANSFER(INT(z"548C2B4B", I4) - TRANSFER(y, 0_I4) / 3, y)
    c = res * res * res * x
    res = res * (1.752319676_RK - c * (1.2509524245_RK - 0.5093818292_RK * c))
    c = ONE - res * res * res * x
    res = res * (ONE + ONETHIRD * c)
  end subroutine invcbrt


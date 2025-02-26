module mod_quatic
  use mod_kinds
  use mod_params, only: RK, IK
  implicit none
  private
  public :: cos_acos
  public :: cosh_acosh
  public :: sinh_asinh
  public :: solve_quartic
  public :: solve_quartic_tuned
  public :: solve_quartic_newton
!
  real(RK), parameter    :: ZERO = 0.0_RK
  real(RK), parameter    :: ONE = 1.0_RK
  real(RK), parameter    :: TWO = 2.0_RK
  real(RK), parameter    :: THREE = 3.0_RK
  real(RK), parameter    :: FOUR = 4.0_RK
  real(RK), parameter    :: SIX = 6.0_RK
  real(RK), parameter    :: EIGHT = 8.0_RK
  real(RK), parameter    :: NINE = 9.0_RK
  real(RK), parameter    :: HALF = 0.5_RK
  real(RK), parameter    :: ONETHIRD = 1.0_RK / 3.0_RK
  real(RK), parameter    :: ONESIX = 1.0_RK / 6.0_RK
  real(RK), parameter    :: ONEQUARTER = 0.25_RK
  real(RK), parameter    :: ONEEIGHT = 0.125_RK
  real(RK), parameter    :: SQRT3 = SQRT(3.0_RK)
#ifdef USE_REAL32
  real(RK), parameter    :: THRESHOLD = 1E-6_RK
  real(RK), parameter    :: DEGENERACY = 1E-4_RK
#else
  real(RK), parameter    :: THRESHOLD = 1E-12_RK
  real(RK), parameter    :: DEGENERACY = 1E-6_RK
#endif
  integer(IK), parameter :: MAXITER = 1000
!
contains
  pure elemental subroutine solve_quartic_tuned(k1, k0, X)
    real(RK), intent(in)        :: k1, k0
    real(RK), intent(inout)     :: x
    real(RK)                    :: D, k11, k22
    real(RK)                    :: c1, c0, minima
    if (ABS(k1) < 1.E-1) then
      ! quasi biquadratic equations.
      k11 = k1 * k1
      k22 = 4._RK
      D = k0 - ONE
      minima = k1 + k0 - ONE
      if (minima >= ZERO) then
        ! multiple root
        x = ONE
      elseif (minima > -1.E-4) then
        c0 = k0 * ONEQUARTER - ONEQUARTER
        if (ABS(c0) < 1.E-8) then
          ! Second order Taylor expansion around the maximum local minima (x=1).
          ! f2(x) = x**2 - (2 - 1/4*k1) * x + (k0 + 3) / 4
          x = ONE - (k1 - SQRT(k11 - minima * 16._RK)) * ONEEIGHT
        else
          ! Third order Taylor expansion around the maximum local minima (x=1).
          c1 = ONE - ONEQUARTER * k1
          call find_a_cubic_root(c1, c0, x)
        end if
      else
        x = SQRT(MAX(ONE + SQRT(MAX(ONE - k0, 0._RK)), 0._RK))
        call newton(k1, k0, X, 5)
      end if
    else
      call exact(k1, k0, X)
    end if
  end subroutine solve_quartic_tuned
!
  pure elemental subroutine solve_quartic_newton(k1, k0, X)
    real(RK), intent(in)        :: k1, k0
    real(RK), intent(inout)     :: x
    real(RK)                    :: D, k11, k22
    real(RK)                    :: c1, c0, minima
    if (ABS(k1) < 1.E-1) then
      ! quasi biquadratic equations.
      k11 = k1 * k1
      k22 = 4._RK
      D = k0 - ONE
      minima = k1 + k0 - ONE
      if (minima >= ZERO) then
        ! multiple root
        x = ONE
      elseif (minima > -1.E-4) then
        c0 = k0 * ONEQUARTER - ONEQUARTER
        if (ABS(c0) < 1.E-8) then
          ! Second order Taylor expansion around the maximum local minima (x=1).
          ! f2(x) = x**2 - (2 - 1/4*k1) * x + (k0 + 3) / 4
          x = ONE - (k1 - SQRT(k11 - minima * 16._RK)) * ONEEIGHT
        else
          ! Third order Taylor expansion around the maximum local minima (x=1).
          c1 = ONE - ONEQUARTER * k1
          call find_a_cubic_root(c1, c0, x)
        end if
      else
        x = SQRT(MAX(ONE + SQRT(MAX(ONE - k0, 0._RK)), 0._RK))
        call newton(k1, k0, X, 5)
      end if
    else
      x = SQRT(MAX(ONE + SQRT(MAX(ONE - k0, 0._RK)), 0._RK))
      call newton(k1, k0, X, 50000)
    end if
  end subroutine solve_quartic_newton
!
  subroutine solve_quartic(k1, k0, X, Z)
    real(RK), intent(in)        :: k1, k0
    real(RK), intent(inout)     :: x
    character(2), intent(inout) :: Z
    real(RK)                    :: D, k11, k22
    real(RK)                    :: c1, c0, minima
    if (ABS(k1) < 1.E-1) then
      ! quasi biquadratic equations.
      k11 = k1 * k1
      k22 = 4._RK
      D = k0 - ONE
      minima = k1 + k0 - ONE
      if (minima >= ZERO) then
        ! multiple root
        x = ONE
        Z = 'A0'
      elseif (minima > -1.E-4) then
        c0 = k0 * ONEQUARTER - ONEQUARTER
        !print *, minima, c0
        if (ABS(c0) < 1.E-8) then
          ! Second order Taylor expansion around the maximum local minima (x=1).
          ! f2(x) = x**2 - (2 - 1/4*k1) * x + (k0 + 3) / 4
          x = ONE - (k1 - SQRT(k11 - minima * 16._RK)) * ONEEIGHT
          Z = 'A1'
        else
          ! Third order Taylor expansion around the maximum local minima (x=1).
          c1 = ONE - ONEQUARTER * k1
          call find_a_cubic_root(c1, c0, x)
          Z = 'A2'
        end if
      else
        x = SQRT(MAX(ONE + SQRT(MAX(ONE - k0, 0._RK)), 0._RK))
        call newton(k1, k0, X, 50)
        Z = 'A5'
      end if
    else
      ! Solve by resolvent cubic.
      call exact(k1, k0, X)
      Z = 'A6'
    end if
  end subroutine solve_quartic
!
  pure subroutine exact(k1, k0, X)
    real(RK), intent(in)    :: k1, k0
    real(RK), intent(inout) :: x
    real(RK)                :: c1, c0
    real(RK)                :: A, S, B
    ! Solve resolvent cubic y**3 - 2 * y**2 + (1 - k0) * y - k1**2 / 8
    ! Must be k1 > epsilon
    c1 = ONE - k0
    c0 = -k1 * k1 * ONEEIGHT
!   c1 = FOUR * (ONE - k0)
!   c0 = -k1 * k1
    call find_a_cubic_root(c1, c0, x)
    X = TWO * X
    if (X >= ZERO) then
      ! if X<0.0, solve the largest real root of x**4 + k2 * x**2 - k1 * x + k0 = 0.
      S = SQRT(X)
      A = FOUR - X
      B = TWO * k1 / S
      if (A > ABS(B)) then
        X = HALF * MAX(S + SQRT(A - B), -S + SQRT(A + B))
      elseif (A > B) then
        X = HALF * (S + SQRT(A - B))
      else
        X = HALF * (-S + SQRT(A + B))
      end if
    else
      ! if X<0.0, solve the smallest real root of x**4 + k2 * x**2 - k1 * x + k0 = 0.
      S = SQRT(-X)
      A = TWO + X
      B = -TWO * k1 / S
      if (A > ABS(B)) then
        X = -HALF * MAX(S - SQRT(A - B), -S - SQRT(A + B))
      elseif (A > B) then
        X = -HALF * (S - SQRT(A - B))
      else
        X = -HALF * (-S - SQRT(A + B))
      end if
    end if
    A = X * X
    B = -TWO + A
    S = A * B + k1 * X + k0
    x = x - S / (k1 + TWO * X * (A + B))
  end subroutine exact
!
  pure subroutine newton(k1, k0, x, maxiter)
    real(RK), intent(in)    :: k1, k0
    real(RK), intent(inout) :: x
    integer, intent(in)     :: maxiter
    real(RK)                :: s, xx, a, f, df
    integer                 :: i
    do i = 1, maxiter
      xx = x * x
      a = -TWO + xx
      f = a * xx + k1 * x + k0
      df = (k1 + (x + x) * (a + xx))
      if (ABS(f) < THRESHOLD) exit
      if (ABS(df) < THRESHOLD) exit
      s = f / df
      x = x - s
      if (ABS(s) < THRESHOLD) exit
    end do
  end subroutine newton
!
! find a positive root of monic cubic equation, x^3 - 2 * x^2 + k1 * x + k0 = 0
  pure subroutine find_a_cubic_root(k1, k0, x)
    real(RK), intent(in)    :: k1, k0
    real(RK), intent(inout) :: x
    real(RK), parameter     :: EIGHTNINE = (8.0_RK / 9.0_RK)
    real(RK), parameter     :: TWOTHIRD = (2.0_RK / 3.0_RK)
    real(RK), parameter     :: FOURTHIRD = (4.0_RK / 3.0_RK)
    real(RK), parameter     :: HALFSQRT3 = HALF * SQRT3
    real(RK)                :: R, Q, H, S
    X = TWOTHIRD
    R = -FOUR * (-(EIGHTNINE - k1) * TWOTHIRD + k0)
    Q = -FOURTHIRD * (k1 - FOURTHIRD)
    if (ABS(Q) < THRESHOLD) then
      if (ABS(R) < THRESHOLD) return
      X = X + SIGN(ONE, R) * (ONEQUARTER * R)**(-ONETHIRD)
    elseif (ABS(R) < THRESHOLD) then
      if (q > ZERO) then
        X = X + HALFSQRT3 * SQRT(q)
      else
        X = X - HALFSQRT3 * SQRT(-q)
      end if
    elseif (Q > ZERO) then
      S = SQRT(Q)
      Q = S * Q
      if (ABS(R) <= Q) then
        H = R / Q
        X = X + S * cos_acos(H)
      else
        H = ABS(Q / R)
        X = X + SIGN(ONE, R) * S * cosh_acosh(H)
      end if
    else
      S = SQRT(-Q)
      H = -R / (Q * S)
      X = X + S * sinh_asinh(H)
    end if
  end subroutine find_a_cubic_root
!
#include "cos_acos.f90"
!
end module mod_quatic

program main
  use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
  use mod_kinds
  use mod_params, only: RK, IK
  use mod_quatic
  implicit none
  call test2()
  call test3()
  call test4()
  call test5()
  call test6()
contains
!
  subroutine test2()
    integer(IK), parameter :: N = 5000000, M = 10
    real(RK)               :: A(N), B(N), C(N)
    real(RK)               :: time_begin_s, time_end_s, time1, time2, tot1, tot2
    real(RK)               :: merr1, merr2, nerr1, nerr2
    integer                :: i
    tot1 = 0.0_RK
    tot2 = 0.0_RK
    merr1 = 0.0_RK
    merr2 = 0.0_RK
    nerr1 = 0.0_RK
    nerr2 = 0.0_RK
    print'(A)', " cos_acos          implement       intrinsic"
    print'(A)', "            --------------------------------"
    do i = 1, M
      call RANDOM_NUMBER(A)
      !A = (A - 0.5_RK) * 1.99_RK
!
      call CPU_TIME(time_begin_s)
      B = cos_acos(A)
      call CPU_TIME(time_end_s)
      time1 = 1000 * (time_end_s - time_begin_s)
      tot1 = tot1 + time1
      merr1 = merr1 + mean_error((A - (4 * B**3 - 3 * B)) / MAX(ABS(A), 1.0E-8_RK))
      nerr1 = nerr1 + norm_error(A - (4 * B**3 - 3 * B))
!
      call CPU_TIME(time_begin_s)
      C = COS(ACOS(A) / 3.0_RK)
      call CPU_TIME(time_end_s)
      time2 = 1000 * (time_end_s - time_begin_s)
      tot2 = tot2 + time2
      merr2 = merr2 + mean_error((A - (4 * C**3 - 3 * C)) / MAX(ABS(A), 1.0E-8_RK))
      nerr2 = nerr2 + norm_error((A - (4 * C**3 - 3 * C)))
      print '(I12,2f16.3)', i, time1, time2
    end do
    print'(A)', "            --------------------------------"
    print'(12X,2f16.3)', tot1, tot2
    print'(A,2f16.9)', "  mean error", merr1 / M, merr2 / M
    print'(A,2f16.9)', "  norm error", nerr1 / M, nerr2 / M
    print *
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
  end subroutine test2
!
  subroutine test3()
    integer(IK), parameter :: N = 5000000, M = 10
    real(RK)               :: A(N), B(N), C(N)
    real(RK)               :: time_begin_s, time_end_s, time1, time2, tot1, tot2
    real(RK)               :: merr1, merr2, nerr1, nerr2
    integer                :: i
    tot1 = 0.0_RK
    tot2 = 0.0_RK
    merr1 = 0.0_RK
    merr2 = 0.0_RK
    nerr1 = 0.0_RK
    nerr2 = 0.0_RK
    print'(A)', " cosh_acosh        implement       intrinsic"
    print'(A)', "            --------------------------------"
    do i = 1, M
      call RANDOM_NUMBER(A)
      A = (A + 0.001_RK) / 1.001_RK
!
      call CPU_TIME(time_begin_s)
      B = cosh_acosh(A)
      call CPU_TIME(time_end_s)
      time1 = 1000 * (time_end_s - time_begin_s)
      tot1 = tot1 + time1
      merr1 = merr1 + mean_error((1.0_RK / A - (4 * B**3 - 3 * B)) / MAX(ABS(1.0_RK / A), 1.0E-8_RK))
      nerr1 = nerr1 + norm_error((1.0_RK / A - (4 * B**3 - 3 * B)))
!
      call CPU_TIME(time_begin_s)
      C = COSH(LOG(ABS(1.0_RK / A) + SQRT(1.0_RK / (A * A) - 1.0_RK)) / 3.0_RK)
      call CPU_TIME(time_end_s)
      time2 = 1000 * (time_end_s - time_begin_s)
      tot2 = tot2 + time2
      merr1 = merr1 + mean_error((1.0_RK / A - (4 * C**3 - 3 * C)) / MAX(ABS(1.0_RK / A), 1.0E-8_RK))
      nerr1 = nerr1 + norm_error((1.0_RK / A - (4 * C**3 - 3 * C)))
      print '(I12,2f16.3)', i, time1, time2
    end do
    print'(A)', "            --------------------------------"
    print'(12X,2f16.3)', tot1, tot2
    print'(A,2f16.9)', "  mean error", merr1 / M, merr2 / M
    print'(A,2f16.9)', "  norm error", nerr1 / M, nerr2 / M
    print *
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
  end subroutine test3
!
  subroutine test4()
    integer(IK), parameter :: N = 5000000, M = 10
    real(RK)               :: A(N), B(N), C(N)
    real(RK)               :: time_begin_s, time_end_s, time1, time2, tot1, tot2
    real(RK)               :: merr1, merr2, nerr1, nerr2
    integer                :: i
    tot1 = 0.0_RK
    tot2 = 0.0_RK
    merr1 = 0.0_RK
    merr2 = 0.0_RK
    nerr1 = 0.0_RK
    nerr2 = 0.0_RK
    print'(A)', " sinh_asinh        implement       intrinsic"
    print'(A)', "            --------------------------------"
    do i = 1, M
      call RANDOM_NUMBER(A)
      A = SIGN(1.0_RK, A) * EXP((A - 0.5_RK) * 10._RK)
!
      call CPU_TIME(time_begin_s)
      B = sinh_asinh(A)
      call CPU_TIME(time_end_s)
      time1 = 1000 * (time_end_s - time_begin_s)
      tot1 = tot1 + time1
      merr1 = merr1 + mean_error((A - (4 * B**3 + 3 * B)) / MAX(ABS(A), 1.0E-8_RK))
      nerr1 = nerr1 + norm_error((A - (4 * B**3 + 3 * B)))
!
      call CPU_TIME(time_begin_s)
      C = SIGN(1.0_RK, A) * SINH(LOG(ABS(A) + SQRT(A * A + 1.0_RK)) / 3.0_RK)
      call CPU_TIME(time_end_s)
      time2 = 1000 * (time_end_s - time_begin_s)
      tot2 = tot2 + time2
      merr2 = merr2 + mean_error((A - (4 * C**3 + 3 * C)) / MAX(ABS(A), 1.0E-8_RK))
      nerr2 = nerr2 + norm_error((A - (4 * C**3 + 3 * C)))
      print '(I12,2f16.3)', i, time1, time2
    end do
    print'(A)', "            --------------------------------"
    print'(12X,2f16.3)', tot1, tot2
    print'(A,2f16.9)', "  mean error", merr1 / M, merr2 / M
    print'(A,2f16.9)', "  norm error", nerr1 / M, nerr2 / M
    print *
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
  end subroutine test4
!
  subroutine test5()
    integer(IK), parameter :: N = 5000000, M = 10
    real(RK)               :: K2, K1(N), K0(N), X(N), Y(N)
    real(RK)               :: time_begin_s, time_end_s, time1, time2, tot1, tot2
    real(RK)               :: merr1, merr2, nerr1, nerr2
    integer                :: i
    tot1 = 0.0_RK
    tot2 = 0.0_RK
    merr1 = 0.0_RK
    merr2 = 0.0_RK
    nerr1 = 0.0_RK
    nerr2 = 0.0_RK
    print'(A)', " solve_quartic     implement          newton"
    print'(A)', "            --------------------------------"
    K2 = -2.0_RK
    do i = 1, M
      call RANDOM_NUMBER(Y)
      call RANDOM_NUMBER(K1)
      Y = 5 * EXP(-ABS(Y) * 0.02_RK)
      K1 = SIGN(1.0_RK, K1 - 0.5_RK) * EXP(-ABS(K1)) * 10._RK
      K0 = -Y**4 - K2 * Y**2 - K1 * Y
!
      call CPU_TIME(time_begin_s)
      call solve_quartic_tuned(K1, K0, X)
      call CPU_TIME(time_end_s)
      time1 = 1000 * (time_end_s - time_begin_s)
      tot1 = tot1 + time1
      merr1 = merr1 + mean_error(X**4 + K2 * X**2 + K1 * X + K0)
      nerr1 = nerr1 + norm_error(X**4 + K2 * X**2 + K1 * X + K0)
!
      call CPU_TIME(time_begin_s)
      call solve_quartic_newton(K1, K0, X)
      call CPU_TIME(time_end_s)
      time2 = 1000 * (time_end_s - time_begin_s)
      tot2 = tot2 + time2
      merr2 = merr2 + mean_error(X**4 + K2 * X**2 + K1 * X + K0)
      nerr2 = nerr2 + norm_error(X**4 + K2 * X**2 + K1 * X + K0)
      print '(I12,2f16.3)', i, time1, time2
    end do
    print'(A)', "            --------------------------------"
    print'(12X,2f16.3)', tot1, tot2
    print'(A,2f16.9)', "  mean error", merr1 / M, merr2 / M
    print'(A,2f16.9)', "  norm error", nerr1 / M, nerr2 / M
    print *
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
  end subroutine test5
!
  subroutine test6()
    real(RK)     :: X, Y, K0, K1, K2
    character(2) :: Z
    integer      :: i, j
    K2 = -2.0_RK
    do j = 1, 10
      do i = -100, 100
        K1 = SIGN(1, i) * 10 * EXP(-ABS(i) * 0.5_RK)
        X = 5 * EXP(-ABS(j) * 0.02_RK)
        Y = X
        K0 = -X**4 - K2 * X**2 - K1 * X
        call solve_quartic(K1, K0, X, Z)
#ifdef USE_REAL32
        if (ABS(X**4 + k2 * X**2 + k1 * X + k0) > 1.E-4) &
          & print '(A,2I4,*(f16.9))', Z, i, j, k2, k1, k0, Y, X, X**4 + k2 * X**2 + k1 * X + k0
#else
        if (ABS(X**4 + k2 * X**2 + k1 * X + k0) > 1.E-8) &
          & print '(A,2I4,*(f16.9))', Z, i, j, k2, k1, k0, Y, X, X**4 + k2 * X**2 + k1 * X + k0
#endif
      end do
    end do
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
  end subroutine test6
!
  pure function mean_error(e) result(res)
    real(RK), intent(in) :: e(:)
    real(RK)             :: res
    res = SQRT(SUM((e)**2) / SIZE(e))
  end function mean_error
!
  pure function norm_error(e) result(res)
    real(RK), intent(in) :: e(:)
    real(RK)             :: res
    res = MAXVAL(ABS(e))
  end function norm_error
end program main


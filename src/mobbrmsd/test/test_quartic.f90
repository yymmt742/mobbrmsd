module mod_quatic
  use mod_kinds
  use mod_params, only: RK, IK
  implicit none
  private
  public :: cos_acos
  public :: cosh_acosh
  public :: solve_quartic
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
  subroutine exact(k1, k0, X)
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
    call newton(k1, k0, X, 2)
  end subroutine exact
!
  subroutine newton(k1, k0, x, maxiter)
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
  subroutine find_a_cubic_root(k1, k0, x)
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
      H = MAX(R / (S * Q), -ONE)
      if (ABS(H) <= ONE) then
        X = X + S * COS(ONETHIRD * ACOS(H))
      else
        X = X + SIGN(ONE, R) * S * COSH(ONETHIRD * LOG(ABS(H) + SQRT(H * H - ONE)))
      end if
    else
      S = SQRT(-Q)
      H = -R / (Q * S)
      X = X + S * SINH(ONETHIRD * LOG(H + SQRT(H * H + ONE)))
    end if
  end subroutine find_a_cubic_root
!
#include "cos_acos.f90"
!
end module mod_quatic

program main
  use mod_kinds
  use mod_params, only: RK, IK
  use mod_quatic
  implicit none
  real(RK)     :: X, Y, K2, K1, K0
  character(2) :: Z
  real(RK)     :: h, a, b, c
  integer      :: i, j
  K2 = -2.0_RK
  print *, RK
  do i = 1, 1000
    h = i * 0.001
    a = cos_acos(h)
    b = cos_acos(-h)
    c = cosh_acosh(h)
    print'(*(f9.3))', h, 1 / h, a, b, c, 4 * a**3 - 3 * a, 4 * b**3 - 3 * b, 1 / (4 * c**3 - 3 * c)
  end do
  stop
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
    if (ABS(X**4 + k2 * X**2 + k1 * X + k0) > 1.E-6) &
      & print '(A,2I4,*(f16.9))', Z, i, j, k2, k1, k0, Y, X, X**4 + k2 * X**2 + k1 * X + k0
#endif
  end do
  end do
end program main


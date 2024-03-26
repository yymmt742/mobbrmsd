!
!| Utility functions for testing.
module mod_testutil
  use mod_params, only: D, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_rotation_matrix
  implicit none
  private
  public :: sample
  public :: covmat
  public :: gcov
  public :: SO3
  public :: eye
!
contains
!
  function sample(d, n) result(res)
    integer(IK), intent(in) :: d, n
    real(RK)                :: cnt(d)
    real(RK)                :: res(d, n)
    integer(IK)             :: i
    call RANDOM_NUMBER(res)
    cnt = SUM(res, 2) / n
    do concurrent(i=1:n)
      res(:, i) = res(:, i) - cnt
    end do
  end function sample
!
  function covmat(d, n) result(res)
    integer(IK), intent(in) :: d, n
    real(RK)                :: res(d, d)
    res(:, :) = MATMUL(sample(d, n), TRANSPOSE(sample(d, n)))
  end function covmat
!
  function gcov(d, n) result(res)
    integer(IK), intent(in) :: d, n
    real(RK)                :: res(d * d + 1)
    real(RK)                :: x(d, n), y(d, n)
    x = sample(d, n)
    y = sample(d, n)
    res(1) = SUM(x * x) + SUM(y * y)
    res(2:) = [MATMUL(y, TRANSPOSE(x))]
  end function gcov
!
  function SO3() result(res)
    real(RK) :: a(4), c, t, s, res(3, 3)
    call RANDOM_NUMBER(a)
    a(:3) = a(:3) / SQRT(DOT_PRODUCT(a(:3), a(:3)))
    a(4) = (a(4) + a(4)) * ACOS(0.0_RK)
    c = COS(a(4))
    t = ONE - c
    s = SIN(a(4))
    res(:, 1) = [c + t * a(1) * a(1), t * a(1) * a(2) - s * a(3), t * a(1) * a(3) + s * a(2)]
    res(:, 2) = [t * a(1) * a(2) + s * a(3), c + t * a(2) * a(2), t * a(2) * a(3) - s * a(1)]
    res(:, 3) = [t * a(1) * a(3) - s * a(2), t * a(2) * a(3) + s * a(1), c + t * a(3) * a(3)]
  end function SO3
!
  pure function eye(d) result(res)
    integer,intent(in) :: d
    real(RK)           :: res(d, d)
    integer            :: i, j
    do concurrent(j=1:d, i=1:d)
      res(i, j) = MERGE(ONE, ZERO, i == j)
    enddo
  end function eye
!
  pure function swp(m, n, per, map, sym, X) result(res)
    integer(IK), intent(in)        :: m, n, per(:), map(:), sym(:, :)
    real(RK), intent(in)           :: X(d, m, n)
    real(RK)                       :: tmp(d, m, n), res(d, m * n)
    integer(IK)                    :: i
    tmp = X
    do i = 1, SIZE(per)
      tmp(:, sym(:, map(i)), per(i)) = X(:, :, i)
    end do
    res = RESHAPE(tmp, [D, m * n])
  end function swp
!
  pure function sd(X, Y) result(res)
    real(RK), intent(in)    :: X(:, :), Y(:, :)
    real(RK)                :: C(D, D), R(D, D), W(100), res
    C = MATMUL(Y, TRANSPOSE(X))
    call estimate_rotation_matrix(SUM(X * X) + SUM(Y * Y), C, R, W)
    res = SUM(X**2) + SUM(Y**2) - 2 * SUM(C * R)
  end function sd
!
end module mod_testutil


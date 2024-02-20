!
!| Utility functions for testing.
module mod_testutil
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: sample
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
end module mod_testutil


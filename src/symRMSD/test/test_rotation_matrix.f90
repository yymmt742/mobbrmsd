program main
  use mod_params, only: D, setup_dimension, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_rotation_matrix
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: z
  real(RK)       :: E2(2, 2), E3(3, 3), E6(6, 6)
!
  E2 = eye(2)
  E3 = eye(3)
  E6 = eye(6)
!
  call z%init('test quartanion d=2')
  call setup_dimension(2)
  call test1(10, 10)
  call test1(20, 10)
  call test1(100, 10)
!
  call z%init('test quartanion d=3')
  call setup_dimension(3)
  call test1(3, 10)
  call test1(5, 10)
  call test1(10, 10)
  call test1(100, 10)
!
  call z%init('test kabsch d=6')
  call setup_dimension(6)
  call test1(1, 2)
  call test1(2, 2)
  call test1(3, 2)
  call test1(5, 4)
  call test1(100, 10)
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1(n, n_test)
    integer, intent(in)   :: n, n_test
    real(RK)              :: Y(d, n), X(d, n), cov(d, d), g
    real(RK)              :: rot(d, d), krot(d, d), sd, kd
    real(RK), allocatable :: w(:)
    integer               :: i
!
    allocate (w(MAX(worksize_rotation_matrix(), worksize_sdmin())))
!
    call RANDOM_NUMBER(X)
!
    do i = 1, N_TEST
      rot = SO(d)
      Y = MATMUL(rot, X)
      g = SUM(X * X) + SUM(Y * Y)
      cov = MATMUL(X, TRANSPOSE(Y))
      call estimate_rotation_matrix(g, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(krot, Y)], ZERO, 'X = YR   ')
      if (d <= n) call z%assert_almost_equal([MATMUL(rot, krot) - eye(d)], ZERO, 'S@RT = I ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye(d)], ZERO, 'R@RT = I ')
      call estimate_sdmin(g, cov, w)
      call z%assert_almost_equal(w(1), ZERO, 'sdmin=0  ')
    end do
!
    do i = 1, N_TEST
!
      call RANDOM_NUMBER(Y)
      cov = MATMUL(X, TRANSPOSE(Y))
      g = SUM(X**2) + SUM(Y**2)
      call estimate_rotation_matrix(g, cov, krot, w)
      call z%assert_greater_equal(SUM(cov * krot), SUM(cov * SO(d)), 'CR >= CQ ')
      call estimate_sdmin(g, cov, w)
!
      sd = SUM((X - MATMUL(krot, Y))**2)
      kd = SUM(cov * krot)
      kd = g - kd - kd
!
      call z%assert_almost_equal(w(1) / sd, ONE, 'sdmin-sd ', place=4)
      call z%assert_almost_equal(w(1) / kd, ONE, 'sdmin-kd ', place=4)
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye(d)], ZERO, 'R@RT = I ')
!
    end do
!
  end subroutine test1
!
  recursive function SO(d) result(res)
    integer(IK), intent(in) :: d
    real(RK)                :: res(d, d)
!
    if (d == 2) then
      block
        real(RK) :: a(1)
        call RANDOM_NUMBER(a)
        res(:, 1) = [COS(a(1)), -SIN(a(1))]
        res(:, 2) = [SIN(a(1)), COS(a(1))]
      end block
    elseif (d == 3) then
      block
        real(RK) :: a(3)
        call RANDOM_NUMBER(a)
        a = a / SQRT(DOT_PRODUCT(a, a))
        res(:, 1) = [a(1) * a(1), a(1) * a(2) - a(3), a(1) * a(3) + a(2)]
        res(:, 2) = [a(1) * a(2) + a(3), a(2) * a(2), a(2) * a(3) - a(1)]
        res(:, 3) = [a(1) * a(3) - a(2), a(2) * a(3) + a(1), a(3) * a(3)]
      end block
    elseif (d == 6) then
      block
        real(RK) :: tmp(d, d)
        res = ZERO
        res(1:2, 3:4) = SO(2)
        res(3:4, 1:2) = SO(2)
        res(5:6, 5:6) = SO(2)
!
        tmp = 0D0
        tmp(4:6, 1:3) = SO(3)
        tmp(1:3, 4:6) = -SO(3)
!
        res = MATMUL(res, tmp)
      end block
    else
      res = ZERO
    end if
!
  end function SO
!
end program main

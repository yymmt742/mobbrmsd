program main
  use mod_params, only: D, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_rotation
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: z
!
! call setup_dimension(2)
!
  call z%init('test rotation')
!
  call test1(1, 10)
  call test1(2, 10)
  call test1(3, 10)
  call test1(10, 10)
  call test1(20, 10)
  call test1(100, 10)
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1(n, n_test)
    integer, intent(in)   :: n, n_test
    real(RK)              :: Y(D, n), X(D, n), cov(D, D), g
    real(RK)              :: rot(D, D), krot(D, D), sd, kd
    real(RK), allocatable :: w(:)
    integer               :: i
!
    allocate (w(MAX(rotation_worksize(), sdmin_worksize())))
!
    X = sample(n)
!
    do i = 1, N_TEST
      rot = SO()
      Y = MATMUL(rot, X)
      g = SUM(X * X) + SUM(Y * Y)
      cov = MATMUL(X, TRANSPOSE(Y))
!
      call estimate_rotation(g, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(krot, Y)], ZERO, 'X = YR   ')
!
      if (D <= n) call z%assert_almost_equal([MATMUL(rot, krot) - eye()], ZERO, 'S@RT = I ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye()], ZERO, 'R@RT = I ')
!
      call estimate_sdmin(g, cov, w)
      call z%assert_almost_equal(w(1), ZERO, 'sdmin=0  ')
    end do
!
    do i = 1, N_TEST
!
      call RANDOM_NUMBER(Y)
      cov = MATMUL(X, TRANSPOSE(Y))
      g = SUM(X**2) + SUM(Y**2)
      call estimate_rotation(g, cov, krot, w)
      call z%assert_greater_equal(SUM(cov * krot), SUM(cov * SO()), 'CR >= CQ ')
      call estimate_sdmin(g, cov, w)
!
      sd = SUM((X - MATMUL(krot, Y))**2)
      kd = SUM(cov * krot)
      kd = g - kd - kd
!
      call z%assert_almost_equal(w(1) / sd, ONE, 'sdmin-sd ', place=4)
      call z%assert_almost_equal(w(1) / kd, ONE, 'sdmin-kd ', place=4)
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye()], ZERO, 'R@RT = I ')
!
    end do
!
  end subroutine test1
!
end program main

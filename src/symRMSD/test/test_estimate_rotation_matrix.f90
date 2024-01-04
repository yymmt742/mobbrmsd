program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_estimate_rotation_matrix
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
! call test1(2, 10, 100)
! call test1(2, 20, 100)
! call test1(2, 100, 100)
!
  call z%init('test quartanion d=3')
  call test1(3, 3, 10)
  call test1(3, 5, 10)
  call test1(3, 10, 10)
  call test1(3, 100, 10)
!
  call z%init('test kabsch d=6')
  call test2(1, 2)
  call test2(2, 2)
  call test2(3, 2)
  call test2(5, 4)
  call test2(100, 10)
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1(d, n, n_test)
    integer, intent(in)   :: d, n, n_test
    real(RK)              :: Y(d, n), X(d, n), cov(d, d), g
    real(RK)              :: rot(d, d), krot(d, d), sd, kd
    real(RK), allocatable :: w(:)
    integer               :: i
!
    call estimate_rotation_matrix(-d, X(1, 1), x, x, x)
    allocate (w(100*NINT(x(1, 1))))
!
    call RANDOM_NUMBER(X)
!
    do i = 1, N_TEST
!
      rot = SO(d)
      Y = MATMUL(rot, X)
      g = SUM(X * X) + SUM(Y * Y)
      cov = MATMUL(X, TRANSPOSE(Y))
      call estimate_rotation_matrix(d, g, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(krot, Y)], ZERO, 'X = YR  ')
      call z%assert_almost_equal([MATMUL(rot, krot) - eye(d)], ZERO, 'S@RT = I')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye(d)], ZERO, 'R@RT = I')
!
      call estimate_sdmin(d, g, cov, w)
      call z%assert_almost_equal(w(1), ZERO, 'sdmin=0 ')
!
    end do
!
    do i = 1, N_TEST
!
      call RANDOM_NUMBER(Y)
      cov = MATMUL(X, TRANSPOSE(Y))
      g = SUM(X**2) + SUM(Y**2)
      call estimate_rotation_matrix(d, g, cov, krot, w)
!
      call z%assert_greater_equal(SUM(cov * krot), SUM(cov * SO(d)), 'CR >= CQ ')
!
      call estimate_sdmin(d, g, cov, w)
!
      sd = SUM((X - MATMUL(krot, Y))**2)
      kd = SUM(cov * krot)
      kd = g - kd - kd
!
      call z%assert_almost_equal(w(1), sd, 'sdmin-sd ')
      call z%assert_almost_equal(w(1), kd, 'sdmin-kd ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye(d)], ZERO, 'R@RT = I ')
!
    end do
!
  end subroutine test1
!
  subroutine test2(n, n_test)
    integer, intent(in)   :: n, n_test
    real(RK)              :: Y(6, n), X(6, n), cov(6, 6)
    real(RK)              :: rot(6, 6), krot(6, 6), g
    real(RK), allocatable :: w(:)
    integer               :: i
!
    call estimate_rotation_matrix(-6, X(1, 1), x, x, x)
    allocate (w(NINT(x(1,1))))
!
    call random_number(X)
!
    do i = 1, N_TEST
!
      rot = SO6()
      Y = MATMUL(rot, X)
      cov = MATMUL(X, TRANSPOSE(Y))
      g = SUM(X * X) + SUM(Y * Y)
      call estimate_rotation_matrix(6, g, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(krot, Y)], 0D0, 'X = RY  ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - E6], 0D0, 'R@RT = I')
!
    enddo
!
  end subroutine test2
!
  function SO(d) result(res)
    integer(IK), intent(in) :: d
    real(RK)                :: a(3), res(d, d)
    call RANDOM_NUMBER(a)
    if (d == 2) then
      res(:, 1) = [COS(a(1)), -SIN(a(1))]
      res(:, 2) = [SIN(a(1)), COS(a(1))]
    elseif (d == 3) then
      a = a / SQRT(DOT_PRODUCT(a, a))
      res(:, 1) = [a(1) * a(1), a(1) * a(2) - a(3), a(1) * a(3) + a(2)]
      res(:, 2) = [a(1) * a(2) + a(3), a(2) * a(2), a(2) * a(3) - a(1)]
      res(:, 3) = [a(1) * a(3) - a(2), a(2) * a(3) + a(1), a(3) * a(3)]
    else
      res = ZERO
    end if
  end function SO
!
  function SO6() result(res)
    real(RK) :: res(6, 6), tmp(6, 6)
!
    res = 0D0
    res(1:2,3:4) = SO(2)
    res(3:4,1:2) = SO(2)
    res(5:6,5:6) = SO(2)
!
    tmp = 0D0
    tmp(4:6,1:3) = SO(3)
    tmp(1:3,4:6) = - SO(3)
!
    res = MATMUL(res, tmp)
!
  end function SO6
!
  pure function eye(d) result(res)
    integer,intent(in) :: d
    real(RK)           :: res(d, d)
    integer            :: i, j
    do concurrent(j=1:d, i=1:d)
      res(i, j) = MERGE(1D0, 0D0, i == j)
    enddo
  end function eye
!
end program main

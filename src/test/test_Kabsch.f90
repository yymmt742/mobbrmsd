program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_Kabsch
  use mod_unittest
  implicit none
  type(unittest) :: z
  real(RK)       :: E2(2, 2), E3(3, 3), E6(6, 6)
!
  E2 = eye(2)
  E3 = eye(3)
  E6 = eye(6)
!
  call z%init('test Kabsch d=3')
  call test1(1, 2)
  call test1(2, 2)
  call test1(3, 2)
  call test1(5, 4)
  call test1(100, 10)
!
  call z%init('test Kabsch d=6')
  call test2(1, 2)
  call test2(2, 2)
  call test2(3, 2)
  call test2(5, 4)
  call test2(100, 10)
!
  call z%init('test Kabsch row major n=3')
  call test3(1, 2)
  call test3(2, 2)
  call test3(3, 2)
  call test3(8, 4)
  call test3(100, 10)
!
  call z%init('test Kabsch row major n=6')
  call test4(1, 2)
  call test4(2, 2)
  call test4(3, 2)
  call test4(8, 4)
  call test4(100, 10)
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1(n, n_test)
    integer, intent(in) :: n, n_test
    integer, parameter  :: d = 3
    real(RK)            :: Y(d, n), X(d, n), cov(d, d)
    real(RK)            :: rot(d, d), krot(d, d)
    real(RK)            :: w(Kabsch_worksize(d))
    integer             :: i
!
    call random_number(X)
!
    do i=1,N_TEST
!
      rot = SO3()
      Y = MATMUL(rot, X)
      cov = MATMUL(X, TRANSPOSE(Y))
      call Kabsch(3, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(krot, Y)], 0D0, 'X = RY  ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - E3], 0D0, 'R@RT = I')
!
    enddo
!
  end subroutine test1
!
  subroutine test2(n, n_test)
    integer, intent(in) :: n, n_test
    real(RK)            :: Y(6, n), X(6, n), cov(6, 6)
    real(RK)            :: rot(6, 6), krot(6, 6)
    real(RK)            :: w(Kabsch_worksize(6))
    integer             :: i
!
    call random_number(X)
!
    do i = 1, N_TEST
!
      rot = SO6()
      Y = MATMUL(rot, X)
      cov = MATMUL(X, TRANSPOSE(Y))
      call Kabsch(6, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(krot, Y)], 0D0, 'X = RY  ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - E6], 0D0, 'R@RT = I')
!
    enddo
!
  end subroutine test2
!
  subroutine test3(d, n_test)
    integer, intent(in) :: d, n_test
    integer, parameter  :: n = 3
    real(RK)            :: Y(d, n), X(d, n), cov(n, n)
    real(RK)            :: rot(n, n), krot(n, n)
    real(RK)            :: w(Kabsch_worksize(n))
    integer             :: i
!
    call random_number(X)
!
    do i=1,N_TEST
!
      rot = SO3()
      Y = MATMUL(X, rot)
      cov = MATMUL(TRANSPOSE(X), Y)
      call Kabsch(n, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(Y, TRANSPOSE(krot))], 0D0, 'X = YR  ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye(n)], 0D0, 'R@RT = I')
!
    enddo
!
  end subroutine test3
!
  subroutine test4(d, n_test)
    integer, intent(in) :: d, n_test
    integer, parameter  :: n = 6
    real(RK)            :: Y(d, n), X(d, n), cov(n, n)
    real(RK)            :: rot(n, n), krot(n, n)
    real(RK)            :: w(Kabsch_worksize(n))
    integer             :: i
!
    call random_number(X)
    do i = 1, n
      X(:, i) = X(:, i) - SUM(X, 2)
    end do
!
    do i=1,N_TEST
!
      rot = SO6()
      Y = MATMUL(X, rot)
      cov = MATMUL(TRANSPOSE(X), Y)
      call Kabsch(n, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(Y, TRANSPOSE(krot))], 0D0, 'X = YR  ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye(n)], 0D0, 'R@RT = I')
!
    enddo
!
  end subroutine test4
!
  function SO2() result(res)
    real(RK) :: a(1), res(2, 2)
    call RANDOM_NUMBER(a)
    res(:, 1) = [ COS(a(1)),-SIN(a(1))]
    res(:, 2) = [ SIN(a(1)), COS(a(1))]
  end function SO2
!
  function SO3() result(res)
    real(RK) :: a(3), res(3, 3)
    call RANDOM_NUMBER(a)
    a = a / SQRT(DOT_PRODUCT(a, a))
    res(:, 1) = [a(1) * a(1),        a(1) * a(2) - a(3), a(1) * a(3) + a(2)]
    res(:, 2) = [a(1) * a(2) + a(3), a(2) * a(2),        a(2) * a(3) - a(1)]
    res(:, 3) = [a(1) * a(3) - a(2), a(2) * a(3) + a(1), a(3) * a(3)]
  end function SO3
!
  function SO6() result(res)
    real(RK) :: res(6, 6), tmp(6, 6)
!
    res = 0D0
    res(1:2,3:4) = SO2()
    res(3:4,1:2) = SO2()
    res(5:6,5:6) = SO2()
!
    tmp = 0D0
    tmp(4:6,1:3) = SO3()
    tmp(1:3,4:6) = - SO3()
!
    res = MATMUL(res, tmp)
!
  end function SO6
!
  function SO6_part() result(res)
    real(RK) :: res(6, 6)
!
    res = eye(6)
    !res(2:4, 2) = [0, 0, 1]
    !res(2:4, 3) = [1, 0, 0]
    !res(2:4, 4) = [0, 1, 0]
    res(3:5,3:5) = SO3()
!
  end function SO6_part
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

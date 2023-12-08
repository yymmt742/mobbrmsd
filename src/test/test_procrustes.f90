program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_pca
  use mod_procrustes
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
  call z%init('test Kabsch row major n=6, partial rotation')
  call test5(1, 1)
  call test5(8, 1)
  call test5(10, 4)
  call test5(300, 10)
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
    real(RK)            :: w(procrustes_worksize(d))
    integer             :: i
!
    call random_number(X)
!
    do i=1,N_TEST
!
      rot = SO3()
      Y = MATMUL(rot, X)
      cov = MATMUL(X, TRANSPOSE(Y))
      call procrustes(3, cov, krot, w)
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
    real(RK)            :: w(procrustes_worksize(6))
    integer             :: i
!
    call random_number(X)
!
    do i = 1, N_TEST
!
      rot = SO6()
      Y = MATMUL(rot, X)
      cov = MATMUL(X, TRANSPOSE(Y))
      call procrustes(6, cov, krot, w)
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
    real(RK)            :: w(procrustes_worksize(n))
    integer             :: i
!
    call random_number(X)
!
    do i=1,N_TEST
!
      rot = SO3()
      Y = MATMUL(X, rot)
      cov = MATMUL(TRANSPOSE(X), Y)
      call procrustes(n, cov, krot, w)
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
    real(RK)            :: w(procrustes_worksize(n))
    integer             :: i
!
    call random_number(X)
!
    do i=1,N_TEST
!
      rot = SO6()
      Y = MATMUL(X, rot)
      cov = MATMUL(TRANSPOSE(X), Y)
      call procrustes(n, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(Y, TRANSPOSE(krot))], 0D0, 'X = YR  ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye(n)], 0D0, 'R@RT = I')
!
    enddo
!
  end subroutine test4
!
  subroutine test5(d, n_test)
    integer, intent(in) :: d, n_test
    integer, parameter  :: n = 12
    real(RK)            :: Y(d, n), X(d, n), cov(n, n)
    real(RK)            :: rot(n, n), krot(n, n)
    real(RK)            :: w(procrustes_worksize(n)+10000)
    integer             :: i
!
    call random_number(X)
!
    do i = 1, N_TEST
      rot = SO12_part()
      Y = MATMUL(X, rot)
      cov = MATMUL(TRANSPOSE(Y), X)
      call get_permutation_matrix(d, n, X, Y, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(Y, krot)], 0D0, 'X = YR  ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye(n)], 0D0, 'R@RT = I')
    enddo
!
  end subroutine test5
!
  function SO2() result(res)
    real(RK) :: a(2), res(2, 2)
    call RANDOM_NUMBER(a)
    res(:, 1) = [ COS(a(1)),-SIN(a(1))]
    res(:, 2) = [ SIN(a(1)), COS(a(1))]
    if (a(2) < 0.5D0) res(2, :) = -res(2, :)
  end function SO2
!
  function SO3() result(res)
    real(RK) :: a(4), res(3, 3), P(3,3)
    call RANDOM_NUMBER(a)
    a(:3) = a(:3) / SQRT(DOT_PRODUCT(a(:3), a(:3)))
    res(:, 1) = [a(1) * a(1),        a(1) * a(2) - a(3), a(1) * a(3) + a(2)]
    res(:, 2) = [a(1) * a(2) + a(3), a(2) * a(2),        a(2) * a(3) - a(1)]
    res(:, 3) = [a(1) * a(3) - a(2), a(2) * a(3) + a(1), a(3) * a(3)]
    P(:,1) = [1,0,0]
    P(:,2) = [0,0,1]
    P(:,3) = [0,1,0]
    if (a(4) < 0.5D0) res = MATMUL(res, P)
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
  function SO12_part() result(res)
    real(RK) :: res(12, 12)
!
    res = eye(12)
    !res(2:4, 2) = [0, 0, 1]
    !res(2:4, 3) = [0, 1, 0]
    !res(2:4, 4) = [1, 0, 0]
    res(1:3,1:3) = SO3()
!
  end function SO12_part
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

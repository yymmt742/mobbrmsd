program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_group_permutation
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test group_permutation')
!
  call test1()
  call test2()
  call test3(2, 2, [1, 2, 2, 1])
  call test3(3, 6, [1, 2, 3, 2, 1, 3, 1, 3, 2, 2, 3, 1, 3, 1, 2, 3, 2, 1])
  call test3(6, 1, [1, 2, 3, 4, 5, 6])
  call test3(6, 1, [4, 5, 6, 1, 2, 3])
  call test3(6, 1, [4, 5, 6, 1, 2, 9])
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    type(group_permutation) :: g
    integer(IK), parameter  :: n = 13
    real(RK)                :: X(n), Y(2, n)
    integer                 :: s(n)
    integer                 :: i
!
    s = [2, 4, 6, 1, 3, 5, 9, 8, 7, 11, 13, 10, 12]
    g = group_permutation(RESHAPE(s, [n, 1]))
    do concurrent(i=1:n)
      X(i) = i
      Y(:, i) = i
    end do
!
    call group_permutation_swap(g%q, 1, 1, X)
    call group_permutation_swap(g%q, 1, 2, Y)
    call u%assert_equal(s, NINT(X, IK),       'swap    s = X ')
    call u%assert_equal(s, NINT(Y(1, :), IK), 'swap    s = Y1')
    call u%assert_equal(s, NINT(Y(2, :), IK), 'swap    s = Y2')
!
    call group_permutation_inverse(g%q, 1, 1, X)
    call group_permutation_inverse(g%q, 1, 2, Y)
    call u%assert_almost_equal(real([(i, i=1, n)], RK), X,       'inverse I = X ')
    call u%assert_almost_equal(real([(i, i=1, n)], RK), Y(1, :), 'inverse I = Y1')
    call u%assert_almost_equal(real([(i, i=1, n)], RK), Y(2, :), 'inverse I = Y2')
!
  end subroutine test1
!
  subroutine test2()
    type(group_permutation) :: g
    integer(IK), parameter  :: n = 8
    real(RK)                :: X(n), Y(2, n)
    integer                 :: s(n, 7), t(n)
    integer                 :: i
!
    s = RESHAPE([2, 3, 4, 5, 6, 7, 8, 1, &
      &          3, 4, 5, 6, 7, 8, 1, 2, &
      &          4, 5, 6, 7, 8, 1, 2, 3, &
      &          5, 6, 7, 8, 1, 2, 3, 4, &
      &          6, 7, 8, 1, 2, 3, 4, 5, &
      &          7, 8, 1, 2, 3, 4, 5, 6, &
      &          8, 1, 2, 3, 4, 5, 6, 7], [8, 7])
    t = [1, 2, 3, 4, 5, 6, 7, 8]
    g = group_permutation(s)
!
    do concurrent(i=1:n)
      X(i) = i
      Y(:, i) = i
    end do
!
    do i = 1, 7
      call group_permutation_swap(g%q, i, 1, X)
      call group_permutation_swap(g%q, i, 2, Y)
      call u%assert_equal(s(:, i), NINT(X, IK),       'swap    s = X ')
      call u%assert_equal(s(:, i), NINT(Y(1, :), IK), 'swap    s = Y1')
      call u%assert_equal(s(:, i), NINT(Y(2, :), IK), 'swap    s = Y2')
!
      call group_permutation_inverse(g%q, i, 1, X)
      call group_permutation_inverse(g%q, i, 2, Y)
      call u%assert_equal(t, NINT(X, IK),       'inverse I = X ')
      call u%assert_equal(t, NINT(Y(1, :), IK), 'inverse I = Y1')
      call u%assert_equal(t, NINT(Y(2, :), IK), 'inverse I = Y2')
    end do
!
  end subroutine test2
!
  subroutine test3(n, m, s)
    type(group_permutation) :: g
    integer(IK), intent(in) :: n, m, s(n * m)
    real(RK)                :: X(n), Y(2, n)
    integer                 :: t(n), v(n, m)
    integer                 :: i
!
    t = [(i, i=1, n)]
    v = RESHAPE(s, [n, m])
    g = group_permutation(v)
!
    do concurrent(i=1:n)
      X(i) = i
      Y(:, i) = i
    end do
!
    do i = 1, group_permutation_nsym(g%q) - 1
      call group_permutation_swap(g%q, i, 1, X)
      call group_permutation_swap(g%q, i, 2, Y)
      call u%assert_equal(v(:, i), NINT(X, IK), 'swap    s = X ')
      call u%assert_equal(v(:, i), NINT(Y(1, :), IK), 'swap    s = Y1')
      call u%assert_equal(v(:, i), NINT(Y(2, :), IK), 'swap    s = Y2')
!
      call group_permutation_inverse(g%q, i, 1, X)
      call group_permutation_inverse(g%q, i, 2, Y)
      call u%assert_equal(t, NINT(X, IK), 'inverse I = X ')
      call u%assert_equal(t, NINT(Y(1, :), IK), 'inverse I = Y1')
      call u%assert_equal(t, NINT(Y(2, :), IK), 'inverse I = Y2')
    end do
!
  end subroutine test3
!
end program main


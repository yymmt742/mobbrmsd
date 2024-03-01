program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_group_permutation
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test group_permutation')
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    type(group_permutation_tuple) :: g
    integer(IK), parameter  :: n = 13
    real(RK)                :: X(n), Y(2, n)
    integer                 :: s(n)
    integer                 :: i
!
    s = [2, 4, 6, 1, 3, 5, 9, 8, 7, 11, 13, 10, 12]
    g = group_permutation_tuple(RESHAPE(s, [n, 1]))
    do concurrent(i=1:n)
      X(i) = i
      Y(:, i) = i
    end do
!
    call group_permutation_swap(g%t, g%w, 1, 1, X)
    call group_permutation_swap(g%t, g%w, 1, 2, Y)
    call u%assert_equal(s, NINT(X, IK),       'swap    s = X ')
    call u%assert_equal(s, NINT(Y(1, :), IK), 'swap    s = Y1')
    call u%assert_equal(s, NINT(Y(2, :), IK), 'swap    s = Y2')
!
    call group_permutation_inverse(g%t, g%w, 1, 1, X)
    call group_permutation_inverse(g%t, g%w, 1, 2, Y)
    call u%assert_almost_equal(real([(i, i=1, n)], RK), X,       'inverse I = X ')
    call u%assert_almost_equal(real([(i, i=1, n)], RK), Y(1, :), 'inverse I = Y1')
    call u%assert_almost_equal(real([(i, i=1, n)], RK), Y(2, :), 'inverse I = Y2')
!
!   g = group_permutation_tuple([1,2,3,4,5,6,7,8,9,10,11,12])
!   g = group_permutation_tuple([12,1,2,3,4,5,6,7,8,9,10,11])
!   g = group_permutation_tuple([7,8,9,10,11,12,1,2,3,4,5,6])
!
  end subroutine test1
!
end program main


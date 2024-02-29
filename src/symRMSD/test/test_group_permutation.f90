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
    s = [2, 4, 6, 1, 3, 5, 9, 8, 7, 11, 13, 10, 12]
    g = group_permutation_tuple(s)
    do concurrent(i=1:n)
      X(i) = i
      Y(:, i) = i
    end do
    call group_permutation_swap(g%t, g%w, 1, X)
    call group_permutation_swap(g%t, g%w, 2, Y)
!
    call u%assert_almost_equal(real(s, RK), X, 's = X')
    call u%assert_almost_equal(real(s, RK), Y(1, :), 's = Y1')
    call u%assert_almost_equal(real(s, RK), Y(2, :), 's = Y2')
!
    g = group_permutation_tuple([1,2,3,4,5,6,7,8,9,10,11,12])
    g = group_permutation_tuple([12,1,2,3,4,5,6,7,8,9,10,11])
    g = group_permutation_tuple([7,8,9,10,11,12,1,2,3,4,5,6])
!
!   s = [3, 2, 1, 6, 5, 4, 8, 7, 9,10, 11, 13, 12]
!   g = group_permutation(s)
!   do concurrent(i=1:n)
!     X(i) = i
!     Y(:, i) = i
!   end do
!   call g%swap(1, X)
!   call g%swap(2, Y)
!   call u%assert_almost_equal(real(s, RK), X, 's = X')
!   call u%assert_almost_equal(real(s, RK), Y(1, :), 's = Y1')
!   call u%assert_almost_equal(real(s, RK), Y(2, :), 's = Y2')
!   call u%assert_equal(g%nfree(), 8, 'nfree')
!   call u%assert_equal(g%free_indices(), [1,3,4,6,7,8,12,13], 'free_indices')
!
  end subroutine test1
!
end program main


program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_diagonalize
  use mod_unittest
  implicit none
  type(unittest) :: z
!
  call z%init('test diagonalize')
  call test1()
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter :: N_TEST=10
    real(RK)           :: X(3, 9), Y(3, 9), E(3, 3)
    real(RK)           :: u(3, 3), v(3), w(100)
    integer            :: i
!
    E(:, 1) = [1, 0, 0]
    E(:, 2) = [0, 1, 0]
    E(:, 3) = [0, 0, 1]
!
!   do i=1,N_TEST
!
!     call random_number(X)
!     call pca(.FALSE., 3, 9, X, u, v, w)
!
!     call z%assert_almost_equal([MATMUL(TRANSPOSE(u), u) - E], 0D0, 'U@UT=I')
!     call z%assert_almost_equal([MATMUL(u, TRANSPOSE(u)) - E], 0D0, 'U@UT=I')
!     Y = MATMUL(transpose(u), X)
!     Y(:,:3) = RESHAPE([v(1), 0D0, 0D0, 0D0, v(2), 0D0, 0D0, 0D0, v(3)], [3, 3]) - MATMUL(Y, TRANSPOSE(Y))
!     call z%assert_almost_equal([Y(:,:3)], 0D0, 'U@S@UT=X')
!
!   enddo
!
  end subroutine test1
!
end program main

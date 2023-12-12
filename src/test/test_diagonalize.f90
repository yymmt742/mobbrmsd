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
    integer, parameter :: n = 10
    integer, parameter :: m = 2
    real(RK)           :: X(n, n, m)
    real(RK)           :: U(n, n)
    real(RK)           :: L(n, m)
    real(RK)           :: rwrk(1000)
    integer            :: i
!
    call random_number(X)
    do i = 1, m
      X(:, n - 1, i) = X(:, n - 2, i)
      X(:, n, i) = X(:, n - 1, i)
      X(:, :, i) = MATMUL(X(:, :, i), TRANSPOSE(X(:, :, i)))
    enddo
print'(10f6.1)', MATMUL(X(:,:,1),X(:,:,2)) + MATMUL(X(:,:,2),X(:,:,1))
return
    call simultaneous_diagonalize(n, m, X, U, L, rwrk)
    do i = 1, m
      print'(10f6.1)', MATMUL(MATMUL(TRANSPOSE(U), X(:, :, i)), U)
      print*
      print'(10f6.1)', L(:,i)
      print*
    end do
!
!   do i=1,N_TEST
!
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

program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_pca
  use mod_unittest
  implicit none
  type(unittest) :: z
!
  call z%init('test pca d=3')
  call test1()
  call test2()
!
! call z%init('test svd d=100')
! call test2(2)
! call test2(4)
! call test2(5)
! call test2(6)
! call test2(7)
! call test2(8)
! call test2(9)
! call test2(10)
! call test2(20)
! call test2(50)
! call test2(100)
!
! call z%init('test svd covariance matrix')
! call test3(3, 100)
! call test3(9, 100)
! call test3(9, 3)
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
    do i=1,N_TEST
!
      call random_number(X)
      call pca(.FALSE., 3, 9, X, u, v, w)
!
      call z%assert_almost_equal([MATMUL(TRANSPOSE(u), u) - E], 0D0, 'U@UT=I')
      call z%assert_almost_equal([MATMUL(u, TRANSPOSE(u)) - E], 0D0, 'U@UT=I')
      Y = MATMUL(transpose(u), X)
      Y(:,:3) = RESHAPE([v(1), 0D0, 0D0, 0D0, v(2), 0D0, 0D0, 0D0, v(3)], [3, 3]) - MATMUL(Y, TRANSPOSE(Y))
      call z%assert_almost_equal([Y(:,:3)], 0D0, 'U@S@UT=X')
!
    enddo
!
  end subroutine test1
!
  subroutine test2()
    integer, parameter :: N_TEST=10
    real(RK)           :: X(3, 9), Y(3, 9), E(9, 9)
    real(RK)           :: u(9, 9), v(9), w(100)
    real(RK)           :: m(9, 9)
    integer            :: i, j, k
!
    do concurrent(i=1:9, j=1:9)
      E(i, j) = MERGE(1, 0, i == j)
    enddo
!
    do i=1,N_TEST
!
      call random_number(X)
      call pca(.TRUE., 3, 9, X, u, v, w)
!
      call z%assert_almost_equal([MATMUL(TRANSPOSE(u), u) - E], 0D0, 'U@UT=I')
      call z%assert_almost_equal([MATMUL(u, TRANSPOSE(u)) - E], 0D0, 'U@UT=I')
      Y = MATMUL(X, u)
      do concurrent(j=1:9, k=1:9)
        m(j, k) = MERGE(v(j), 0D0, j == k)
      enddo
      m = m - MATMUL(TRANSPOSE(Y), Y)
      call z%assert_almost_equal([m], 0D0, 'U@S@UT=X')
!
    enddo
!
  end subroutine test2
!
! subroutine test2(d)
!   integer,intent(in)    :: d
!   integer, parameter    :: N_TEST = 10
!   real(RK)              :: Y(d * d), X(d * d), E(d, d), Q(d, d)
!   real(RK)              :: s(d), u(d, d), vt(d, d)
!   real(RK), allocatable :: w(:)
!   integer               :: lw, i, j
!
!   lw = svd_worksize(d)
!   allocate(w(lw))
!
!   Q = 0D0
!   do concurrent(j=1:d,i = 1:d)
!     E(i, j) = MERGE(1D0, 0D0, i == j)
!   end do
!
!   do i = 1, N_TEST
!
!     call RANDOM_NUMBER(X)
!     Y = X; call svd(d, Y, s, u, vt, w)
!
!     Y = [MATMUL(u, TRANSPOSE(u)) - E]
!     call z%assert_almost_equal(Y, 0D0, 'U@UT=I')
!
!     Y = [MATMUL(TRANSPOSE(vt), vt) - E]
!     call z%assert_almost_equal(Y, 0D0, 'V@VT=I')
!
!     do j = 1, d
!       Q(j, j) = s(j)
!     end do
!     Y = [MATMUL(MATMUL(u, Q), vt)] - X
!     call z%assert_almost_equal(Y, 0D0, 'U@S@VT=X')
!
!   end do
!
! end subroutine test2
!
! subroutine test3(d, n)
!   integer,intent(in)    :: d, n
!   real(RK)              :: V(d, n), W(d, n), X(d, d), Y(d, d), E(d, d), Q(d, d)
!   real(RK)              :: s(d), u(d, d), vt(d, d)
!   real(RK), allocatable :: wk(:)
!   integer               :: lw, i, j
!
!   lw = svd_worksize(d)
!   allocate(wk(lw))
!
!   Q = 0D0
!   do concurrent(j=1:d,i = 1:d)
!     E(i, j) = MERGE(1D0, 0D0, i == j)
!   end do
!
!   call RANDOM_NUMBER(V)
!   call RANDOM_NUMBER(W)
!   X = MATMUL(V, TRANSPOSE(W))
!   Y = X
!   call svd(d, Y, s, u, vt, wk)
!
!   Y = MATMUL(u, TRANSPOSE(u)) - E
!   call z%assert_almost_equal([Y], 0D0, 'U@UT=I')
!
!   Y = MATMUL(TRANSPOSE(vt), vt) - E
!   call z%assert_almost_equal([Y], 0D0, 'V@VT=I')
!
!   do j = 1, d
!     Q(j, j) = s(j)
!   end do
!   Y = MATMUL(MATMUL(u, Q), vt) - X
!   call z%assert_almost_equal([Y], 0D0, 'U@S@VT=X')
!
! end subroutine test3
!
end program main

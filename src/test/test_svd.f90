program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_svd
  use mod_unittest
  implicit none
  type(unittest) :: z
!
  call test1()
  call test2()
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter :: N_TEST=10
    real(RK)           :: Y(9), X(9), E(3, 3)
    real(RK)           :: s(3), u(3, 3), vt(3, 3), w(100)
    integer            :: i
!
    call z%init('test svd d=3')
!
    E(:, 1) = [1, 0, 0]
    E(:, 2) = [0, 1, 0]
    E(:, 3) = [0, 0, 1]
!
    do i=1,N_TEST
!
      call random_number(X)
      Y=X ; call svd(3, Y, s, u, vt, w)
!
      Y = [MATMUL(u, TRANSPOSE(u)) - E]
      call z%assert_almost_equal(Y, 0D0, 'U@UT=I')
!
      Y = [MATMUL(TRANSPOSE(vt), vt) - E]
      call z%assert_almost_equal(Y, 0D0, 'V@VT=I')
!
      Y = [MATMUL(MATMUL(u, RESHAPE([s(1), 0D0, 0D0, 0D0, s(2), 0D0, 0D0, 0D0, s(3)], [3, 3])), vt)] - X
      call z%assert_almost_equal(Y, 0D0, 'U@S@VT=X')
!
    enddo
!
  end subroutine test1
!
  subroutine test2()
    integer, parameter    :: N_TEST = 10
    integer, parameter    :: d = 10
    real(RK)              :: Y(d * d), X(d * d), E(d, d), Q(d, d)
    real(RK)              :: s(d), u(d, d), vt(d, d)
    real(RK), allocatable :: w(:)
    integer               :: lw, i, j
!
    call z%init('test svd d=100')
!
    lw = svd_worksize(d)
    allocate(w(lw))
!
    Q = 0D0
    E = 0D0
    do i = 1, d
      E(i, i) = 1D0
    end do
!
    do i = 1, N_TEST
!
      call RANDOM_NUMBER(X)
      Y = X; call svd(d, Y, s, u, vt, w)
!
      Y = [MATMUL(u, TRANSPOSE(u)) - E]
      call z%assert_almost_equal(Y, 0D0, 'U@UT=I')
!
      Y = [MATMUL(TRANSPOSE(vt), vt) - E]
      call z%assert_almost_equal(Y, 0D0, 'V@VT=I')
!
      do j = 1, d
        Q(j, j) = s(j)
      end do
      Y = [MATMUL(MATMUL(u, Q), vt)] - X
      call z%assert_almost_equal(Y, 0D0, 'U@S@VT=X')
!
    end do
!
  end subroutine test2
!
end program main

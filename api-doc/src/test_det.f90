program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_det
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call test1()
  call test2()
  call test3()
  call test4()
  call test5()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter :: N_TEST=10
    real(RK)           :: W(1), X(1)
    integer            :: i
!
    call u%init('test det d=1')
!
    do i=1,N_TEST
      call random_number(X)
      W=X ; call det(1, W)
      call u%assert_almost_equal(W(1), X(1), 'det d=1')
      W=X ; call det_sign(1, W)
      call u%assert_almost_equal(W(1), SIGN(ONE, X(1)), 'det_sign d=1')
    enddo
!
  end subroutine test1
!
  subroutine test2()
    integer, parameter :: N_TEST=10
    real(RK)           :: W(4), X(4)
    real(RK)           :: da
    integer            :: i
!
    call u%init('test det d=2')
!
    do i=1,N_TEST
      call random_number(X)
      da = X(1)*X(4)-X(2)*X(3)
      W=X ; call det(2, W)
      call u%assert_almost_equal(W(1),            da, 'det d=2')
      W=X ; call det_sign(2, W)
      call u%assert_almost_equal(W(1), SIGN(ONE, da), 'det_sign d=2')
    enddo
!
  end subroutine test2
!
  subroutine test3()
    integer, parameter :: N_TEST=10
    real(RK)           :: W(9), X(9)
    real(RK)           :: da
    integer            :: i
!
    call u%init('test det d=3')
!
    do i=1,N_TEST
      call random_number(X)
      da = X(1)*X(5)*X(9)+X(2)*X(6)*X(7)+X(3)*X(4)*X(8) &
          -X(3)*X(5)*X(7)-X(2)*X(4)*X(9)-X(1)*X(6)*X(8)
      W=X ; call det(3, W)
      call u%assert_almost_equal(W(1),            da, 'det d=3')
      W=X ; call det_sign(3, W)
      call u%assert_almost_equal(W(1), SIGN(ONE, da), 'det_sign d=3')
    enddo
!
  end subroutine test3
!
  subroutine test4()
    integer, parameter :: N_TEST=20
    real(RK)           :: W(16), Y(16), X(16)
    integer            :: i
!
    call u%init('test det d=4')
!
    do i=1,N_TEST
      call random_number(X)
      W=X ; call det(4, W)
      Y=X ; call det_sign(4, Y)
      call u%assert_almost_equal(Y(1), SIGN(ONE, W(1)), 'det_sign d=4')
    enddo
!
  end subroutine test4
!
  subroutine test5()
    integer, parameter :: N_TEST=20
    integer, parameter :: d = 100
    real(RK)           :: W(d*d), Y(d*d), X(d*d)
    integer            :: i
!
    call u%init('test det d=100')
!
    do i=1,N_TEST
      call random_number(X)
      W=X ; call det(d, W)
      Y=X ; call det_sign(d, Y)
      call u%assert_almost_equal(Y(1), SIGN(ONE, W(1)), 'det_sign d=100')
    enddo
!
  end subroutine test5
!
end program main

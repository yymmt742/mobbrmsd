program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_Kabsch
  use mod_unittest
  implicit none
  type(unittest) :: z
!
  call test1()
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter :: N_TEST=10
    real(RK)           :: Y(18), X(18), rot(3,3)
    integer            :: i
!
    call z%init('test Kabsch d=3')
!
    do i=1,N_TEST
      call random_number(X)
      call random_number(Y)
      call Kabsch(3, 6, X, Y, rot)
      print'(3f9.3)',rot
!     call z%assert_almost_equal(res, X(1), 'det d=1')
!     call z%assert_almost_equal(sig, SIGN(ONE, X(1)), 'det_sign d=1')
    enddo
!
  end subroutine test1
!
end program main

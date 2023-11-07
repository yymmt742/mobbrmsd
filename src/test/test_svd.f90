program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_svd
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
    real(RK)           :: Y(9), X(9)
    real(RK)           :: s(3), u(3,3), vt(3,3), w(100)
    integer            :: i
!
    call z%init('test svd d=1')
!
    do i=1,N_TEST
      call random_number(X)
      Y=X ; call svd(3, Y, s, u(1,1), vt(1,1), w)
      print'(3f9.3)',matmul(u,transpose(u)),matmul(vt,transpose(vt)),s
      print'(3f9.3)',[matmul(matmul(u,reshape([s(1),0d0,0d0,0d0,s(2),0d0,0d0,0d0,s(3)],[3,3])),vt)]-X
      print*
!     call z%assert_almost_equal(res, X(1), 'det d=1')
!     call z%assert_almost_equal(sig, SIGN(ONE, X(1)), 'det_sign d=1')
    enddo
!
  end subroutine test1
!
end program main

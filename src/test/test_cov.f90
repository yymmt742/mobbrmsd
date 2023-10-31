program main
  use mod_params, only: RK, IK, RZERO
  use mod_cov
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
!
    real(RK), parameter    :: X(12) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    real(RK), parameter    :: Y(12) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    real(RK)               :: C1(1, 1)
    real(RK)               :: C2(2, 2)
    real(RK)               :: C3(3, 3)
!
    call u%init('test cov')
!
    C1 = RZERO
    call cov(1, 12, X, Y, C1)
    call u%assert_almost_equal([C1], [650D0], 'cov 1d')
!
    C2 = RZERO
    call cov(2, 6, X, Y, C2)
    call u%assert_almost_equal([C2], [286D0, 322D0, &
   &                                   322D0, 364D0], &
   &                                   'cov 2d')
!
    C3 = RZERO
    call cov(3, 4, X, Y, C3)
    call u%assert_almost_equal([C3], [166D0, 188D0, 210D0,&
   &                                   188D0, 214D0, 240D0,&
   &                                   210D0, 240D0, 270D0],&
   &                                   'cov 3d')
!
  end subroutine test1
!
end program main

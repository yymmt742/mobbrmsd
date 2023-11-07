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
    call u%assert_almost_equal([C1], [650D0], 'cov  1d')
    C1 = RZERO
    call cov1(1, 12, X, Y, C1)
    call u%assert_almost_equal([C1], [650D0], 'cov1 1d')
    C1 = RZERO
    call cov2(1, 12, X, Y, C1)
    call u%assert_almost_equal([C1], [650D0], 'cov2 1d')
    C1 = RZERO
    call cov3(1, 12, X, Y, C1)
    call u%assert_almost_equal([C1], [650D0], 'cov3 1d')
!
    C2 = RZERO
    call cov(2, 6, X, Y, C2)
    call u%assert_almost_equal([C2], [286D0, 322D0, &
   &                                   322D0, 364D0], &
   &                                   'cov  2d')
    C2 = RZERO
    call cov1(2, 6, X, Y, C2)
    call u%assert_almost_equal([C2], [286D0, 322D0, &
   &                                   322D0, 364D0], &
   &                                   'cov1 2d')
    C2 = RZERO
    call cov2(2, 6, X, Y, C2)
    call u%assert_almost_equal([C2], [286D0, 322D0, &
   &                                   322D0, 364D0], &
   &                                   'cov2 2d')
    C2 = RZERO
    call cov3(2, 6, X, Y, C2)
    call u%assert_almost_equal([C2], [286D0, 322D0, &
   &                                   322D0, 364D0], &
   &                                   'cov3 2d')
!
    C3 = RZERO
    call cov(3, 4, X, Y, C3)
    call u%assert_almost_equal([C3], [166D0, 188D0, 210D0,&
   &                                   188D0, 214D0, 240D0,&
   &                                   210D0, 240D0, 270D0],&
   &                                   'cov  3d')
    C3 = RZERO
    call cov1(3, 4, X, Y, C3)
    call u%assert_almost_equal([C3], [166D0, 188D0, 210D0,&
   &                                   188D0, 214D0, 240D0,&
   &                                   210D0, 240D0, 270D0],&
   &                                   'cov1 3d')
    C3 = RZERO
    call cov2(3, 4, X, Y, C3)
    call u%assert_almost_equal([C3], [166D0, 188D0, 210D0,&
   &                                  188D0, 214D0, 240D0,&
   &                                  210D0, 240D0, 270D0],&
   &                                  'cov2 3d')
    C3 = RZERO
    call cov3(3, 4, X, Y, C3)
    call u%assert_almost_equal([C3], [166D0, 188D0, 210D0,&
   &                                  188D0, 214D0, 240D0,&
   &                                  210D0, 240D0, 270D0],&
   &                                   'cov3 3d')
!
  end subroutine test1
!
  pure subroutine cov1(d, n, x, y, res)
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                RK => REAL64,  &
    &                IK => INT32
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: x(d, n), y(d, n)
    real(RK), intent(inout) :: res(d, d)
!
    res = MATMUL(x, TRANSPOSE(y))
!
  end subroutine cov1
!
  pure subroutine cov2(d, n, x, y, res)
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                RK => REAL64,  &
    &                IK => INT32
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: x(d, n), y(d, n)
    real(RK), intent(inout) :: res(d, d)
    integer(IK)             :: i, j, k
!
    do k = 1, n
      do concurrent(j=1:d, i=1:d)
        res(i, j) = res(i, j) + x(i, k) * y(j, k)
      end do
    end do
!
  end subroutine cov2
!
  pure subroutine cov3(d, n, x, y, res)
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                RK => REAL64,  &
    &                IK => INT32
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: x(d, n), y(d, n)
    real(RK), intent(inout) :: res(d, d)
!
    interface
      include 'dgemm.h'
    end interface
!
    call DGEMM('N', 'T', d, d, n, 1D0, x(1:,1:), d, y(1:,1:), d, 0D0, res(1:,1:), d)
!
  end subroutine cov3
!
end program main

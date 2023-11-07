  pure subroutine cov(d, n, x, y, res)
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                RK => REAL64,  &
    &                IK => INT32
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: x(d, n), y(d, n)
    real(RK), intent(inout) :: res(d, d)
!
    res = MATMUL(x, TRANSPOSE(y))
!
  end subroutine cov

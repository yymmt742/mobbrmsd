  pure subroutine cov(d, n, x, y, res)
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
  end subroutine cov

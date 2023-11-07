  pure subroutine cov(d, n, x, y, res)
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
  end subroutine cov

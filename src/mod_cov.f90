module mod_cov
  use mod_params
  implicit none
  private
  public :: cov
!
contains
!
  pure subroutine cov(d, n, x, y, res)
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
!
end module mod_cov

module mod_cov
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  implicit none
  private
  public :: cov
!
  interface
    include 'dgemm.h'
  end interface
!
contains
!
  pure subroutine cov(d, n, x, y, res)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: x(d, n), y(d, n)
    real(RK), intent(inout) :: res(d, d)
!
    call DGEMM('N', 'T', d, d, n, ONE, x(1:,1:), d, y(1:,1:), d, ZERO, res(1:,1:), d)
!
  end subroutine cov
!
end module mod_cov

module mod_rot
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  implicit none
  private
  public :: rot
!
  interface
    include 'dgemm.h'
  end interface
!
contains
!
  pure subroutine rot(d, n, u, x, res)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: u(d, d), x(d, n)
    real(RK), intent(inout) :: res(d, n)
!
    call DGEMM('N', 'N', d, n, d, ONE, u(1:,1:), d, x(1:,1:), d, ZERO, res(1:,1:), d)
!
  end subroutine rot
!
end module mod_rot

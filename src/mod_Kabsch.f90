module mod_Kabsch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_cov
  use mod_det
  implicit none
  private
  public :: rot
!
contains
!
  pure subroutine rot_Kabsch(d, n, x, y, res)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: x(d, n), y(d, n)
    real(RK), intent(inout) :: res(d, d)
!
  end subroutine rot_Kabsch
!
end module mod_Kabsch

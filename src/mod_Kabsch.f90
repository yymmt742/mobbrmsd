module mod_Kabsch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_svd
  use mod_cov
  use mod_det
  implicit none
  private
  public :: Kabsch
!
contains
!
  pure subroutine Kabsch(d, n, x, y, rot)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: x(d, n), y(d, n)
    real(RK), intent(inout) :: rot(d, d)
    real(RK)                :: m(d, d)
    real(RK)                :: u(d, d), vt(d, d), s(d), w(4*d**2)
!
    call cov(d, n, x, y, m)
    call svd(d, m, s, u, vt, w)
    m = MATMUL(u, TRANSPOSE(vt))
    call det_sign(d, m, s(1))
    vt(:, d) = vt(:, d) * s(1)
    rot = MATMUL(vt, u)
!
  end subroutine Kabsch
!
end module mod_Kabsch

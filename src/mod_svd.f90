module mod_svd
  use mod_params, only : IK, RK, ONE => RONE, ZERO => RZERO
  implicit none
  private
  public :: svd
!
  interface
    include 'dgesvd.h'
  end interface
!
contains
!
   pure subroutine svd(d, x, s, u, vt, w)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: x(d, d)
    real(RK), intent(inout) :: s(d), u(d, d), vt(d, d), w(:)
    integer(IK)             :: info
!
      call DGESVD('A', 'A', d, d, x, d, s, u, d, vt, d, w, SIZE(w), info)
!
  end subroutine svd
!
end module mod_svd

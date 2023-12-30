!| Calculate the rotation matrix that minimizes |X-RY|^2 using the Kabsch-Umeyama algorithm.
!  Here, RR^T=I and det(R)=1 are satisfied.
module mod_Kabsch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_svd
  use mod_det
  implicit none
  private
  public :: Kabsch_worksize, Kabsch
!
  real(RK), parameter :: DEF_threshold = 1.0E-8_RK
!
  interface
    include 'dgemm.h'
  end interface
!
contains
!
!| Calculate work array size for d*d matrix.
  pure elemental function Kabsch_worksize(d) result(res)
    integer(IK), intent(in)           :: d
    !! matrix collumn dimension.
    integer(IK)                       :: res
!
    if (d < 1) then
      res = 0
    else
      res = svd_worksize(d) + d * d * 3 + d
    end if
!
  end function Kabsch_worksize
!
!| Calculate the rotation matrix from covariance matrix.
  pure subroutine Kabsch(d, cov, rot, w, ldcov)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    real(RK), intent(in)          :: cov(*)
    !! target d*n array
    real(RK), intent(inout)       :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout)       :: w(*)
    !! work array, must be larger than Kabsch_worksize(d)
    !! if row_major, must be larger than Kabsch_worksize(n)
    integer(IK), intent(in), optional :: ldcov
    integer(IK)                   :: dd, m, s, u, vt, iw
!
    if (d < 1) RETURN
    if (d == 1)then
      rot(1) = ONE
      RETURN
    endif
!
    dd = d * d
    m = 1
    u = m + dd
    vt = u + dd
    s = vt + dd
    iw = s + d
!
    call copy(dd, cov, w(m))
    call svd(d, w(m), w(s), w(u), w(vt), w(iw), ldx=ldcov)
!
    call DGEMM('N', 'N', d, d, d, ONE, w(u), d, w(vt), d, ZERO, w(s), d)
    call det_sign(d, w(s:s + dd - 1))
    if (w(s) < ZERO) w(u + dd - d:u + dd - 1) = -w(u + dd - d:u + dd - 1)
    call DGEMM('N', 'N', d, d, d, ONE, w(u), d, w(vt), d, ZERO, w(s), d)
!
    rot(:dd) = w(s:s+dd-1)
!
  end subroutine Kabsch
!
  pure subroutine copy(d, src, dst)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: src(d)
    real(RK), intent(inout) :: dst(d)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dst(i) = src(i)
    end do
  end subroutine copy
!
end module mod_Kabsch

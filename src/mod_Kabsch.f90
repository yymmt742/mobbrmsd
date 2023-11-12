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
  interface
    include 'dgemm.h'
  end interface
!
contains
!
!| Calculate work array size for d*d matrix.
  pure elemental function Kabsch_worksize(d) result(res)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    integer(IK)                   :: res
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
  pure subroutine Kabsch(d, n, cov, rot, w)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)       :: n
    !! matrix row dimension.
    real(RK), intent(in)          :: cov(*)
    !! target d*n array
    real(RK), intent(inout)       :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout)       :: w(*)
    !! work array, must be larger than Kabsch_worksize(d)
    !! if row_major, must be larger than Kabsch_worksize(n)
    integer(IK)                   :: dd, m, s, u, vt, iw, l
!
    if (d < 1) RETURN
!
    l = MIN(d, n)
    dd = d * d
    m = 1
    u = m + dd
    vt = u + dd
    s = vt + dd
    iw = s + d
!
    w(:dd) = cov(:dd)
!
    call svd(d, w(m), w(s), w(u), w(vt), w(iw))
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
end module mod_Kabsch

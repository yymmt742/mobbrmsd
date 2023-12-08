!| Calculate the rotation matrix that minimizes |X-RY|^2 using the Kabsch-Umeyama algorithm.
!  Here, RR^T=I and det(R)=1 are satisfied.
module mod_Kabsch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_svd
  use mod_pca
  use mod_det
  implicit none
  private
  public :: Kabsch_worksize, Kabsch
  public :: get_rotation_matrix
!
  real(RK), parameter :: DEF_threshold = 1.0E-8_RK
!
  interface
    include 'dsyev.h'
    include 'dgemm.h'
  end interface
!
contains
!
!| Calculate work array size for d*d matrix.
  pure elemental function Kabsch_worksize(d, n) result(res)
    integer(IK), intent(in)           :: d
    !! matrix collumn dimension.
    integer(IK), intent(in), optional :: n
    !! matrix row dimension.
    integer(IK)                       :: res
!
    if (d < 1) then
      res = 0
    else
      if (PRESENT(n)) then
        res = MAX(d * d + 2 * MAX(d - 1, 0) * n + MAX(d - 1, 0)**2 * 3 + svd_worksize(MAX(d - 1, 0)), &
       &          d * d * 3 + 2 * d * n + d + MAX(svd_worksize(d) + d * d * 3 + d, pca_worksize(d)))
        else
        res = svd_worksize(d) + d * d * 3 + d
      end if
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
    w(:dd) = cov(:dd)
!
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
!| Calculate the rotation matrix from covariance matrix.
  pure subroutine get_rotation_matrix(d, n, X, Y, C, R, w, threshold)
    integer(IK), intent(in)        :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)        :: n
    !! matrix row dimension.
    real(RK), intent(in)           :: X(*)
    !! reference d*n array
    real(RK), intent(in)           :: Y(*)
    !! target d*n array
    real(RK), intent(in)           :: C(*)
    !! covariance matrix
    real(RK), intent(inout)        :: R(*)
    !! rotation d*d matrix
    real(RK), intent(inout)        :: W(*)
    !! work array, must be larger than Kabsch_worksize(d, n)
    real(RK), intent(in), optional :: threshold
    !! threshold
    real(RK)                       :: threshold_
    integer(IK)                    :: z, s, u, iw, ds, dp
!
    if (d < 1) RETURN
    if (d == 1)then
      R(1) = ONE
      RETURN
    endif
!
    ds = MIN(d, n)
    u = 1
    s = u + d * d
    z = s + d
    iw = z + d * n
    call pca_rotation(d, n, X, Y, w(z), w(u), w(s), w(iw))
!
    threshold_ = DEF_threshold
    if (PRESENT(threshold)) threshold_ = threshold
!
    dp = COUNT(w(s:s + ds - 1) > threshold_)
!
    if (dp < d) then
      block
        integer(IK) :: ux, yu, cov, uu, vt
        ux = z
        yu = ux + dp * n
        cov = yu + dp * n
        uu = cov + dp * dp
        vt = uu + dp * dp
        iw = vt + dp * dp
        call partial_Kabsch(d, dp, n, C, w(U), R, w(ux), w(yu), w(cov), w(uu), w(vt), w(iw))
      end block
    else
      call Kabsch(d, C, R, w)
    endif
!
  end subroutine get_rotation_matrix
!
  pure subroutine pca_rotation(d, n, X, Y, Z, U, S, W)
    integer(IK), intent(in)      :: d, n
    real(RK), intent(in)         :: X(d, n), Y(d, n)
    real(RK), intent(inout)      :: Z(d, n), U(d, d), S(d), W(*)
    integer(IK)                  :: i, j
!
      do concurrent(i=1:d, j=1:n)
        Z(i, j) = X(i, j) - Y(i, j)
      end do
!
      call pca(.FALSE., d, n, Z, U, S, W)
!
  end subroutine pca_rotation
!
  pure subroutine partial_Kabsch(d, dp, n, CC, U, R, UX, YU, C, UU, VT, W)
    integer(IK), intent(in) :: d, dp, n
    real(RK), intent(in)    :: CC(d, d), U(d, d)
    real(RK), intent(inout) :: R(d, d), UX(dp, n), YU(n, dp), C(dp, dp)
    real(RK), intent(inout) :: UU(dp,dp), VT(dp,dp), W(*)
    integer(IK)             :: i, j, dd
!
    do concurrent(i=1:d, j=1:d)
      R(i, j) = MERGE(ONE, ZERO, i == j)
    end do
!
    if (dp < 1) return
!
    dd = dp * dp
    call DGEMM('T', 'N', dp, d, d, ONE, U, d, CC, d, ZERO, W, dp)
    call DGEMM('N', 'N', dp, dp, d, ONE, W, dp,  U, d, ZERO, C, dp)
!
    if (dp == 1) then
      UU(1, 1) = SIGN(ONE, C(1, 1))
      VT(1, 1) = ONE
    else
      call svd(dp, C, w, UU, VT, w(dp+1))
      call copy(dd, UU, w)
      call det_sign(dp, w)
      call copy(dd, VT, w(2))
      call det_sign(dp, w(2:dd+1))
      if (w(1) * w(2) < ZERO) UU(:, dp) = -UU(:, dp)
    end if
!
    call DGEMM('N', 'N', dp, dp, dp, ONE, UU, dp, VT, dp, ZERO, R, d)
    call DGEMM('N', 'N', d, d, d, ONE, U, d, R, d, ZERO, W, d)
    call DGEMM('N', 'T', d, d, d, ONE, W, d, U, d, ZERO, R, d)
!
  end subroutine partial_Kabsch
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

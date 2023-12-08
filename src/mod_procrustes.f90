!| Calculate the rotation matrix that minimizes |X-RY|^2 using the procrustes-Umeyama algorithm.
!  Here, RR^T=I is satisfied.
module mod_procrustes
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_svd
  use mod_pca
  implicit none
  private
  public :: procrustes_worksize, procrustes
  public :: get_permutation_matrix
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
  pure elemental function procrustes_worksize(d) result(res)
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
  end function procrustes_worksize
!
!| Calculate the rotation matrix from covariance matrix.
  pure subroutine procrustes(d, cov, rot, w, ldcov)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    real(RK), intent(in)          :: cov(*)
    !! target d*n array
    real(RK), intent(inout)       :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout)       :: w(*)
    !! work array, must be larger than procrustes_worksize(d)
    !! if row_major, must be larger than procrustes_worksize(n)
    integer(IK), intent(in), optional :: ldcov
    integer(IK)                   :: dd, m, s, u, vt, iw
!
    if (d < 1) RETURN
    if (d == 1)then
      rot(1) = SIGN(ONE, cov(1))
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
    call DGEMM('N', 'N', d, d, d, ONE, w(u), d, w(vt), d, ZERO, w(s), d)
    rot(:dd) = w(s:s+dd-1)
!
  end subroutine procrustes
!
!| Calculate the rotation matrix from covariance matrix.
  pure subroutine get_permutation_matrix(d, n, X, Y, C, R, w, threshold)
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
    !! work array, must be larger than procrustes_worksize(d, n)
    real(RK), intent(in), optional :: threshold
    !! threshold
    real(RK)                       :: threshold_
    integer(IK)                    :: s, u, iw, dn, np
!
    if (d < 1) RETURN
    if (n == 1)then
      R(1) = SIGN(ONE, C(1))
      RETURN
    endif
!
    dn = MIN(d, n)
    u = 1
    s = u + n * n
    iw = s + n
    call pca_rotation(d, n, d * n, X, Y, w(u), w(s), w(iw))
!
    threshold_ = DEF_threshold
    if (PRESENT(threshold)) threshold_ = threshold
!
    np = COUNT(w(s:s + dn - 1) > threshold_)
!
    if (np < n) then
      block
        integer(IK) :: ux, yu, cov, uu, vt
        ux = s
        yu = ux + d * np
        cov = yu + d * np
        uu = cov + np * np
        vt = uu + np * np
        iw = vt + np * np
        call partial_procrustes(d, n, np, C, w(u), R, w(ux), w(yu), w(cov), w(uu), w(vt), w(iw))
      end block
    else
      call procrustes(d, C, R, w)
    endif
!
  end subroutine get_permutation_matrix
!
  pure subroutine pca_rotation(d, n, dn, X, Y, U, S, W)
    integer(IK), intent(in)      :: d, n, dn
    real(RK), intent(in)         :: X(dn), Y(dn)
    real(RK), intent(inout)      :: U(n, n), S(n), W(*)
    integer(IK)                  :: i
!
      do concurrent(i=1:dn)
        W(i) = X(i) - Y(i)
      end do
!
      call pca(.true., d, n, W, U, S, W(dn + 1))
!
  end subroutine pca_rotation
!
  pure subroutine partial_procrustes(d, n, np, CC, U, R, UX, YU, C, UU, VT, W)
    integer(IK), intent(in) :: d, n, np
    real(RK), intent(in)    :: CC(n, n), U(n, n)
    real(RK), intent(inout) :: R(n, n), UX(d, np), YU(np, d), C(np, np)
    real(RK), intent(inout) :: UU(np, np), VT(np, np), W(*)
    integer(IK)             :: i, j, nn
!
    do concurrent(i=1:n, j=1:n)
      R(i, j) = MERGE(ONE, ZERO, i == j)
    end do
!
    if (np < 1) return
!
    nn = np * np
    call DGEMM('T', 'N', np, n, n, ONE, U, n, CC, n, ZERO, W, np)
    call DGEMM('N', 'N', np, np, n, ONE, W, np,  U, n, ZERO, C, np)
!
    if (np == 1) then
      UU(1, 1) = SIGN(ONE, C(1, 1))
      VT(1, 1) = ONE
    else
      call svd(np, C, W, UU, VT, W(np+1))
    end if
!
    call DGEMM('N', 'N', np, np, np, ONE, UU, np, VT, np, ZERO, R, n)
    call DGEMM('N', 'N', n, n, n, ONE, U, n, R, n, ZERO, W, n)
    call DGEMM('N', 'T', n, n, n, ONE, W, n, U, n, ZERO, R, n)
!
  end subroutine partial_procrustes
!
end module mod_procrustes

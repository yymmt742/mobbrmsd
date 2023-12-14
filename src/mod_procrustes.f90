!| Calculate the rotation matrix that minimizes |X-RY|^2 using the procrustes-Umeyama algorithm.
!  Here, RR^T=I is satisfied.
module mod_procrustes
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_optarg
  use mod_svd
  implicit none
  private
  public :: procrustes_worksize, procrustes
! public :: get_permutation_matrix
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
    call copy(d, optarg(ldcov, d), cov, w(m))
    call svd(d, w(m), w(s), w(u), w(vt), w(iw))
    call DGEMM('N', 'N', d, d, d, ONE, w(u), d, w(vt), d, ZERO, w(s), d)
    rot(:dd) = w(s:s+dd-1)
!
  end subroutine procrustes
!
  pure subroutine copy(d, lx, X, Z)
    integer(IK), intent(in) :: d, lx
    real(RK), intent(in)    :: X(lx, *)
    real(RK), intent(inout) :: Z(d, *)
    integer(IK)             :: i, j
    do concurrent(i=1:d, j=1:d)
      Z(i, j) = X(i, j)
    end do
  end subroutine copy
!
end module mod_procrustes

!| Calculate the rotation matrix that minimizes |X-RY|^2 using the Kabsch-Umeyama algorithm.
!  Here, RR^T=I and det(R)=1 are satisfied.
module mod_rotation
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, HALF => RHALF
  use mod_params, only: gemm, gesvd
  use mod_det
  implicit none
  private
  public :: sdmin_worksize
  public :: estimate_sdmin
  public :: rotation_worksize
  public :: estimate_rotation
!
  real(RK), parameter    :: THRESHOLD = 1E-8_RK
  integer(IK), parameter :: MAXITER = 100000
!
contains
!
!| Inquire function for memory size of estimate_sdmin.
  pure elemental function sdmin_worksize() result(res)
    integer(IK) :: res
    if (D == 1) then
      res = 1
    else
      res = worksize_Kabsch() + DD + 1
    endif
  end function sdmin_worksize
!
!| Compute the least-squares sum_i^n |x_i-Ry_i|^2 from cov = YX^T and g = tr[XX^T] + tr[YY^T].
  pure subroutine estimate_sdmin(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_sdmin().
!
    if (D == 1) then
      W(1) = g - cov(1) - cov(1)
    else
      call Kabsch(cov, w(2), w(d * d + 2))
      w(1) = dot(DD, cov, w(2))
      w(1) = w(1) + w(1)
      w(1) = g - w(1)
    end if
!
  end subroutine estimate_sdmin
!
!| Inquire function for memory size of rotation.
  pure elemental function rotation_worksize() result(res)
    integer(IK) :: res
    if (d == 1) then
      res = 0
    else
      res = worksize_Kabsch()
    endif
  end function rotation_worksize
!
  !| work array size for Kabsch algorithm.
  pure elemental function worksize_Kabsch() result(res)
    real(RK)    :: w(1)
    integer(IK) :: res, info
!
    call gesvd('A', 'A', D, D, w, D, w, w, D, w, D, w, -1, info)
    res = NINT(w(1)) + DD * 3 + D
!
  end function worksize_Kabsch
!
!| Compute the transpose rotation matrix for minimize tr[CR] from cov = YX^T and g = tr[XX^T] + tr[YY^T].
  pure subroutine estimate_rotation(g, cov, rot, w)
    real(RK), intent(in)    :: g
    !! g = tr[XX^T] + tr[YY^T]
    real(RK), intent(in)    :: cov(*)
    !! covariance dxd matrix, YX^T
    real(RK), intent(inout) :: rot(*)
    !! rotation dxd matrix
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_rotation().
!
    if (d == 1) then
      rot(1) = ONE
    else
      call Kabsch(cov, rot, w)
    end if
!
  end subroutine estimate_rotation
!
!| Calculate the rotation matrix R^T from covariance matrix.
  pure subroutine Kabsch(cov, rot, w)
    real(RK), intent(in)          :: cov(*)
    !! target d*n array
    real(RK), intent(inout)       :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout)       :: w(*)
    !! work array, must be larger than Kabsch_worksize(d)
    !! if row_major, must be larger than Kabsch_worksize(n)
    integer(IK), parameter        :: m = 1
    integer(IK)                   :: s, u, vt, iw, lw, info
!
    u = m + dd
    vt = u + dd
    s = vt + dd
    iw = s + d
!
    call gesvd('A', 'A', d, d, w(m), d, w(s), w(u), d, w(vt), d, w(iw), -1, info)
    lw = NINT(w(iw))
!
    call copy(dd, cov, w(m))
    call gesvd('A', 'A', d, d, w(m), d, w(s), w(u), d, w(vt), d, w(iw), lw, info)
!
    call gemm('N', 'N', d, d, d, ONE, w(u), d, w(vt), d, ZERO, w(s), d)
    call det_sign(w(s:s + dd - 1))
    if (w(s) < ZERO) w(u + dd - d:u + dd - 1) = -w(u + dd - d:u + dd - 1)
    call gemm('N', 'N', d, d, d, ONE, w(u), d, w(vt), d, ZERO, w(s), d)
    call copy(dd, w(s), rot(1))
!
  end subroutine Kabsch
!
  pure function dot(N, X, Y) result(res)
    integer(IK), intent(in) :: N
    real(RK), intent(in)    :: X(*), Y(*)
    real(RK)                :: res
    integer(IK)             :: i
    res = ZERO
    do i = 1, N
      res = res + X(i) * Y(i)
    end do
  end function dot
!
  pure subroutine copy(N, X, Y)
    integer(IK), intent(in) :: N
    real(RK), intent(in)    :: X(*)
    real(RK), intent(inout) :: Y(*)
    integer(IK)             :: i
    do concurrent(i=1:N)
      Y(i) = X(i)
    end do
  end subroutine copy
!
end module mod_rotation


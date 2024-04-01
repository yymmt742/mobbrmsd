!| Calculate the rotation matrix that minimizes |X-RY|^2 for D=2. <br>
!  Here, RR^T=I and det(R)=1 are satisfied. <br>
!  This code is based on the method of Theobald. doi : 10.1107/S0108767305015266
module mod_rotation
  use mod_kinds, only: IK, RK
  implicit none
  private
  public :: sdmin_worksize
  public :: estimate_sdmin
  public :: rotation_worksize
  public :: estimate_rotation
!
  real(RK), parameter    :: ZERO = 0.0_RK
  real(RK), parameter    :: HALF = 0.5_RK
  real(RK), parameter    :: ONE = 1.0_RK
  real(RK), parameter    :: THRESHOLD = 1E-14_RK
!
contains
!
!| Inquire function for memory size of estimate_sdmin.
  pure elemental function sdmin_worksize() result(res)
    integer(IK) :: res
    res = 2
  end function sdmin_worksize
!
!| Compute the least-squares sum_i^n |x_i-Ry_i|^2 from cov = YX^T and g = tr[XX^T] + tr[YY^T].
  pure subroutine estimate_sdmin(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*d array
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_sdmin().
!
    if (g < THRESHOLD) then
      w(1) = ZERO
    else
      w(1) = cov(1) + cov(4)
      w(1) = w(1) * w(1)
      w(2) = cov(2) - cov(3)
      w(2) = w(2) * w(2)
      w(1) = SQRT(w(1) + w(2))
      w(1) = w(1) + w(1)
      w(1) = g - w(1)
    end if
!
  end subroutine estimate_sdmin
!
!| Inquire function for memory size of rotation.
  pure elemental function rotation_worksize() result(res)
    integer(IK) :: res
    res = 0
  end function rotation_worksize
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
    if (g < THRESHOLD) then
      rot(1) = ONE; rot(2) = ZERO
      rot(3) = ZERO; rot(4) = ONE
      return
    end if
!
    rot(1) = cov(1) + cov(4)
    rot(2) = cov(3) - cov(2)
    rot(3) = rot(1) * rot(1) + rot(2) * rot(2)
!
    if (rot(3) < THRESHOLD) then
      rot(1) = ONE;  rot(2) = ZERO
      rot(3) = ZERO; rot(4) = ONE
      return
    end if
!
    rot(3) = SQRT(ONE / rot(3))
    rot(1) = rot(1) * rot(3)      ! cos theta
    rot(3) = HALF + HALF * rot(1) ! cos^2 theta/2
    rot(4) = rot(3) - ONE         ! - sin^2 theta/2
    rot(1) = rot(4) + rot(3)      ! cos^2 theta/2 - sin^2 theta/2
    rot(4) = rot(4) * rot(3)      ! cos^2 theta/2 * sin^2 theta/2
    rot(4) = rot(4) + rot(4)      ! 2 * cos^2 theta/2 * sin^2 theta/2
    rot(4) = rot(4) + rot(4)      ! 4 * cos^2 theta/2 * sin^2 theta/2
    rot(2) = rot(4) * rot(2)
    rot(2) = SIGN(ONE, rot(2)) * SQRT(ABS(rot(4)))
    rot(3) = -rot(2)
    rot(4) = rot(1)
!
  end subroutine estimate_rotation
!
end module mod_rotation


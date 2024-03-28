!| Calculate the rotation matrix that minimizes |X-RY|^2 for D=3. <br>
!  Here, RR^T=I and det(R)=1 are satisfied.
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
  integer(IK), parameter :: MAXITER = 100000
!
contains
!
!| Inquire function for memory size of estimate_sdmin.
  pure elemental function sdmin_worksize() result(res)
    integer(IK) :: res
    res = 10
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
    integer(IK), parameter  :: k11 = 1, k22 = 2, k33 = 3, k44 = 4
    integer(IK), parameter  :: k21 = 5, k31 = 6, k41 = 7
    integer(IK), parameter  :: k32 = 8, k42 = 9, k43 = 10
    integer(IK), parameter  :: a1 = 5, a2 = 6, b2 = 7, b8 = 8, c1 = 9
    integer(IK), parameter  :: l1 = 1, l0 = 2, l2 = 3, l3 = 4
    integer(IK)             :: i
!
    if (g < THRESHOLD) then
      w(1) = ZERO
      return
    end if
!
    w(k11) = cov(1) + cov(5) + cov(9)
    w(k21) = cov(8) - cov(6)
    w(k31) = cov(3) - cov(7)
    w(k41) = cov(4) - cov(2)
    w(k22) = cov(1) - cov(5) - cov(9)
    w(k32) = cov(4) + cov(2)
    w(k42) = cov(3) + cov(7)
    w(k33) = -cov(1) + cov(5) - cov(9)
    w(k43) = cov(8) + cov(6)
    w(k44) = -cov(1) - cov(5) + cov(9)
!
    w(c1) = (w(k11) * w(k22) - w(k21) * w(k21)) * (w(k33) * w(k44) - w(k43) * w(k43)) &
   &      + (w(k21) * w(k31) - w(k11) * w(k32)) * (w(k32) * w(k44) - w(k43) * w(k42)) &
   &      + (w(k11) * w(k42) - w(k21) * w(k41)) * (w(k32) * w(k43) - w(k33) * w(k42)) &
   &      + (w(k21) * w(k32) - w(k22) * w(k31)) * (w(k31) * w(k44) - w(k43) * w(k41)) &
   &      + (w(k22) * w(k41) - w(k21) * w(k42)) * (w(k31) * w(k43) - w(k41) * w(k33)) &
   &      + (w(k31) * w(k42) - w(k41) * w(k32))**2
!
    w(a1) = cov(1) * cov(1) + cov(2) * cov(2) + cov(3) * cov(3) &
   &      + cov(4) * cov(4) + cov(5) * cov(5) + cov(6) * cov(6) &
   &      + cov(7) * cov(7) + cov(8) * cov(8) + cov(9) * cov(9)
    w(a2) = w(a1) + w(a1)
!
    w(b2) = cov(1) * (cov(8) * cov(6) - cov(5) * cov(9)) &
   &      + cov(4) * (cov(2) * cov(9) - cov(8) * cov(3)) &
   &      + cov(7) * (cov(5) * cov(3) - cov(2) * cov(6))
!
    w(b2) = w(b2) + w(b2)
    w(b8) = w(b2) + w(b2)
    w(b8) = w(b8) + w(b8)
!
    if (ABS(w(b8)) < THRESHOLD) then
      w(l2) = w(c1) + w(c1)
      w(l2) = w(l2) + w(l2)
      w(l1) = SQRT(ABS(HALF * (SQRT(ABS(w(a2) * w(a2) - w(l2))) + w(a2))))
    else
      w(l1) = g
      do i = 1, MAXITER
        w(l2) = w(l1) * w(l1)
        w(l3) = w(l2) * w(l1)
        w(l0) = w(l2) * w(l2)
        w(l2) = w(l0) - w(a2) * w(l2) + w(b8) * w(l1) + w(c1)
        w(l3) = w(l3) - w(a1) * w(l1) + w(b2)
        w(l3) = w(l3) + w(l3)
        w(l3) = w(l3) + w(l3)
        w(l0) = w(l1)
        w(l1) = w(l1) - w(l2) / w(l3)
        if (ABS(w(l0) - w(l1)) < THRESHOLD) exit
      end do
    end if
!
    w(l1) = w(l1) + w(l1)
    w(l1) = g - w(l1)
!
  end subroutine estimate_sdmin
!
!| Inquire function for memory size of rotation.
  pure elemental function rotation_worksize() result(res)
    integer(IK) :: res
    res = 28
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
    integer(IK), parameter  :: k11 = 1, k22 = 2, k33 = 3, k44 = 4
    integer(IK), parameter  :: k21 = 5, k31 = 6, k41 = 7
    integer(IK), parameter  :: k32 = 8, k42 = 9, k43 = 10
    integer(IK), parameter  :: a1 = 11, a2 = 12, b2 = 13, b8 = 14, c1 = 15
    integer(IK), parameter  :: l0 = 16, l1 = 17, l2 = 18, l3 = 19
    integer(IK), parameter  :: q0 = 11, q1 = 12, q2 = 13, q3 = 14, q4 = 15
!
    if (g < THRESHOLD) then
      rot(1) =  ONE; rot(2) = ZERO; rot(3) = ZERO
      rot(4) = ZERO; rot(5) =  ONE; rot(6) = ZERO
      rot(7) = ZERO; rot(8) = ZERO; rot(9) =  ONE
      return
    end if
    w(9) = (ONE + ONE) / g
    w(1) = w(9) * cov(1); w(2) = w(9) * cov(2); w(3) = w(9) * cov(3)
    w(4) = w(9) * cov(4); w(5) = w(9) * cov(5); w(6) = w(9) * cov(6)
    w(7) = w(9) * cov(7); w(8) = w(9) * cov(8); w(9) = w(9) * cov(9)
!
    call quartenion_rotmatrix_d3(g, w(1), rot, w(10))
!
  end subroutine estimate_rotation
!
  pure subroutine quartenion_rotmatrix_d3(g, s, r, w)
    real(RK), intent(in)    :: g, s(*)
    !! target d*n array
    real(RK), intent(inout) :: r(*)
    !! rotation d*d matrix
    real(RK), intent(inout) :: w(*)
    integer(IK), parameter  :: k11 = 1, k22 = 2, k33 = 3, k44 = 4
    integer(IK), parameter  :: k21 = 5, k31 = 6, k41 = 7
    integer(IK), parameter  :: k32 = 8, k42 = 9, k43 = 10
    integer(IK), parameter  :: a1 = 11, a2 = 12, b2 = 13, b8 = 14, c1 = 15
    integer(IK), parameter  :: l0 = 16, l1 = 17, l2 = 18, l3 = 19
    integer(IK), parameter  :: q0 = 11, q1 = 12, q2 = 13, q3 = 14, q4 = 15
!
    w(a1) = s(1) * s(1) + s(2) * s(2) + s(3) * s(3) &
   &      + s(4) * s(4) + s(5) * s(5) + s(6) * s(6) &
   &      + s(7) * s(7) + s(8) * s(8) + s(9) * s(9)
    w(a2) = w(a1) + w(a1)
!
    w(b2) = s(1) * (s(8) * s(6) - s(5) * s(9)) &
   &      + s(4) * (s(2) * s(9) - s(8) * s(3)) &
   &      + s(7) * (s(5) * s(3) - s(2) * s(6))
    w(b2) = w(b2) + w(b2)
    w(b8) = w(b2) + w(b2)
    w(b8) = w(b8) + w(b8)
!
    w(k11) =  s(1) + s(5) + s(9)
    w(k21) =  s(8) - s(6)
    w(k31) =  s(3) - s(7)
    w(k41) =  s(4) - s(2)
    w(k22) =  s(1) - s(5) - s(9)
    w(k32) =  s(4) + s(2)
    w(k42) =  s(3) + s(7)
    w(k33) = -s(1) + s(5) - s(9)
    w(k43) =  s(8) + s(6)
    w(k44) = -s(1) - s(5) + s(9)
!
    w(c1) = (w(k11) * w(k22) - w(k21) * w(k21)) * (w(k33) * w(k44) - w(k43) * w(k43)) &
   &      + (w(k21) * w(k31) - w(k11) * w(k32)) * (w(k32) * w(k44) - w(k43) * w(k42)) &
   &      + (w(k11) * w(k42) - w(k21) * w(k41)) * (w(k32) * w(k43) - w(k33) * w(k42)) &
   &      + (w(k21) * w(k32) - w(k22) * w(k31)) * (w(k31) * w(k44) - w(k43) * w(k41)) &
   &      + (w(k22) * w(k41) - w(k21) * w(k42)) * (w(k31) * w(k43) - w(k41) * w(k33)) &
   &      + (w(k31) * w(k42) - w(k41) * w(k32)) **2
!
    if (ABS(w(b8)) < THRESHOLD) then
      w(l2) = w(c1) + w(c1)
      w(l2) = w(l2) + w(l2)
      w(l1) = SQRT(ABS(HALF * (SQRT(ABS(w(a2) * w(a2) - w(l2))) + w(a2))))
    else
      w(l1) = ONE
      do
        w(l2) = w(l1) * w(l1)
        w(l3) = w(l2) * w(l1)
        w(l0) = w(l2) * w(l2)
        w(l2) = w(l0) - w(a2) * w(l2) + w(b8) * w(l1) + w(c1)
        w(l3) = w(l3) - w(a1) * w(l1) + w(b2)
        w(l3) = w(l3) + w(l3)
        w(l3) = w(l3) + w(l3)
        w(l0) = w(l1)
        w(l1) = w(l1) - w(l2) / w(l3)
        if (ABS(w(l0) - w(l1)) < THRESHOLD) exit
      end do
    end if
!
    w(k11) = w(k11) - w(l1)
    w(k22) = w(k22) - w(l1)
    w(k33) = w(k33) - w(l1)
    w(k44) = w(k44) - w(l1)
!
    w(q1) = w(k22) * (w(k33) * w(k44) - w(k43) * w(k43)) &
   &      + w(k32) * (w(k43) * w(k42) - w(k32) * w(k44)) &
   &      + w(k42) * (w(k32) * w(k43) - w(k33) * w(k42))
    w(q2) = w(k21) * (w(k43) * w(k43) - w(k33) * w(k44)) &
   &      + w(k31) * (w(k32) * w(k44) - w(k43) * w(k42)) &
   &      + w(k41) * (w(k33) * w(k42) - w(k32) * w(k43))
    w(q3) = w(k21) * (w(k32) * w(k44) - w(k43) * w(k42)) &
   &      + w(k22) * (w(k43) * w(k41) - w(k31) * w(k44)) &
   &      + w(k42) * (w(k31) * w(k42) - w(k32) * w(k41))
    w(q4) = w(k21) * (w(k33) * w(k42) - w(k32) * w(k43)) &
   &      + w(k22) * (w(k31) * w(k43) - w(k33) * w(k41)) &
   &      + w(k32) * (w(k32) * w(k41) - w(k31) * w(k42))
!
    w(k11) = w(q1) * w(q1)
    w(k22) = w(q2) * w(q2)
    w(k33) = w(q3) * w(q3)
    w(k44) = w(q4) * w(q4)
    w(q0) = ONE / (w(k11) + w(k22) + w(k33) + w(k44))
    w(k11) = w(k11) * w(q0)
    w(k22) = w(k22) * w(q0)
    w(k33) = w(k33) * w(q0)
    w(k44) = w(k44) * w(q0)
    w(q0) = w(q0) + w(q0)
    w(k21) = w(q1) * w(q2) * w(q0)
    w(k31) = w(q1) * w(q3) * w(q0)
    w(k41) = w(q1) * w(q4) * w(q0)
    w(k32) = w(q2) * w(q3) * w(q0)
    w(k42) = w(q2) * w(q4) * w(q0)
    w(k43) = w(q3) * w(q4) * w(q0)
!
    r(1) = w(k11) + w(k22) - w(k33) - w(k44)
    r(2) = w(k32) - w(k41)
    r(3) = w(k42) + w(k31)
    r(4) = w(k32) + w(k41)
    r(5) = w(k11) - w(k22) + w(k33) - w(k44)
    r(6) = w(k43) - w(k21)
    r(7) = w(k42) - w(k31)
    r(8) = w(k43) + w(k21)
    r(9) = w(k11) - w(k22) - w(k33) + w(k44)
!
  end subroutine quartenion_rotmatrix_d3
!
end module mod_rotation


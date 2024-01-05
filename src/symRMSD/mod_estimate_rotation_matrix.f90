!| Calculate the rotation matrix that minimizes |X-RY|^2 using the Kabsch-Umeyama algorithm.
!  Here, RR^T=I and det(R)=1 are satisfied.
module mod_estimate_rotation_matrix
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, HALF => RHALF
  use mod_det
  implicit none
  private
  public :: estimate_sdmin
  public :: estimate_rotation_matrix
!
  interface
    include 'ddot.h'
    include 'dcopy.h'
    include 'dgemm.h'
    include 'dgesvd.h'
  end interface
!
  real(RK), parameter :: THRESHOLD = 1E-8_RK
!
contains
!
  pure subroutine estimate_sdmin(d, g, cov, w)
    integer(IK), intent(in) :: d
    !! spatial dimension
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize(d).
!
    if (d == 0) then
      W(1) = 1
    elseif (d == -1) then
      W(1) = 1
    elseif (d == -2) then
      W(1) = 2
    elseif (d == -3) then
      W(1) = 10
    elseif (d == 1) then
      W(1) = g - cov(1) - cov(1)
    elseif (d == 2) then
      if (g < THRESHOLD) then
        w(1) = ZERO
        return
      end if
      call quartenion_sdmin_d2(g, cov, w)
    elseif (d == 3) then
      if (g < THRESHOLD) then
        w(1) = ZERO
        return
      end if
      call quartenion_sdmin_d3(g, cov, w)
    elseif (d < -3) then
      call Kabsch(d, cov, w, w)
      w(1) = w(1) + d * d + 1
    else
      call Kabsch(d, cov, w(2), w(d * d + 2))
      w(1) = ddot(d * d, cov, 1, w(2), 1)
      w(1) = w(1) + w(1)
      w(1) = g - w(1)
    end if
!
  end subroutine estimate_sdmin
!
  pure subroutine quartenion_sdmin_d2(g, c, w)
    real(RK), intent(in)    :: g, c(*)
    real(RK), intent(inout) :: w(*)
!
    w(1) = c(1) + c(4)
    w(1) = w(1) * w(1)
    w(2) = c(2) - c(3)
    w(2) = w(2) * w(2)
    w(1) = SQRT(w(1) + w(2))
    w(1) = w(1) + w(1)
    w(1) = g - w(1)
!
  end subroutine quartenion_sdmin_d2
!
  pure subroutine quartenion_sdmin_d3(g, c, w)
    real(RK), intent(in)    :: g, c(*)
    real(RK), intent(inout) :: w(*)
    integer(IK), parameter  :: k11 = 1, k22 = 2, k33 = 3, k44 = 4
    integer(IK), parameter  :: k21 = 5, k31 = 6, k41 = 7
    integer(IK), parameter  :: k32 = 8, k42 = 9, k43 = 10
    integer(IK), parameter  :: a1 = 5, a2 = 6, b2 = 7, b8 = 8, c1 = 9
    integer(IK), parameter  :: l1 = 1, l0 = 2, l2 = 3, l3 = 4
!
    w(k11) =  c(1) + c(5) + c(9)
    w(k21) =  c(8) - c(6)
    w(k31) =  c(3) - c(7)
    w(k41) =  c(4) - c(2)
    w(k22) =  c(1) - c(5) - c(9)
    w(k32) =  c(4) + c(2)
    w(k42) =  c(3) + c(7)
    w(k33) = -c(1) + c(5) - c(9)
    w(k43) =  c(8) + c(6)
    w(k44) = -c(1) - c(5) + c(9)
!
    w(c1) = (w(k11) * w(k22) - w(k21) * w(k21)) * (w(k33) * w(k44) - w(k43) * w(k43)) &
   &      + (w(k21) * w(k31) - w(k11) * w(k32)) * (w(k32) * w(k44) - w(k43) * w(k42)) &
   &      + (w(k11) * w(k42) - w(k21) * w(k41)) * (w(k32) * w(k43) - w(k33) * w(k42)) &
   &      + (w(k21) * w(k32) - w(k22) * w(k31)) * (w(k31) * w(k44) - w(k43) * w(k41)) &
   &      + (w(k22) * w(k41) - w(k21) * w(k42)) * (w(k31) * w(k43) - w(k41) * w(k33)) &
   &      + (w(k31) * w(k42) - w(k41) * w(k32)) **2
!
    w(a1) = c(1) * c(1) + c(2) * c(2) + c(3) * c(3) &
   &      + c(4) * c(4) + c(5) * c(5) + c(6) * c(6) &
   &      + c(7) * c(7) + c(8) * c(8) + c(9) * c(9)
    w(a2) = w(a1) + w(a1)
!
    w(b2) = c(1) * (c(8) * c(6) - c(5) * c(9)) &
   &      + c(4) * (c(2) * c(9) - c(8) * c(3)) &
   &      + c(7) * (c(5) * c(3) - c(2) * c(6))
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
    w(l1) = w(l1) + w(l1)
    w(l1) = g - w(l1)
!
  end subroutine quartenion_sdmin_d3
!
  pure subroutine estimate_rotation_matrix(d, g, cov, rot, w)
    integer(IK), intent(in) :: d
    !! spatial dimension
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than Kabsch_worksize(d)
    !! if row_major, must be larger than Kabsch_worksize(n)
!
    if (d == 0) then
      W(1) = 0
    elseif (d == -1) then
      W(1) = 0
    elseif (d == -2) then
      W(1) = 0
    elseif (d == -3) then
      W(1) = 28
    elseif (d == 1) then
      rot(1) = ONE
    elseif (d == 2) then
      if (g < THRESHOLD) then
        rot(1) =  ONE; rot(2) = ZERO
        rot(3) = ZERO; rot(4) = ONE
        return
      end if
      call quartenion_rotmatrix_d2(cov, rot)
    elseif (d == 3) then
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
      call quartenion_rotmatrix_d3(g, w(1), rot, w(10))
    else
      call Kabsch(d, cov, rot, w)
    end if
!
  end subroutine estimate_rotation_matrix
!
  pure subroutine quartenion_rotmatrix_d2(c, rot)
    real(RK), intent(in)    :: c(*)
    real(RK), intent(inout) :: rot(*)
!
    rot(1) = c(1) + c(4)
    rot(2) = c(3) - c(2)
    rot(3) = rot(1) * rot(1) + rot(2) * rot(2)
!
    if (rot(3) < Threshold) then
      rot(1) = ONE;  rot(2) = ZERO
      rot(3) = ZERO; rot(4) = ONE
      return
    end if
!
    rot(3) = SQRT(ONE / rot(3))
    rot(1) = rot(1) * rot(3)       ! cos theta
    rot(3) = HALF + HALF * rot(1)  ! cos^2 theta/2
    rot(4) = rot(3) - ONE          ! - sin^2 theta/2
    rot(1) = rot(4) + rot(3)       ! cos^2 theta/2 - sin^2 theta/2
    rot(4) = rot(4) * rot(3)       ! cos^2 theta/2 * sin^2 theta/2
    rot(4) = rot(4) + rot(4)       ! 2 * cos^2 theta/2 * sin^2 theta/2
    rot(4) = rot(4) + rot(4)       ! 4 * cos^2 theta/2 * sin^2 theta/2
    rot(2) = rot(4) * rot(2)
    rot(2) = SIGN(ONE, rot(2)) * SQRT(ABS(rot(4)))
    rot(3) = -rot(2)
    rot(4) = rot(1)
!
  end subroutine quartenion_rotmatrix_d2
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
!| Calculate the rotation matrix R^T from covariance matrix.
  pure subroutine Kabsch(d, cov, rot, w)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    real(RK), intent(in)          :: cov(*)
    !! target d*n array
    real(RK), intent(inout)       :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout)       :: w(*)
    !! work array, must be larger than Kabsch_worksize(d)
    !! if row_major, must be larger than Kabsch_worksize(n)
    integer(IK), parameter        :: m = 1
    integer(IK)                   :: dd, s, u, vt, iw, lw, info
!
    if (d == 0) then
      w(1) = ZERO
      return
    elseif (d == 1)then
      rot(1) = ONE
      RETURN
    elseif (d < 0) then
      s = ABS(d)
      call DGESVD('A', 'A', s, s, w, s, w, w, s, w, s, w, -1, info)
      w(1) = w(1) + s * s * 3 + s
      return
    end if
!
    dd = d * d
    u = m + dd
    vt = u + dd
    s = vt + dd
    iw = s + d
!
    call dgesvd('A', 'A', d, d, w(m), d, w(s), w(u), d, w(vt), d, w(iw), -1, info)
    lw = NINT(w(iw))
!
    call dcopy(dd, cov, 1, w(m), 1)
    call dgesvd('A', 'A', d, d, w(m), d, w(s), w(u), d, w(vt), d, w(iw), lw, info)
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
end module mod_estimate_rotation_matrix

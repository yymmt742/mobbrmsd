!| Calculate the rotation matrix that minimizes |X-RY|^2 using the Kabsch-Umeyama algorithm.
!  Here, RR^T=I and det(R)=1 are satisfied.
module mod_Kabsch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, HALF => RHALF
  use mod_det
  implicit none
  private
  public :: Kabsch_worksize, Kabsch, quartenion_operation
!
  interface
    include 'dcopy.h'
    include 'dgemm.h'
    include 'dgesvd.h'
  end interface
!
  real(RK), parameter :: THRESHOLD = 1D-8
contains
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
    if(d==3)then
    else
    endif
!
  end subroutine estimate_rotation_matrix
!
!| Calculate work array size for d*d matrix.
  pure elemental function Kabsch_worksize(d) result(res)
    integer(IK), intent(in)           :: d
    !! matrix collumn dimension.
    real(RK)                          :: dum(1)
    integer(IK)                       :: res, info
!
    if (d < 1) then
      res = 0
    else
      call DGESVD('A', 'A', d, d, dum, d, dum, dum, d, dum, d, dum, -1, info)
      res =  NINT(dum(1)) + d * d * 3 + d
    end if
!
  end function Kabsch_worksize
!
!| Calculate the rotation matrix from covariance matrix.
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
    integer(IK)                   :: dd, m, s, u, vt, iw, lw, info
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
  pure subroutine quartenion_operation(g, cov, rot, w)
    real(RK), intent(in)    :: g
    !! target d*n array
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than Kabsch_worksize(d)
    !! if row_major, must be larger than Kabsch_worksize(n)
    integer(IK), parameter  :: k11 = 1, k22 = 2, k33 = 3, k44 = 4
    integer(IK), parameter  :: k21 = 5, k31 = 6, k41 = 7
    integer(IK), parameter  :: k32 = 8, k42 = 9, k43 = 10
    integer(IK), parameter  :: a1 = 11, a2 = 12, b2 = 13, b8 = 14, c1 = 15
    integer(IK), parameter  :: s1 = 16, s2 = 17, s3 = 18
    integer(IK), parameter  :: s4 = 19, s5 = 20, s6 = 21
    integer(IK), parameter  :: s7 = 22, s8 = 23, s9 = 24
    integer(IK), parameter  :: l0 = 16, l1 = 17, l2 = 18
    integer(IK), parameter  :: l3 = 19, l4 = 20, lf = 21, lg = 22
    integer(IK), parameter  :: q0 = 11, q1 = 12, q2 = 13, q3 = 14, q4 = 15
!
    w(s9) = (ONE + ONE) / g
    w(s1) = w(s9) * cov(1)
    w(s2) = w(s9) * cov(2)
    w(s3) = w(s9) * cov(3)
    w(s4) = w(s9) * cov(4)
    w(s5) = w(s9) * cov(5)
    w(s6) = w(s9) * cov(6)
    w(s7) = w(s9) * cov(7)
    w(s8) = w(s9) * cov(8)
    w(s9) = w(s9) * cov(9)
!
    w(a1) = w(s1) * w(s1) + w(s2) * w(s2) + w(s3) * w(s3) &
   &      + w(s4) * w(s4) + w(s5) * w(s5) + w(s6) * w(s6) &
   &      + w(s7) * w(s7) + w(s8) * w(s8) + w(s9) * w(s9)
    w(a2) = w(a1) + w(a1)
!
    w(b2) = w(s1) * (w(s8) * w(s6) - w(s5) * w(s9)) &
   &      + w(s4) * (w(s2) * w(s9) - w(s8) * w(s3)) &
   &      + w(s7) * (w(s5) * w(s3) - w(s2) * w(s6))
    w(b2) = w(b2) + w(b2)
    w(b8) = w(b2) + w(b2)
    w(b8) = w(b8) + w(b8)
!
    w(k11) =  w(s1) + w(s5) + w(s9)
    w(k21) =  w(s8) - w(s6)
    w(k31) =  w(s3) - w(s7)
    w(k41) =  w(s4) - w(s2)
    w(k22) =  w(s1) - w(s5) - w(s9)
    w(k32) =  w(s4) + w(s2)
    w(k42) =  w(s3) + w(s7)
    w(k33) = -w(s1) + w(s5) - w(s9)
    w(k43) =  w(s8) + w(s6)
    w(k44) = -w(s1) - w(s5) + w(s9)
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
      w(l3) = w(a2) * w(a2)
      w(l1) = SQRT(HALF * (SQRT(w(a2) * w(a2) - w(l2)) + w(a2)))
    else
      w(l1) = ONE
      do
        w(l0) = w(l1)
        w(l2) = w(l1) * w(l1)
        w(l3) = w(l2) * w(l1)
        w(l4) = w(l2) * w(l2)
        w(lf) = w(l4) - w(a2) * w(l2) + w(b8) * w(l1) + w(c1)
        w(lg) = w(l3) - w(a1) * w(l1) + w(b2)
        w(lg) = w(lg) + w(lg)
        w(lg) = w(lg) + w(lg)
        w(l1) = w(l1) - w(lf) / w(lg)
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
    rot(1) = w(k11) + w(k22) - w(k33) - w(k44)
    rot(2) = w(k32) - w(k41)
    rot(3) = w(k42) + w(k31)
    rot(4) = w(k32) + w(k41)
    rot(5) = w(k11) - w(k22) + w(k33) - w(k44)
    rot(6) = w(k43) - w(k21)
    rot(7) = w(k42) - w(k31)
    rot(8) = w(k43) + w(k21)
    rot(9) = w(k11) - w(k22) - w(k33) + w(k44)
!
  end subroutine quartenion_operation
!
end module mod_Kabsch

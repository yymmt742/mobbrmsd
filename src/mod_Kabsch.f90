!| Calculate the rotation matrix that minimizes |X-RY|^2 using the Kabsch-Umeyama algorithm.
!  Here, RR^T=I and det(R)=1 are satisfied.
module mod_Kabsch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, HALF => RHALF
  use mod_det
  implicit none
  private
  public :: Kabsch_worksize, Kabsch, estimate_rotation_matrix
!
  interface
    include 'dcopy.h'
    include 'dgemm.h'
    include 'dgesvd.h'
  end interface
!
  real(RK), parameter :: THRESHOLD = 1D-8
!
contains
!
  subroutine estimate_rotation_matrix(d, g, cov, rot, w)
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
    if(d==0)then
      W(1) = ZERO
    elseif(d==-1)then
      W(1) = ZERO
    elseif(d==-2)then
      W(1) = 10
    elseif(d==-3)then
      W(1) = 28
    elseif(d==1)then
      rot(1) = ONE
    elseif(d==2)then
      w(4) = (ONE + ONE) / g
      w(1) = w(4) * cov(1); w(2) = w(4) * cov(2)
      w(3) = w(4) * cov(3); w(4) = w(4) * cov(4)
      call quartenion_operation_d2(w(1), rot, w(5))
    elseif(d==3)then
      w(9) = (ONE + ONE) / g
      w(1) = w(9) * cov(1); w(2) = w(9) * cov(2); w(3) = w(9) * cov(3)
      w(4) = w(9) * cov(4); w(5) = w(9) * cov(5); w(6) = w(9) * cov(6)
      w(7) = w(9) * cov(7); w(8) = w(9) * cov(8); w(9) = w(9) * cov(9)
      call quartenion_operation_d3(w(1), rot, w(10))
    else
      call Kabsch(d, cov, rot, w)
    endif
!
  end subroutine estimate_rotation_matrix
!
  pure subroutine quartenion_operation_d2(s, rot, w)
    real(RK), intent(in)    :: s(*)
    !! target d*n array
    real(RK), intent(inout) :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout) :: w(*)
    integer(IK), parameter  :: c0 = 1
    integer(IK), parameter  :: k11 = 2, k41 = 3
    integer(IK), parameter  :: k22 = 1, k32 = 2
    integer(IK), parameter  :: q1 = 2, q4 = 1
    integer(IK), parameter  :: a = 4, c = 5, l = 6
    integer(IK), parameter  :: q0 = 3, p11 = 4, p44 = 5, p41 = 6
!
    w(k22) = s(1) - s(4)
    w(k32) = s(3) + s(2)
    w(c0) = w(k22) * w(k22) + w(k32) * w(k32)
    w(k11) = s(1) + s(4)
    w(k41) = s(3) - s(2)
!
    w(a) = s(1) * s(1) + s(2) * s(2) + s(3) * s(3) + s(4) * s(4)
    w(a) = w(a) + w(a)
    w(c) = (w(k11) * w(k11) + w(k41) * w(k41)) * w(c0)
    w(l) = w(c) + w(c)
    w(l) = w(l) + w(l)
    w(l) = SQRT(HALF * (SQRT(w(a) * w(a) - w(l)) + w(a)))
!
    w(c0) = w(c0) - w(l) * w(l)
    w(q1) = (w(k11) + w(l)) * w(c0)
    w(q4) = w(k41) * w(c0)
!
    w(p11) = w(q1) * w(q1)
    w(p44) = w(q4) * w(q4)
    w(q0) = ONE / (w(p11) + w(p44))
    w(p11) = w(p11) * w(q0)
    w(p44) = w(p44) * w(q0)
    w(q0) = w(q0) + w(q0)
    w(p41) = w(q1) * w(q4) * w(q0)
!
    rot(1) = w(p11) - w(p44)
    rot(2) = -w(p41)
    rot(3) = -rot(2)
    rot(4) = rot(1)
!
  end subroutine quartenion_operation_d2
!
  pure subroutine quartenion_operation_d3(s, rot, w)
    real(RK), intent(in)    :: s(*)
    !! target d*n array
    real(RK), intent(inout) :: rot(*)
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
      w(l1) = SQRT(HALF * (SQRT(w(a2) * w(a2) - w(l2)) + w(a2)))
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
  end subroutine quartenion_operation_d3
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
    elseif (d < 0) then
      s = ABS(d)
      call DGESVD('A', 'A', s, s, w, s, w, w, s, w, s, w, -1, info)
      w(1) = w(1) + s * s * 3 + s
      return
    end if
!
    if (d == 1)then
      rot(1) = ONE
      RETURN
    endif
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
end module mod_Kabsch

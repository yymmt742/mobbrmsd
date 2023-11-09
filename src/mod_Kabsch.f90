!| Calculate the rotation matrix that minimizes |X-RY|^2 using the Kabsch-Umeyama algorithm.
module mod_Kabsch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_svd
  use mod_cov
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
  pure elemental function optarg(l) result(res)
    logical, intent(in), optional :: l
    logical                       :: res
    if (PRESENT(l)) then
      res = l
    else
      res = .false.
    end if
  end function optarg
!
!| Calculate work array size for d*n matrix.
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
!| Calculate the rotation matrix.
  subroutine Kabsch(d, n, x, y, rot, w, row_major)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)       :: n
    !! matrix row dimension.
    real(RK), intent(in)          :: x(*)
    !! reference d*n array
    real(RK), intent(in)          :: y(*)
    !! target d*n array
    real(RK), intent(inout)       :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout)       :: w(*)
    !! work array, must be larger than Kabsch_worksize(d)
    !! if row_major, must be larger than Kabsch_worksize(n)
    logical, intent(in), optional :: row_major
    logical                       :: l
    integer(IK)                   :: d1, dd, s, u, vt, iw
!
    if (d < 1 .or. n < 1) RETURN
!
    l = optarg(row_major)
    if (l) then
!
      call cov_row_major(d, n, x, y, w)
      d1 = n
!
    else
!
      call cov(d, n, x, y, w)
      d1 = d
!
    end if
!
    dd = d1 * d1
    s = dd + 1
    u = s + d1
    vt = u + dd
    iw = vt + dd
!
    call svd(d1, w, w(s), w(u), w(vt), w(iw))
!
    call det_sign(d1, w(u:u+dd-1), w)
    w(s) = w(1)
    call det_sign(d1, w(vt:vt+dd-1), w)
    w(s) = w(s) * NINT(w(1))
!
    if (w(s) < ZERO) w(u + dd - d1:u + dd - 1) = -w(u + dd - d1:u + dd - 1)
!
    if (l) then
      call DGEMM('T', 'T', d1, d1, d1, ONE, w(vt), d1, w(u), d1, ZERO, w, d1)
    else
      call DGEMM('N', 'N', d1, d1, d1, ONE, w(u), d1, w(vt), d1, ZERO, w, d1)
    endif
!
    rot(:dd) = w(:dd)
!
  end subroutine Kabsch
!
end module mod_Kabsch

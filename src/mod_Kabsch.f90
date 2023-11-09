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
      res = svd_worksize(d) + d * d * 3 + d + 100
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
!   integer(IK)                   :: d1, dd, s, u, vt, iw
    real(RK)                      :: St(n, n)
    real(RK)                      :: s(d), u(d, d), vt(d, d), m(d, d)
    integer(IK)                   :: i,j
!
    if (d < 1 .or. n < 1) RETURN
!
!   if (optarg(row_major)) then
!
!     w(:n * n) = [MATMUL(TRANSPOSE(RESHAPE(x(:d * n), [d, n])), RESHAPE(y(:d * n), [d, n]))]
!     call cov_row_major(d, n, x, y, w)
!     d1 = n
!     dd = d1 * d1
!     u = dd + 1
!     vt = u + dd
!     s = vt + dd
!     iw = s + d1
!
!     call svd(d1, w, w(s), w(u), w(vt), w(iw))
!
!     call DGEMM('N', 'N', d1, d1, d1, ONE, w(u), d1, w(vt), d1, ZERO, w, d1)
!     call det_sign(d1, w)
!
!     if (w(1) < ZERO) w(u + dd - d1:u + dd - 1) = -w(u + dd - d1:u + dd - 1)
!     call DGEMM('N', 'N', d1, d1, d1, ONE, w(u), d1, w(vt), d1, ZERO, w, d1)
!     !rot(:dd) = [TRANSPOSE(RESHAPE(w(:dd), [d1, d1]))]
!     rot(:dd) = w(:dd)
!
!   else
!
!     d1 = d
!     dd = d1 * d1
!     u = dd + 1
!     vt = u + dd
!     s = vt + dd
!     iw = s + d1
!
print*,d,n
      call cov(d, n, x, y, m)
      !call cov(d1, n, x, y, w)
!
      call svd(d, m, s, u, vt, w)
      !call svd(d1, w, w(s), w(u), w(vt), w(iw))
!
!     call DGEMM('N', 'N', d1, d1, d1, ONE, w(u), d1, w(vt), d1, ZERO, w, d1)
!     w(:dd) = w(u:u+dd-1)
!
      w(:d * d) = [u]
      call det_sign(d, w)
      w(d * d + 1) = w(1)
      w(:d * d) = [vt]
      call det_sign(d, w)
      w(1) = w(d * d + 1) * w(1)
      print*,w(1)
      if (w(1) < ZERO) u(:, d) = -u(:, d)
      print'(6f9.3)',MATMUL(vt, u)
      rot(:d * d) = [MATMUL(vt, u)]
      print'(6f9.3)',rot(:d*d)
!     call det_sign(d1, w)
!     w(s) = w(1)
!     w(:dd) = w(vt:vt+dd-1)
!     call det_sign(d1, w)
!     w(s) = w(s) * w(1)
!
!     if (w(s) < ZERO) w(u + dd - d1:u + dd - 1) = -w(u + dd - d1:u + dd - 1)
!     call DGEMM('N', 'N', d1, d1, d1, ONE, w(u), d1, w(vt), d1, ZERO, w, d1)
!     rot(:dd) = w(:dd)
!
!   end if
!
  end subroutine Kabsch
!
end module mod_Kabsch

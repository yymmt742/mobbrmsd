!| Calculate the rotation matrix that minimizes \(|\mathbf{X}-\mathbf{R}\mathbf{Y}|^2\). <br>
!  Here, \(\mathbf{R}\mathbf{R}^\top=\mathbf{I}\) and \(\det(\mathbf{R})=1\) are satisfied. <br>
!  This code is based on the Kabsch algorithm.
!  doi : [10.1107/S0567739476001873](https://scripts.iucr.org/cgi-bin/paper?S0567739476001873)
module mod_rotation
  use mod_kinds, only: IK, RK
  use mod_dimspec_functions, only: D, DD
  implicit none
  private
  public :: sdmin_worksize
  public :: estimate_rcmax
  public :: estimate_sdmin
  public :: rotation_worksize
  public :: estimate_rotation
!
  interface
#ifdef USE_REAL32
    include 'sgemm.h'
    include 'sgesvd.h'
    include 'sgetrf.h'
#else
    include 'dgemm.h'
    include 'dgesvd.h'
    include 'dgetrf.h'
#endif
  end interface
!
  real(RK), parameter    :: ZERO = 0.0_RK
  real(RK), parameter    :: ONE = 1.0_RK
!
contains
!
!| Inquire function for memory size of estimate_sdmin.
  pure elemental function sdmin_worksize() result(res)
    integer(IK) :: res
    res = worksize_Kabsch() + DD + 1
  end function sdmin_worksize
!
!| Compute \(\max_{\mathbf{R}}\text{tr}[\mathbf{R}\mathbf{C}]\).
  pure subroutine estimate_rcmax(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_sdmin().
    call Kabsch(cov, w(2), w(DD + 2))
    w(1) = dot(DD, cov, w(2))
  end subroutine estimate_rcmax
!
!| Compute \(\min_{\mathbf{R}}(g-2\text{tr}[\mathbf{R}\mathbf{C}])\),
!  where \(g = tr[\mathbf{X}\mathbf{X}^\top] + tr[\mathbf{Y}\mathbf{Y}^\top]\)
  pure subroutine estimate_sdmin(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_sdmin().
!
    call Kabsch(cov, w(2), w(DD + 2))
    w(1) = dot(DD, cov, w(2))
    w(1) = w(1) + w(1)
    w(1) = g - w(1)
!
  end subroutine estimate_sdmin
!
!| Inquire function for memory size of estimate_rotation().
  pure elemental function rotation_worksize() result(res)
    integer(IK) :: res
    res = worksize_Kabsch()
  end function rotation_worksize
!
!| Compute the transpose of rotation matrix \(\mathbf{R}\)
!  that maximize \(\text{tr}[\mathbf{R}\mathbf{Y}\mathbf{X}^\top]\).
  pure subroutine estimate_rotation(g, cov, rot, w)
    real(RK), intent(in)    :: g
    !! \(g=\text{tr}[\mathbf{X}\mathbf{X}^\top}]+\text{tr}[\mathbf{X}\mathbf{X}^\top}]\).
    real(RK), intent(in)    :: cov(*)
    !! covariance matrix \(\mathbf{C}=\mathbf{Y}\mathbf{X}^\top\in\mathbb R^{D\times D}\)
    !! rotation
    real(RK), intent(inout) :: rot(*)
    !! rotation matrix, \(\mathbf{R}\in\mathbb R^{D\times D}\).
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_rotation().
    call Kabsch(cov, rot, w)
  end subroutine estimate_rotation
!
!| work array size for Kabsch algorithm.
  pure elemental function worksize_Kabsch() result(res)
    real(RK)    :: w(1)
    integer(IK) :: res, info
#ifdef USE_REAL32
    call SGESVD('A', 'A', D, D, w, D, w, w, D, w, D, w, -1, info)
#else
    call DGESVD('A', 'A', D, D, w, D, w, w, D, w, D, w, -1, info)
#endif
    res = NINT(w(1)) + DD * 3 + D
  end function worksize_Kabsch
!
  pure subroutine Kabsch(cov, rot, w)
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than Kabsch_worksize(d)
    !! if row_major, must be larger than Kabsch_worksize(n)
    integer(IK), parameter  :: m = 1
    integer(IK)             :: s, u, vt, iw, lw, info
!
    u = m + DD
    vt = u + DD
    s = vt + DD
    iw = s + D
#ifdef USE_REAL32
    call SGESVD('A', 'A', D, D, w(m), D, w(s), w(u), D, w(vt), D, w(iw), -1, info)
    lw = NINT(w(iw))
!
    call copy(DD, cov, w(m))
    call SGESVD('A', 'A', D, D, w(m), D, w(s), w(u), D, w(vt), D, w(iw), lw, info)
!
    call SGEMM('N', 'N', D, D, D, ONE, w(u), D, w(vt), D, ZERO, w(s), D)
    call det_sign(w(s))
    if (w(s) < ZERO) call neg(d, w(u + DD - D))
    call SGEMM('N', 'N', D, D, D, ONE, w(u), D, w(vt), D, ZERO, w(s), D)
    call copy(DD, w(s), rot(1))
#else
    call DGESVD('A', 'A', D, D, w(m), D, w(s), w(u), D, w(vt), D, w(iw), -1, info)
    lw = NINT(w(iw))
!
    call copy(DD, cov, w(m))
    call DGESVD('A', 'A', D, D, w(m), D, w(s), w(u), D, w(vt), D, w(iw), lw, info)
!
    call DGEMM('N', 'N', D, D, D, ONE, w(u), D, w(vt), D, ZERO, w(s), D)
    call det_sign(w(s))
    if (w(s) < ZERO) call neg(d, w(u + DD - D))
    call DGEMM('N', 'N', D, D, D, ONE, w(u), D, w(vt), D, ZERO, w(s), D)
    call copy(DD, w(s), rot(1))
#endif
  end subroutine Kabsch
!
!| calculate determinant sign of square matrix x, with leading dimension.
  pure subroutine det_sign(x)
    real(RK), intent(inout) :: x(*)
     !! square matrix, on exit, x(1) is assigned the determinant sign of x, <br>
     !! and the other elements are undefined.
!
    if (D < 1) then
      return
    elseif (D == 1) then
      x(1) = SIGN(ONE, x(1))
    elseif (D == 2) then
      x(1) = SIGN(ONE, x(1) * x(4) - x(2) * x(3))
    elseif (D == 3) then
      x(1) = SIGN(ONE, x(1) * (x(5) * x(9) - x(8) * x(6)) +&
        &              x(4) * (x(8) * x(3) - x(2) * x(9)) +&
        &              x(7) * (x(2) * x(6) - x(5) * x(3)))
    else
      block
        integer(IK) :: i, j, k, ipiv(D)
#ifdef USE_REAL32
        call SGETRF(D, D, x, D, ipiv, j)
#else
        call DGETRF(D, D, x, D, ipiv, j)
#endif
        ipiv(1) = COUNT([(ipiv(i) == i, i=1, D)])
        j = 1
        k = D + 1
        do i = 1, D
          if (x(j) <= ZERO) ipiv(1) = ipiv(1) + 1
          j = j + k
        end do
        if (MODULO(ipiv(1), 2) == 0) then
          x(1) = ONE
        else
          x(1) = -ONE
        end if
      end block
    end if
!
  end subroutine det_sign
!
  pure subroutine neg(N, X)
    integer(IK), intent(in) :: N
    real(RK), intent(inout) :: X(*)
    integer(IK)             :: i
    do concurrent(i=1:N)
      X(i) = -X(i)
    end do
  end subroutine neg
!
  pure function dot(N, X, Y) result(res)
    integer(IK), intent(in) :: N
    real(RK), intent(in)    :: X(*), Y(*)
    real(RK)                :: res, tmp
    integer(IK)             :: i
    res = ZERO
    tmp = ZERO
    do i = 2, N, 2
      res = res + X(i - 1) * Y(i - 1)
      tmp = tmp + X(i - 0) * Y(i - 0)
    end do
    res = res + tmp
    if (MODULO(N, 2) == 1) res = res + X(N) * Y(N)
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


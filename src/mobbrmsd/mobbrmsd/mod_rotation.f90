!| Calculate the rotation matrix that minimizes |X-RY|^2 for D=3. <br>
!  Here, RR^T=I and det(R)=1 are satisfied.
module mod_rotation
  use mod_kinds, only: IK, RK
  use blas_lapack_interface, only: D, DD, DGEMM
  implicit none
  private
  public :: sdmin_worksize
  public :: estimate_sdmin
  public :: rotation_worksize
  public :: estimate_rotation
!
  interface
    include 'dgesvd.h'
    include 'dgetrf.h'
  end interface
!
  real(RK), parameter    :: ZERO = 0.0_RK
  real(RK), parameter    :: ONE = 1.0_RK
  real(RK), parameter    :: THRESHOLD = 1E-8_RK
  integer(IK), parameter :: MAXITER = 100000
!
contains
!
!| Inquire function for memory size of estimate_sdmin.
  pure elemental function sdmin_worksize() result(res)
    integer(IK) :: res
    res = worksize_Kabsch() + DD + 1
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
    call Kabsch(cov, w(2), w(DD + 2))
    w(1) = dot(DD, cov, w(2))
    w(1) = w(1) + w(1)
    w(1) = g - w(1)
!
  end subroutine estimate_sdmin
!
!| Inquire function for memory size of rotation.
  pure elemental function rotation_worksize() result(res)
    integer(IK) :: res
    res = worksize_Kabsch()
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
    call Kabsch(cov, rot, w)
!
  end subroutine estimate_rotation
!
  !| work array size for Kabsch algorithm.
  pure elemental function worksize_Kabsch() result(res)
    real(RK)    :: w(1)
    integer(IK) :: res, info
!
    call DGESVD('A', 'A', D, D, w, D, w, w, D, w, D, w, -1, info)
    res = NINT(w(1)) + DD * 3 + D
!
  end function worksize_Kabsch
!
!| Calculate the rotation matrix R^T from covariance matrix.
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
!
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
!
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
         call DGETRF(D, D, x, D, ipiv, j)
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
  pure subroutine neg(N, X)
    integer(IK), intent(in) :: N
    real(RK), intent(inout) :: X(*)
    integer(IK)             :: i
    do concurrent(i=1:N)
      X(i) = -X(i)
    end do
  end subroutine neg
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


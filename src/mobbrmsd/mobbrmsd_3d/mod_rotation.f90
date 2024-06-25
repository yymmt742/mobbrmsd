!| Calculate the rotation matrix that minimizes \(|\mathbf{X}-\mathbf{R}\mathbf{Y}|^2\) for \(D=3\). <br>
!  Here, \(\mathbf{R}\mathbf{R}^\top=\mathbf{I}\) and \(\det(\mathbf{R})=1\) are satisfied. <br>
!  This code is based on the method of Coutsias et.al.
!  doi : [10.1002/jcc.25802](https://onlinelibrary.wiley.com/doi/10.1002/jcc.25802)
module mod_rotation
  use mod_kinds, only: IK, RK, R8
  implicit none
  private
  public :: sdmin_worksize
  public :: estimate_rcmax
  public :: estimate_sdmin
  public :: rotation_worksize
  public :: estimate_rotation
!
  real(R8), parameter    :: ZERO = 0.0_R8
  real(R8), parameter    :: HALF = 0.5_R8
  real(R8), parameter    :: ONE = 1.0_R8
  real(R8), parameter    :: TWO = 2.0_R8
  real(R8), parameter    :: FOUR = 4.0_R8
  real(R8), parameter    :: EIGHT = 8.0_R8
  real(R8), parameter    :: THRESHOLD = 1E-14_R8
  real(R8), parameter    :: DEGENERACY = 1E-6_R8
  integer(IK), parameter :: MAXITER = 100000
!
contains
!
!| Inquire function for memory size of estimate_sdmin.
  pure elemental function sdmin_worksize() result(res)
    integer(IK) :: res
#ifdef USE_REAL32
    res = 1
#else
    res = 9
#endif
  end function sdmin_worksize
!
!| Compute \(\min_{R}\text{tr}[\mathbf{R}\mathbf{C}]\).
  pure subroutine estimate_rcmax(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_sdmin().
#ifdef USE_REAL32
    real(R8)                :: w8(9)
    call find_lambda_max(g, cov, w8)
    w(1) = w8(1)
#else
    call find_lambda_max(g, cov, w)
#endif
  end subroutine estimate_rcmax
!
!| Compute the least-squares sum_i^n |x_i-Ry_i|^2 from cov = YX^T and g = tr[XX^T] + tr[YY^T].
  pure subroutine estimate_sdmin(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_sdmin().
#ifdef USE_REAL32
    real(R8)                :: w8(9)
    call find_lambda_max(g, cov, w8)
    w(1) = real(real(g, R8) - (w8(1) + w8(1)), RK)
#else
    call find_lambda_max(g, cov, w)
    w(1) = g - w(1) - w(1)
#endif
  end subroutine estimate_sdmin
!
!| Inquire function for memory size of rotation.
  pure elemental function rotation_worksize() result(res)
    integer(IK) :: res
#ifdef USE_REAL32
    res = 1
#else
    res = 18
#endif
  end function rotation_worksize
!
!| Compute the transpose rotation matrix for minimize tr[CR] from cov = YX^T and g = tr[XX^T] + tr[YY^T].
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  pure subroutine estimate_rotation(g, cov, rot, w)
    real(RK), intent(in)    :: g
    !! g = tr[XX^T] + tr[YY^T]
    real(RK), intent(in)    :: cov(*)
    !! covariance dxd matrix, YX^T
    real(RK), intent(inout) :: rot(*)
    !! rotation dxd matrix
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_rotation().
#ifdef USE_REAL32
    real(R8)                :: w8(18)
#endif
    if (g < THRESHOLD) then
      rot(1) = ONE; rot(2) = ZERO; rot(3) = ZERO
      rot(4) = ZERO; rot(5) = ONE; rot(6) = ZERO
      rot(7) = ZERO; rot(8) = ZERO; rot(9) = ONE
      return
    end if
#ifdef USE_REAL32
    call find_rotmatrix(g, cov, w8, rot)
#else
    call find_rotmatrix(g, cov, w, rot)
#endif
  end subroutine estimate_rotation
!
!| Compute rotation matrix. <br>
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  pure subroutine find_rotmatrix(g, cov, w, rot)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(R8), intent(inout) :: w(*)
    !! work array
    real(RK), intent(inout) :: rot(*)
    !! rotation matrix
    integer(IK), parameter  :: l2 = 1, l3 = 2, l4 = 3
    integer(IK), parameter  :: y2 = 4, y3 = 5, y4 = 6
    integer(IK), parameter  :: dg1 = 1, dg2 = 2, dg3 = 3, dg4 = 18
    integer(IK), parameter  :: aa1 = 1, aa2 = 2, aa3 = 2, aa4 = 1
    integer(IK), parameter  :: a11 = 4, a21 = 5, a31 = 6, a41 = 12
    integer(IK), parameter  :: a22 = 13, a32 = 14, a42 = 15, a33 = 16, a43 = 17, a44 = 18
    integer(IK), parameter  :: s22 = 7, s23 = 8, s24 = 9, s33 = 10, s34 = 11, s44 = 12
    integer(IK), parameter  :: v1 = 3, v2 = 4, v3 = 5, v4 = 6
    integer(IK), parameter  :: v11 = 7, v21 = 8, v31 = 9, v41 = 10
    integer(IK), parameter  :: v22 = 11, v32 = 12, v42 = 13, v33 = 14, v43 = 15, v44 = 16
!
    call find_lambda_max(g, cov, w)
!
!   A = L - lambda_max * I
!
!   L = (R11+R22+R33  R23-R32      R31-R13      R12-R21    )
!       (R23-R32      R11-R22-R33  R12+R21      R13+R31    )
!       (R31-R13      R12+R21     -R11+R22-R33  R23+R32    )
!       (R12-R21      R13+R31      R23+R32     -R11-R22+R33)
!
    w(dg4) = cov(9) - w(1)
    w(dg3) = cov(9) + w(1)
    w(dg2) = cov(1) - cov(5)
    w(dg1) = cov(1) + cov(5)
    w(a11) = w(dg1) + w(dg4)
    w(a21) = cov(8) - cov(6)
    w(a31) = cov(3) - cov(7)
    w(a41) = cov(4) - cov(2)
    w(a22) = w(dg2) - w(dg3)
    w(a32) = cov(4) + cov(2)
    w(a42) = cov(7) + cov(3)
    w(a33) = -w(dg2) - w(dg3)
    w(a43) = cov(8) + cov(6)
    w(a44) = -w(dg1) + w(dg4)
!
    w(aa1) = ABS(w(a11))
    w(aa3) = ABS(w(a33))
!
    if (w(aa1) > w(aa3)) then
      w(aa2) = ABS(w(a22))
      if (w(aa1) > w(aa2)) then
        w(l4) = ONE / w(a11)
        w(l2) = w(a21) * w(l4)
        w(l3) = w(a31) * w(l4)
        w(l4) = w(a41) * w(l4)
        w(s22) = w(a22) - w(a21) * w(l2)
        w(s23) = w(a32) - w(a21) * w(l3)
        w(s33) = w(a33) - w(a31) * w(l3)
        w(s24) = w(a42) - w(a21) * w(l4)
        w(s34) = w(a43) - w(a31) * w(l4)
        w(s44) = w(a44) - w(a41) * w(l4)
        call find_null_vector(w)
        w(v1) = -w(l2) * w(y2) - w(l3) * w(Y3) - w(l4) * w(y4)
      else
        w(l4) = ONE / w(a22)
        w(l2) = w(A21) * w(l4)
        w(l3) = w(A32) * w(l4)
        w(l4) = w(A42) * w(l4)
        w(s22) = w(a11) - w(a21) * w(l2)
        w(s23) = w(a31) - w(a21) * w(l3)
        w(s33) = w(a33) - w(a32) * w(l3)
        w(s24) = w(a41) - w(a21) * w(l4)
        w(s34) = w(a43) - w(a32) * w(l4)
        w(s44) = w(a44) - w(a42) * w(l4)
        call find_null_vector(w)
        w(l2) = -w(l2) * w(y2) - w(l3) * w(Y3) - w(l4) * w(y4)
        w(v1) = w(y2)
        w(v2) = w(l2)
      end if
    else
      w(aa4) = ABS(w(a44))
      if (w(aa3) > w(aa4)) then
        w(l4) = ONE / w(a33)
        w(l2) = w(a31) * w(l4)
        w(l3) = w(a32) * w(l4)
        w(l4) = w(a43) * w(l4)
        w(s22) = w(a11) - w(a31) * w(l2)
        w(s23) = w(a21) - w(a31) * w(l3)
        w(s33) = w(a22) - w(a32) * w(l3)
        w(s24) = w(a41) - w(a31) * w(l4)
        w(s34) = w(a42) - w(a32) * w(l4)
        w(s44) = w(a44) - w(a43) * w(l4)
        call find_null_vector(w)
        w(l2) = -w(l2) * w(y2) - w(l3) * w(y3) - w(l4) * w(y4)
        w(v1) = w(y2)
        w(v2) = w(y3)
        w(v3) = w(l2)
      else
        w(l4) = ONE / w(a44)
        w(l2) = w(a41) * w(l4)
        w(l3) = w(a42) * w(l4)
        w(l4) = w(a43) * w(l4)
        w(s22) = w(a11) - w(a41) * w(l2)
        w(s23) = w(a21) - w(a41) * w(l3)
        w(s33) = w(a22) - w(a42) * w(l3)
        w(s24) = w(a31) - w(a41) * w(l4)
        w(s34) = w(a32) - w(a42) * w(l4)
        w(s44) = w(a33) - w(a43) * w(l4)
        call find_null_vector(w)
        w(l2) = -w(l2) * w(y2) - w(l3) * w(Y3) - w(l4) * w(y4)
        w(v1) = w(y2)
        w(v2) = w(y3)
        w(v3) = w(y4)
        w(v4) = w(l2)
      end if
    end if
!
    w(v11) = w(v1) * w(v1)
    w(v22) = w(v2) * w(v2)
    w(v33) = w(v3) * w(v3)
    w(v44) = w(v4) * w(v4)
    w(v21) = w(v1) * w(v2)
    w(v31) = w(v1) * w(v3)
    w(v41) = w(v1) * w(v4)
    w(v32) = w(v2) * w(v3)
    w(v42) = w(v2) * w(v4)
    w(v43) = w(v3) * w(v4)
!
    w(l2) = ONE / (w(v11) + w(v22) + w(v33) + w(v44))
    w(l3) = w(l2) + w(l2)
!
    rot(1) = w(l2) * (w(v11) + w(v22) - w(v33) - w(v44))
    rot(2) = w(l3) * (w(v32) - w(v41))
    rot(3) = w(l3) * (w(v42) + w(v31))
    rot(4) = w(l3) * (w(v32) + w(v41))
    rot(5) = w(l2) * (w(v11) - w(v22) + w(v33) - w(v44))
    rot(6) = w(l3) * (w(v43) - w(v21))
    rot(7) = w(l3) * (w(v42) - w(v31))
    rot(8) = w(l3) * (w(v43) + w(v21))
    rot(9) = w(l2) * (w(v11) - w(v22) - w(v33) + w(v44))
  end subroutine find_rotmatrix
!
!| Compute maximum eigen value of S. <br>
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  pure subroutine find_lambda_max(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(R8), intent(inout) :: w(*)
    integer(IK), parameter  :: k1 = 2, k0 = 3, k2 = 4
    integer(IK), parameter  :: xk = 1, s = 5, xx = 6
    integer(IK), parameter  :: a = 8, f = 7, df = 8, gt = 9
    integer(IK), parameter  :: d11 = 3, d22 = 4, d33 = 5, d21 = 6, d31 = 7, d32 = 8
    integer(IK), parameter  :: a1 = 6, a2 = 5
    integer(IK)             :: k
!
    if (g < THRESHOLD) then
      w(1) = ZERO
      return
    end if
!
!   K1 = - 8 det|R|
!&<
    w(k1) = - cov(1) * (cov(5) * cov(9) - cov(8) * cov(6)) &
   &        - cov(4) * (cov(8) * cov(3) - cov(2) * cov(9)) &
   &        - cov(7) * (cov(2) * cov(6) - cov(5) * cov(3))
    w(k1) = EIGHT * w(k1)
!>&
!   D = RR^T
!
    w(d11) = cov(1) * cov(1) + cov(2) * cov(2) + cov(3) * cov(3)
    w(d22) = cov(4) * cov(4) + cov(5) * cov(5) + cov(6) * cov(6)
    w(d33) = cov(7) * cov(7) + cov(8) * cov(8) + cov(9) * cov(9)
    w(d21) = cov(4) * cov(1) + cov(5) * cov(2) + cov(6) * cov(3)
    w(d31) = cov(7) * cov(1) + cov(8) * cov(2) + cov(9) * cov(3)
    w(d32) = cov(4) * cov(7) + cov(5) * cov(8) + cov(6) * cov(9)
!
!   A1 = D11 * (D22+D33) + D22 * D33 - D12**2 - D13**2 - D23**2
!   A2 = tr[D]
!
    w(a1) = w(d11) * (w(d22) + w(d33)) + w(d22) * w(d33) - w(d21)**2 - w(d31)**2 - w(d32)**2
    w(a1) = FOUR * w(a1)
    w(a2) = w(d11) + w(d22) + w(d33)
!
!   K0 = 16 * (A2**2 - 4*A1)
!   K2 = -8 * A2
!
    w(k0) = w(a2) * w(a2) - w(a1)
    w(k2) = -TWO * w(a2)
!
    w(gt) = ABS(g * THRESHOLD)
    if (ABS(w(k1)) < w(gt)) then
!
!     find solution of x**4 + K2 * x**2 + K0 = 0
!     that is x**2 = (K2 + sqrt{K2**2 - 4 * K0}) / 2
!
      w(xk) = w(k2) * w(k2) - FOUR * w(k0)
!
      if (w(xk) < ZERO) then
        if (w(k2) < ZERO) then
          w(xk) = SQRT(-HALF * w(k2))
        else
          w(xk) = ZERO
        end if
      else
        w(xk) = SQRT(w(xk)) - w(k2)
        if (w(xk) < ZERO) then
          w(xk) = ZERO
        else
          w(xk) = SQRT(HALF * w(xk))
        end if
      end if
!
    else
!
!     find solution of x**4 + K2 * x**2 + K1 * x + K0 = 0
!
      w(s) = ONE
      w(xk) = HALF * g
!
      do k = 1, MAXITER
        w(xx) = w(xk) * w(xk)
        w(a) = w(k2) + w(xx)
        w(f) = w(a) * w(xx) + w(k1) * w(xk) + w(k0)
        w(df) = w(k1) + (w(xk) + w(xk)) * (w(a) + w(xx))
        if (ABS(w(df)) < w(gt) .and. ABS(w(f)) < w(gt)) exit
        w(s) = w(f) / w(df)
        w(xk) = w(xk) - w(s)
        if (w(s) < w(gt) * ABS(w(xk))) exit
      end do
    end if
!
  end subroutine find_lambda_max
!
!| Find null vector of S. <br>
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  pure subroutine find_null_vector(w)
    real(R8), intent(inout) :: w(*)
    integer(IK), parameter  :: y2 = 4, y3 = 5, y4 = 6
    integer(IK), parameter  :: s22 = 7, s23 = 8, s24 = 9
    integer(IK), parameter  :: s33 = 10, s34 = 11, s44 = 12
    integer(IK), parameter  :: m22 = 13, m23 = 14, m24 = 15
    integer(IK), parameter  :: m33 = 16, m34 = 17, m44 = 18
    integer(IK), parameter  :: mm2 = 4, mm3 = 5, mm4 = 6
!
    w(m22) = w(s33) * w(s44) - w(s34) * w(s34)
    w(m23) = w(s34) * w(s24) - w(s23) * w(s44)
    w(m24) = w(s23) * w(s34) - w(s33) * w(s24)
!
    w(mm2) = w(m22) * w(m22)
    w(mm3) = w(m23) * w(m23)
    w(mm4) = w(m24) * w(m24)
    w(m33) = w(mm2) + w(mm3) + w(mm4)
!
    if (w(m33) < DEGENERACY) then
      w(m44) = w(m33) - w(mm2)
      w(m33) = w(s22) * w(s44) - w(mm4)
      w(m34) = w(s22) * w(s34) - w(s23) * w(s24)
      if (w(m44) < DEGENERACY) then
        w(mm4) = w(m44) - w(mm3)
        w(m44) = w(s22) * w(s33) - w(s23) * w(s23)
        if (w(mm4) < DEGENERACY) then
          ! double degeneracy
          if (ABS(w(s22)) > DEGENERACY) then
            w(y2) = -w(s23)
            w(y3) = w(s22)
            w(y4) = ZERO
          elseif (ABS(w(s33)) > DEGENERACY) then
            w(y2) = w(s33)
            w(y3) = -w(s23)
            w(y4) = ZERO
          elseif (ABS(w(s44)) > DEGENERACY) then
            w(y2) = w(s44)
            w(y3) = ZERO
            w(y4) = -w(s24)
          else
            w(y2) = ONE
            w(y3) = ZERO
            w(y4) = ZERO
          end if
        else
          ! Triple degeneracy
          w(y2) = ZERO
          w(y3) = ZERO
          w(y4) = w(m44)
        end if
      else
        w(y2) = ZERO
        w(y3) = w(m33)
        w(y4) = w(m34)
      end if
    else
      w(y2) = w(m22)
      w(y3) = w(m23)
      w(y4) = w(m24)
    end if
!
  end subroutine find_null_vector
!
end module mod_rotation


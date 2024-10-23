!| Calculate the rotation matrix that minimizes \(|\mathbf{X}-\mathbf{R}\mathbf{Y}|^2\) for \(D=3\). <br>
!  Here, \(\mathbf{R}\mathbf{R}^\top=\mathbf{I}\) and \(\det(\mathbf{R})=1\) are satisfied. <br>
!  This code is based on the method of Coutsias et.al.
!  doi : [10.1002/jcc.25802](https://onlinelibrary.wiley.com/doi/10.1002/jcc.25802)
module mod_rotation
  use mod_kinds, only: IK, RK
  implicit none
  private
  public :: sdmin_worksize
  public :: estimate_rcmax
  public :: estimate_sdmin
  public :: rotation_worksize
  public :: estimate_rotation
!
  real(RK), parameter    :: ZERO = 0.0_RK
  real(RK), parameter    :: QUARTER = 0.25_RK
  real(RK), parameter    :: HALF = 0.5_RK
  real(RK), parameter    :: ONE = 1.0_RK
  real(RK), parameter    :: TWO = 2.0_RK
  real(RK), parameter    :: THREE = 3.0_RK
  real(RK), parameter    :: FOUR = 4.0_RK
  real(RK), parameter    :: SIX = 6.0_RK
  real(RK), parameter    :: EIGHT = 8.0_RK
  real(RK), parameter    :: SIXTEEN = 16.0_RK
  real(RK), parameter    :: ONETHIRD = 1.0_RK / 3.0_RK
  real(RK), parameter    :: ONESIX = 1.0_RK / 6.0_RK
  real(RK), parameter    :: ONEQUARTER = 0.25_RK
  real(RK), parameter    :: ONEEIGHT = 0.125_RK
  real(RK), parameter    :: EIGHTNINE = (8.0_RK / 9.0_RK)
  real(RK), parameter    :: TWOTHIRD = (2.0_RK / 3.0_RK)
  real(RK), parameter    :: FOURTHIRD = (4.0_RK / 3.0_RK)
  real(RK), parameter    :: SQRT3 = SQRT(3.0_RK)
  real(RK), parameter    :: HALFSQRT3 = HALF * SQRT3
#ifdef USE_REAL32
  real(RK), parameter    :: THRESHOLD = 1E-6_RK
  real(RK), parameter    :: DEGENERACY = 1E-4_RK
  real(RK), parameter    :: DEGENERACY1 = 1E-3_RK
  real(RK), parameter    :: DEGENERACY2 = 1E-6_RK
#else
  real(RK), parameter    :: THRESHOLD = 1E-12_RK
  real(RK), parameter    :: DEGENERACY = 1E-6_RK
  real(RK), parameter    :: DEGENERACY1 = 1E-4_RK
  real(RK), parameter    :: DEGENERACY2 = 1E-12_RK
#endif
  integer(IK), parameter :: MAXITER = 1000
!
contains
!
!| Inquire function for memory size of estimate_sdmin.
  pure elemental function sdmin_worksize() result(res)
    integer(IK) :: res
    res = 9
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
    call find_lambda_max(g, cov, w)
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
    call find_lambda_max(g, cov, w)
    !w(1) = g * (ONE - w(1))
    w(1) = g - TWO * w(1)
  end subroutine estimate_sdmin
!
!| Inquire function for memory size of rotation.
  pure elemental function rotation_worksize() result(res)
    integer(IK) :: res
    res = 18
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
    integer(IK), parameter  :: l2 = 1, l3 = 2, l4 = 3
    integer(IK), parameter  :: y2 = 4, y3 = 5, y4 = 6
    integer(IK), parameter  :: dg1 = 1, dg2 = 2, dg3 = 3, dg4 = 18, nrm = 7
    integer(IK), parameter  :: aa1 = 1, aa2 = 2, aa3 = 2, aa4 = 1
    integer(IK), parameter  :: a11 = 4, a21 = 5, a31 = 6, a41 = 12
    integer(IK), parameter  :: a22 = 13, a32 = 14, a42 = 15, a33 = 16, a43 = 17, a44 = 18
    integer(IK), parameter  :: s22 = 7, s23 = 8, s24 = 9, s33 = 10, s34 = 11, s44 = 12
    integer(IK), parameter  :: v1 = 3, v2 = 4, v3 = 5, v4 = 6
    integer(IK), parameter  :: v11 = 7, v21 = 8, v31 = 9, v41 = 10
    integer(IK), parameter  :: v22 = 11, v32 = 12, v42 = 13, v33 = 14, v43 = 15, v44 = 16
    if (g < THRESHOLD) then
      rot(1) = ONE; rot(2) = ZERO; rot(3) = ZERO
      rot(4) = ZERO; rot(5) = ONE; rot(6) = ZERO
      rot(7) = ZERO; rot(8) = ZERO; rot(9) = ONE
      return
    end if
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
    w(l3) = w(l2) * TWO
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
  end subroutine estimate_rotation
!
!| Compute maximum eigen value of S. <br>
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  pure subroutine find_lambda_max(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    integer(IK), parameter  :: x = 1, a = 2, b = 3
    integer(IK), parameter  :: w1 = 5, w2 = 6, w3 = 7, w4 = 8, w5 = 9
    integer(IK), parameter  :: d11 = 4, d22 = 5, d33 = 6, d21 = 7, d31 = 8, d32 = 9
    integer(IK), parameter  :: p = 1, r = 2, q = 3, sr = 4, rr = 5
!   integer(IK), parameter  :: xk = 1, k1 = 2, k0 = 3, sa = 4
!   integer(IK), parameter  :: a1 = 6, a2 = 5
!   integer(IK), parameter  :: m0 = 5
!   integer(IK), parameter  :: r = 5, q = 6, h = 7, s = 8
!   integer(IK), parameter  :: xx = 5, f = 6, df = 7
!
    if (g < THRESHOLD) then
      w(x) = ZERO
      return
    end if
!
!   K1 = - 8 det|R|
!   A2 = tr[D]
    call det3(g, cov, w(q))
!    w(k1) = -EIGHT * w(1)
!
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
!
    w(p) = (w(d22) + w(d33))
    w(r) = w(d11) + w(p)
    if (w(r) < THRESHOLD) then
      w(x) = ZERO
      return
    end if
    w(p) = w(p) * w(d11) + w(d22) * w(d33) - w(d21)**2 - w(d31)**2 - w(d32)**2
    w(sr) = SQRT(w(r))
    w(rr) = FOUR / (w(r) * w(r))
    w(a) = -TWO * w(q) * w(rr) * w(sr)
    w(b) = ONE - w(p) * w(rr)
    call find_a_quartic_root(w(a), w(b), w(x), w(w1), w(w2), w(w3), w(w4), w(w5))
    w(x) = w(x) * w(sr)
!
!   w(a1) = w(d11) * (w(d22) + w(d33)) + w(d22) * w(d33) - w(d21)**2 - w(d31)**2 - w(d32)**2
!   w(a1) = FOUR * w(a1)
!
!   w(a2) = w(d11) + w(d22) + w(d33)
!
!   K2 = -2 * A2
!   K1 = -8 det|R|
!   K0 = A2**2 - 4*A1
!
!   Rescale to X_  = X * SQRT(A2), 0 <= X_ <= g / (2 * SQRT(A2))
!              K2_ = K2 / A2                 = -2
!              K1_ = SQRT(A2) * K1 / (A2**2) = - 8 * det|R| / (A2**(3/2))
!              K0_ = K0 / A2**2              = 1 - 4 * A1 / A2**2
!
!   w(sa) = SQRT(w(a2))
!   w(k0) = ONE / (w(a2) * w(a2))
!
!   w(k1) = w(k0) * w(sa) * w(k1)
!   w(k0) = ONE - w(k0) * w(a1)
!
!   normalize
!
!   w(m0) = w(k1)**2
!   w(m0) = (-27._RK * w(m0) - 288._RK * w(k0) * 32._RK) * w(m0) + 256._RK * w(k0) * (w(k0) - ONE)**2
!   if (ABS(w(m0)) < DEGENERACY1) then
!     ! Third order Taylor expansion around x=1.
!     ! f3(x) / 4 = x**3 - 2 * x**2 + (1 + k1/4) * x + (k0 - 1)/4
!     if (ABS(w(k0) - ONE) < DEGENERACY2) then
!       ! Solve x**3 - 2 * x**2 + (1 + k1/4) * x = 0
!       w(xk) = ONE + HALF * SQRT(MAX(-w(k1), ZERO))
!     else
!       ! Solve f3(x) = 0
!       w(k1) = ONE + ONEQUARTER * w(k1)
!       w(k0) = w(k0) * ONEQUARTER - ONEQUARTER
!       call find_a_cubic_root(w(k1), w(k0), w(xk), w(r), w(q), w(h), w(s))
!     end if
!   else
!     w(xk) = HALF * g / w(sa)
!     call newton_quartic(w(k1), w(k0), w(xk), w(xx), w(f), w(df))
!   end if
!   w(xk) = w(xk) * w(sa)
  end subroutine find_lambda_max
!
  pure elemental subroutine find_a_quartic_root(a, b, x, w1, w2, w3, w4, w5)
    real(RK), intent(in)    :: a, b
    real(RK), intent(inout) :: x, w1, w2, w3, w4, w5
    real(RK), parameter     :: ONESQRT3 = ONE / SQRT3
    real(RK), parameter     :: ONESIX = ONE / SIX
    real(RK), parameter     :: QUARTERSQRT3 = QUARTER * SQRT3
    real(RK), parameter     :: TWOTHIRD = TWO / THREE
    real(RK), parameter     :: CONST = 5.0_RK * SQRT3 / 36.0_RK
    integer                 :: i
    if (ABS(a) < THRESHOLD) then
      x = SQRT(ONE + SQRT(MAX(ONE - b, ZERO)))
      return
    end if
    call find_a_cubic_root_(-ONE, QUARTER * a, x, w1, w2, w3, w4)

    w5 = x * x
    w5 = w1 * (w1 - TWO) + a * x + b

    if (x < ONESQRT3) then
      call find_a_cubic_root_(QUARTERSQRT3 * a - TWOTHIRD, QUARTER * (a + SQRT3 * b) - CONST, x, w1, w2, w3, w4)
    else
      w1 = ONE / x                               ! 1/x0
      w2 = HALF * x - ONESIX * w1                ! g = x0/2 - (6 x0)**(-1)
      w3 = w2 * w2                               ! g**2
      call find_a_cubic_root_(-THREE * w3, QUARTER * (w5 * w1 + EIGHT * w2 * w3), x, w1, w2, w3, w4)
    end if
!
    if (THRESHOLD + w5 > ZERO) return
!
    do i = 1, MAXITER
      w1 = x * x
      w2 = -TWO + w1
      w3 = w2 * w1 + a * x + b
      w2 = (a + TWO * x * (w2 + w1))
      if (ABS(w3) < THRESHOLD) exit
      if (ABS(w2) < THRESHOLD) exit
      w1 = w1 / w2
      x = x - w1
      if (ABS(w1) < THRESHOLD) exit
    end do
!
  end subroutine find_a_quartic_root
!
! pure elemental subroutine newton_quartic(k1, k0, x, xx, f, df)
!   real(RK), intent(in)    :: k1, k0
!   real(RK), intent(inout) :: x, xx, f, df
!   integer                 :: i
!   do i = 1, MAXITER
!     xx = x * x
!     df = -TWO + xx
!     f = df * xx + k1 * x + k0
!     df = (k1 + (x + x) * (df + xx))
!     if (ABS(f) < THRESHOLD) exit
!     if (ABS(df) < THRESHOLD) exit
!     f = f / df
!     x = x - f
!     if (ABS(f) < THRESHOLD) exit
!   end do
! end subroutine newton_quartic
!
!| find a positive root of cubic equation x^3 + p x + q = 0, using Viete's solution.
  pure subroutine find_a_cubic_root_(p, q, x, w1, w2, w3, w4)
    real(RK), intent(in)    :: p, q
    real(RK), intent(inout) :: x, w1, w2, w3, w4
    if (ABS(p) < THRESHOLD) then
      if (ABS(q) < THRESHOLD) then
        x = ZERO ! X^3 = 0.
      else
        ! x^3 = -q
        call invcbrt(q, x, w1)
        x = -ONE / x
      end if
      return
    elseif (ABS(q) < THRESHOLD) then
      ! q = r**(-1/3)
      if (p < ZERO) then
        x = SQRT(-p)
      else
        x = SQRT(p)
      end if
    elseif (p < ZERO) then
      w1 = TWO * SQRT(p * (-ONETHIRD))
      w2 = -ONETHIRD * p * w1
      if (w2 < ABS(q)) then
        if (q < ZERO) then
          w2 = -w2 / q
          call cosh_acosh(w2, x, w3, w4) ! q = cos(arccos(h)/3)
          x = w1 * x
        else
          w2 = w2 / q
          call cosh_acosh(w2, x, w3, w4) ! q = cos(arccos(h)/3)
          x = -w1 * x
        end if
      else
        w2 = -q / w2
        call cos_acos(w2, x, w3, w4) ! q = cos(arccos(h)/3)
        x = w1 * x
      end if
    else
      w1 = TWO * SQRT(ONETHIRD * p)
      w2 = THREE * q / (p * w1)
      call sinh_asinh(w2, x, w3)
      x = -w1 * x
    end if
  end subroutine find_a_cubic_root_
!
!| Find null vector of S. <br>
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  pure subroutine find_null_vector(w)
    real(RK), intent(inout) :: w(*)
    integer(IK), parameter  :: y2 = 4, y3 = 5, y4 = 6
    integer(IK), parameter  :: s22 = 7, s23 = 8, s24 = 9
    integer(IK), parameter  :: s33 = 10, s34 = 11, s44 = 12
    integer(IK), parameter  :: m22 = 13, m23 = 14, m24 = 15
    integer(IK), parameter  :: m33 = 16, m34 = 17, m44 = 18
    integer(IK), parameter  :: mm = 4
!
    w(m22) = w(s33) * w(s44) - w(s34) * w(s34)
    w(m23) = w(s34) * w(s24) - w(s23) * w(s44)
    w(m24) = w(s23) * w(s34) - w(s33) * w(s24)
    w(mm) = ABS(w(m22)) + ABS(w(m23)) + ABS(w(m24))
    if (w(mm) < DEGENERACY) then
      w(m33) = w(s22) * w(s44) - w(s24) * w(s24)
      w(m34) = w(s22) * w(s34) - w(s23) * w(s24)
      w(mm) = ABS(w(m33)) + ABS(w(m34))
      if (w(mm) < DEGENERACY) then
        w(m44) = w(s22) * w(s33) - w(s23) * w(s23)
        w(mm) = ABS(w(m44))
        if (w(mm) < DEGENERACY) then
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
            ! Triple degeneracy
            w(y2) = ONE
            w(y3) = ZERO
            w(y4) = ZERO
          end if
        else
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
  end subroutine find_null_vector
!
  pure subroutine det3(g, cov, w)
    real(RK), intent(in)    :: g, cov(*)
    real(RK), intent(inout) :: w(*)
    w(1) = ABS(cov(1))
    w(2) = ABS(cov(2))
    w(3) = g * THRESHOLD
    if (w(1) > w(2)) then
      w(2) = ABS(cov(3))
      if (w(1) > w(2)) then
! |cov(1)| > |cov(2)| .and. |cov(1)| > |cov(3)|
        if (w(1) < w(3)) then
! |cov(2)|, |cov(3)|, |cov(1)| < threshold
          w(1) = ZERO
          return
        end if
        w(1) = ABS(cov(5))
        w(2) = ABS(cov(6))
        if (w(1) > w(2)) then
!   147
!   258
!   369, pivot = 0
          if (w(1) < w(3)) then
            w(1) = cov(4) * (cov(2) * cov(9) - cov(8) * cov(3))
          else
            w(1) = ((cov(1) * cov(5) - cov(2) * cov(4)) &
                & * (cov(1) * cov(9) - cov(3) * cov(7)) &
                & - (cov(1) * cov(8) - cov(2) * cov(7)) &
                & * (cov(1) * cov(6) - cov(3) * cov(4)) &
                   ) / cov(1)
          end if
        else
!   147
!   369
!   258, pivot = 1
          if (w(2) < w(3)) then
            w(1) = -cov(4) * (cov(3) * cov(8) - cov(9) * cov(2))
            return
          else
            w(1) = -((cov(1) * cov(6) - cov(3) * cov(4)) &
                & * (cov(1) * cov(8) - cov(2) * cov(7)) &
                & - (cov(1) * cov(9) - cov(3) * cov(7)) &
                & * (cov(1) * cov(5) - cov(2) * cov(4)) &
                   ) / cov(1)
          end if
        end if
      else
! |cov(3)| > |cov(1)| > |cov(2)|
        if (w(2) < w(3)) then
          w(1) = ZERO
          return
        end if
        w(1) = ABS(cov(4))
        w(2) = ABS(cov(5))
        if (w(1) > w(2)) then
!   369
!   147
!   258, pivot = 2
          if (w(1) < w(3)) then
            w(1) = cov(6) * (cov(1) * cov(8) - cov(7) * cov(2))
          else
            w(1) = ((cov(3) * cov(4) - cov(1) * cov(6)) &
                & * (cov(3) * cov(8) - cov(2) * cov(9)) &
                & - (cov(3) * cov(7) - cov(1) * cov(9)) &
                & * (cov(3) * cov(5) - cov(2) * cov(6)) &
                   ) / cov(3)
          end if
        else
!   369
!   258
!   147, pivot = 1
          if (w(2) < w(3)) then
            w(1) = cov(4) * (cov(2) * cov(9) - cov(8) * cov(3))
          else
            w(1) = -((cov(3) * cov(5) - cov(2) * cov(6)) &
                & * (cov(3) * cov(7) - cov(1) * cov(9)) &
                & - (cov(3) * cov(8) - cov(2) * cov(9)) &
                & * (cov(3) * cov(4) - cov(1) * cov(6)) &
                    ) / cov(3)
          end if
        end if
      end if
    else
      w(1) = ABS(cov(3))
      if (w(1) > w(2)) then
! |cov(3)| > |cov(2)| >= |cov(1)|
        if (w(1) < w(3)) then
          w(1) = ZERO
          return
        end if
        w(1) = ABS(cov(4))
        w(2) = ABS(cov(5))
        if (w(1) > w(2)) then
!   369
!   147
!   258, pivot = 2
          if (w(1) < w(3)) then
            w(1) = cov(6) * (cov(1) * cov(8) - cov(7) * cov(2))
          else
            w(1) = ((cov(3) * cov(4) - cov(1) * cov(6)) &
                & * (cov(3) * cov(8) - cov(2) * cov(9)) &
                & - (cov(3) * cov(7) - cov(1) * cov(9)) &
                & * (cov(3) * cov(5) - cov(2) * cov(6)) &
                   ) / cov(3)
          end if
        else
!   369
!   258
!   147, pivot = 1
          if (w(2) < w(3)) then
            w(1) = cov(4) * (cov(2) * cov(9) - cov(8) * cov(3))
          else
            w(1) = -((cov(3) * cov(5) - cov(2) * cov(6)) &
                & * (cov(3) * cov(7) - cov(1) * cov(9)) &
                & - (cov(3) * cov(8) - cov(2) * cov(9)) &
                & * (cov(3) * cov(4) - cov(1) * cov(6)) &
                    ) / cov(3)
          end if
        end if
      else
! |cov(2)| >= |cov(1)| .and. |cov(2)| >= |cov(3)|
        if (w(2) < w(3)) then
          w(1) = ZERO
          return
        end if
        w(1) = ABS(cov(4))
        w(2) = ABS(cov(6))
        if (w(1) > w(2)) then
!   258
!   147
!   369, pivot = 1
          if (w(1) < w(3)) then
            w(1) = -cov(5) * (cov(1) * cov(9) - cov(7) * cov(3))
          else
            w(1) = -((cov(2) * cov(4) - cov(1) * cov(5)) &
                & * (cov(2) * cov(9) - cov(3) * cov(8)) &
                & - (cov(2) * cov(7) - cov(1) * cov(8)) &
                & * (cov(2) * cov(6) - cov(3) * cov(5)) &
                   ) / cov(2)
          end if
        else
!   258
!   369
!   147, pivot = 2
          if (w(2) < w(3)) then
            w(1) = cov(5) * (cov(3) * cov(7) - cov(9) * cov(1))
          else
            w(1) = ((cov(2) * cov(6) - cov(3) * cov(5)) &
                & * (cov(2) * cov(7) - cov(1) * cov(8)) &
                & - (cov(2) * cov(9) - cov(3) * cov(8)) &
                & * (cov(2) * cov(4) - cov(1) * cov(5)) &
                  ) / cov(2)
          end if
        end if
      end if
    end if
  end subroutine det3
!
!| find a positive root of monic cubic equation, x^3 - 2 * x^2 + k1 * x + k0 = 0
! pure subroutine find_a_cubic_root(k1, k0, x, r, q, h, s)
!   real(RK), intent(in)    :: k1, k0
!   real(RK), intent(inout) :: x, r, q, h, s
!   R = -FOUR * ((k1 - EIGHTNINE) * TWOTHIRD + k0)
!   Q = FOURTHIRD * (FOURTHIRD - k1)
!   if (ABS(Q) < THRESHOLD) then
!     if (ABS(R) < THRESHOLD) then
!       x = TWOTHIRD
!     else
!       r = ONEQUARTER * r
!       call invcbrt(r, q, s) ! q = r**(-1/3)
!       x = TWOTHIRD + q
!     end if
!     return
!   elseif (ABS(r) < THRESHOLD) then
!     if (q > ZERO) then
!       x = TWOTHIRD + HALFSQRT3 * SQRT(q)
!     else
!       x = TWOTHIRD - HALFSQRT3 * SQRT(-q)
!     end if
!   elseif (Q > ZERO) then
!     s = SQRT(q)
!     q = s * q
!     if (ABS(R) <= Q) then
!       h = r / q
!       call cos_acos(h, q, r, x) ! q = cos(arccos(h)/3)
!       X = TWOTHIRD + s * q
!     else
!       h = ABS(q / r)
!       s = SIGN(ONE, r) * s
!       call cosh_acosh(H, q, r, x) ! q = cosh(arccosh(1/h)/3)
!       X = TWOTHIRD + s * q
!     end if
!   else
!     s = SQRT(-q)
!     h = -r / (q * s)
!     call sinh_asinh(h, q, r)
!     x = x + s * q
!   end if
! end subroutine find_a_cubic_root
!
#include "cos_acos.f90"
#include "cosh_acosh.f90"
#include "sinh_asinh.f90"
!
!| https://www.mdpi.com/1996-1073/14/4/1058
  pure elemental subroutine invcbrt(x, res, c)
    use mod_kinds, only: I4, R4
    implicit none
    real(RK), intent(in)    :: x
    real(RK), intent(inout) :: res, c
    real(R4)                :: y
    y = ABS(x)
    res = SIGN(ONE, x) * TRANSFER(INT(z"548C2B4B", I4) - TRANSFER(y, 0_I4) / 3, y)
    c = res * res * res * x
    res = res * (1.752319676_RK - c * (1.2509524245_RK - 0.5093818292_RK * c))
    c = ONE - res * res * res * x
    res = res * (ONE + ONETHIRD * c)
  end subroutine invcbrt
!
end module mod_rotation


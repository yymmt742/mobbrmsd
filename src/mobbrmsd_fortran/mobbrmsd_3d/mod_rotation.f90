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
  real(RK), parameter    :: SQRT2 = SQRT(2.0_RK)
  real(RK), parameter    :: SQRT3 = SQRT(3.0_RK)
  real(RK), parameter    :: HALFSQRT3 = HALF * SQRT3
#ifdef USE_REAL32
  real(RK), parameter    :: THRESHOLD = 1E-6_RK
  real(RK), parameter    :: DEGENERACY = 1E-3_RK
#else
  real(RK), parameter    :: THRESHOLD = 1E-12_RK
  real(RK), parameter    :: DEGENERACY = 1E-6_RK
#endif
  integer(IK), parameter :: MAXITER = 10
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
    !pure subroutine estimate_rcmax(g, cov, w)
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
    !pure subroutine estimate_sdmin(g, cov, w)
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
    res = 17
  end function rotation_worksize
!
!| Compute the transpose rotation matrix for minimize tr[CR] from cov = YX^T and g = tr[XX^T] + tr[YY^T].
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  pure subroutine estimate_rotation(g, cov, rot, w)
    !pure subroutine estimate_rotation(g, cov, rot, w)
    real(RK), intent(in)    :: g
    !! g = tr[XX^T] + tr[YY^T]
    real(RK), intent(in)    :: cov(*)
    !! covariance dxd matrix, YX^T
    real(RK), intent(inout) :: rot(*)
    !! rotation dxd matrix
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_rotation().
    integer(IK), parameter  :: q0 = 1, q1 = 2, q2 = 3, q3 = 4, ww = 15
    integer(IK), parameter  :: a00 = 5
    integer(IK), parameter  :: a10 = 6, a11 = 7
    integer(IK), parameter  :: a20 = 8, a21 = 9, a22 = 10
    integer(IK), parameter  :: a30 = 11, a31 = 12, a32 = 13, a33 = 14
    integer(IK), parameter  :: p01 = 1, m01 = 2, p2d = 3, m2d = 4
    integer(IK), parameter  :: abs00 = 1, abs11 = 2, abs22 = 2, abs33 = 1
    integer(IK), parameter  :: v11 = 5, v21 = 6, v31 = 7, v41 = 8
    integer(IK), parameter  :: v22 = 9, v32 = 10, v42 = 11, v33 = 12, v43 = 13, v44 = 14
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
!   L = (R11+R22+R33  R23-R32      R31-R13      R12-R21    )
!       (R23-R32      R11-R22-R33  R12+R21      R13+R31    )
!       (R31-R13      R12+R21     -R11+R22-R33  R23+R32    )
!       (R12-R21      R13+R31      R23+R32     -R11-R22+R33)
!
    w(p2d) = cov(9) - w(1)
    w(m2d) = -cov(9) - w(1)
    w(p01) = cov(1) + cov(5)
    w(m01) = cov(1) - cov(5)
    w(a00) = w(p2d) + w(p01)
    w(a11) = w(m2d) + w(m01)
    w(a22) = w(m2d) - w(m01)
    w(a33) = w(p2d) - w(p01)
    w(a10) = cov(8) - cov(6)
    w(a20) = cov(3) - cov(7); w(a21) = cov(2) + cov(4)
    w(a30) = cov(4) - cov(2); w(a31) = cov(7) + cov(3); w(a32) = cov(6) + cov(8)
!
    w(abs00) = ABS(w(a00))
    w(abs22) = ABS(w(a22))

    if (w(abs00) > w(abs22)) then
      w(abs11) = ABS(w(a11))
      if (w(abs00) > w(abs11)) then
        call find_null_vector(w(a00), w(a10), w(a11), w(a20), w(a21), w(a22), w(a30), w(a31), w(a32), w(a33), &
            &                 w(q0), w(q1), w(q2), w(q3), w(ww))
      else
        call find_null_vector(w(a11), w(a10), w(a00), w(a21), w(a20), w(a22), w(a31), w(a30), w(a32), w(a33), &
            &                 w(q1), w(q0), w(q2), w(q3), w(ww))
      end if
    else
      w(abs33) = ABS(a33)
      if (w(abs22) > w(abs33)) then
        call find_null_vector(w(a22), w(a21), w(a11), w(a20), w(a10), w(a00), w(a32), w(a31), w(a30), w(a33), &
            &                 w(q2), w(q1), w(q0), w(q3), w(ww))
      else
        call find_null_vector(w(a33), w(a31), w(a11), w(a32), w(a21), w(a22), w(a30), w(a10), w(a20), w(a00), &
            &                 w(q3), w(q1), w(q2), w(q0), w(ww))
      end if
    end if
!
    w(v11) = w(q0) * w(q0)
    w(v22) = w(q1) * w(q1)
    w(v33) = w(q2) * w(q2)
    w(v44) = w(q3) * w(q3)
    w(q0) = w(q0) * SQRT2
    w(q1) = w(q1) * SQRT2
    w(q2) = w(q2) * SQRT2
    w(q3) = w(q3) * SQRT2
    w(v21) = w(q0) * w(q1)
    w(v31) = w(q0) * w(q2)
    w(v41) = w(q0) * w(q3)
    w(v32) = w(q1) * w(q2)
    w(v42) = w(q1) * w(q3)
    w(v43) = w(q2) * w(q3)
!
    rot(1) = w(v11) + w(v22) - w(v33) - w(v44)
    rot(2) = w(v32) - w(v41)
    rot(3) = w(v42) + w(v31)
    rot(4) = w(v32) + w(v41)
    rot(5) = w(v11) - w(v22) + w(v33) - w(v44)
    rot(6) = w(v43) - w(v21)
    rot(7) = w(v42) - w(v31)
    rot(8) = w(v43) + w(v21)
    rot(9) = w(v11) - w(v22) - w(v33) + w(v44)
  end subroutine estimate_rotation
!
!| Compute maximum eigen value of S. <br>
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  pure subroutine find_lambda_max(g, cov, w)
    !pure subroutine find_lambda_max(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    integer(IK), parameter  :: x = 1, a = 2, b = 3, sr = 4
    integer(IK), parameter  :: d11 = 4, d22 = 5, d33 = 6, d21 = 7, d31 = 8, d32 = 9
    integer(IK), parameter  :: p = 3, r = 2, q = 1, rr = 5
    integer(IK), parameter  :: w1 = 5, w2 = 6, w3 = 7, w4 = 8, w5 = 9
!
    if (g < THRESHOLD) then
      w(x) = ZERO
      return
    end if
!
!   K1 = - 8 det|R|
!   A2 = tr[D]
!    call det3(g, cov, w(q))
    call det3(cov(1), cov(2), cov(3), cov(4), cov(5), cov(6), cov(7), cov(8), cov(9), w(q))
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
    w(p) = w(d22) + w(d33)
    w(r) = w(d11) + w(p)
    if (w(r) < THRESHOLD) then
      w(x) = ZERO
      return
    end if
    w(p) = w(p) * w(d11) + w(d22) * w(d33) - w(d21) * w(d21) - w(d31) * w(d31) - w(d32) * w(d32)
    w(sr) = SQRT(w(r))
    w(rr) = FOUR / (w(r) * w(r))
    w(a) = -TWO * w(q) * w(rr) * w(sr)
    w(b) = ONE - w(p) * w(rr)
    call find_a_quartic_root(w(a), w(b), w(x), w(w1), w(w2), w(w3), w(w4), w(w5))
    w(x) = w(x) * w(sr)
!
  end subroutine find_lambda_max
!
  pure elemental subroutine find_a_quartic_root(a, b, x, w1, w2, w3, w4, w5)
    real(RK), intent(in)    :: a, b
    real(RK), intent(inout) :: x, w1, w2, w3, w4, w5
    real(RK), parameter     :: ONESQRT3 = ONE / SQRT3
    real(RK), parameter     :: ONESIX = ONE / SIX
    real(RK), parameter     :: QUARTERSQRT3 = QUARTER * SQRT3
    real(RK), parameter     :: TWOTHIRD = TWO / THREE
    real(RK), parameter     :: CONST = -5.0_RK * SQRT3 / 36.0_RK + ONESQRT3
    integer                 :: i
    if (ABS(a) < THRESHOLD) then
      x = SQRT(ONE + SQRT(MAX(ONE - b, ZERO)))
      return
    end if
    call find_a_cubic_root(-ONE, QUARTER * a, x, w1, w2, w3, w4)

    w5 = x * x
    w5 = w5 * (w5 - TWO) + a * x + b
    i = MERGE(1, 0, THRESHOLD + w5 > ZERO)

    if (x < ONESQRT3) then
      call find_a_cubic_root(QUARTERSQRT3 * a - TWOTHIRD, QUARTER * (a + SQRT3 * b) + CONST, x, w1, w2, w3, w4)
    else
      w1 = ONE / x                ! 1/x0
      w2 = HALF * x - ONESIX * w1 ! g = x0/2 - (6 x0)**(-1)
      w3 = w2 * w2                ! g**2
      w4 = w5 * w1                ! x - g
      w5 = x - w2                 ! x - g
      call find_a_cubic_root(-THREE * w3, QUARTER * (w4 + EIGHT * w2 * w3), x, w1, w2, w3, w4)
      x = x + w5
    end if
!
    if (i > 0) return
!
    do i = 1, MAXITER
      w1 = x * x
      w2 = -TWO + w1
      w3 = w2 * w1 + a * x + b
      if (ABS(w3) < THRESHOLD) exit
      w2 = a + TWO * x * (w2 + w1)
      if (ABS(w2) < THRESHOLD) exit
      w1 = w3 / w2
      x = x - w1
      if (ABS(w1) < THRESHOLD) exit
    end do
!
  end subroutine find_a_quartic_root
!
!| find a positive root of cubic equation x^3 + p x + q = 0, using Viete's solution.
  pure subroutine find_a_cubic_root(p, q, x, w1, w2, w3, w4)
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
  end subroutine find_a_cubic_root
!
!| Find null vector of S. <br>
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  pure subroutine find_null_vector(a00, a10, a11, a20, a21, a22, a30, a31, a32, a33, q0, q1, q2, q3, w)
    real(RK), intent(inout) :: a00
    real(RK), intent(inout) :: a10, a11
    real(RK), intent(inout) :: a20, a21, a22
    real(RK), intent(inout) :: a30, a31, a32, a33
    real(RK), intent(inout) :: q0, q1, q2, q3, w(*)
    integer(IK), parameter  :: d3 = 1
    integer(IK), parameter  :: l1 = 1, l2 = 2, l3 = 3
!
    if (ABS(a00) < THRESHOLD) then
      a00 = ZERO
      call det3(a11, a21, a31, a21, a22, a32, a31, a32, a33, w(d3)); q0 = w(d3); a00 = a00 + w(d3) * w(d3)
      call det3(a10, a21, a31, a30, a32, a33, a20, a22, a32, w(d3)); q1 = w(d3); a00 = a00 + w(d3) * w(d3)
      call det3(a10, a11, a31, a20, a21, a32, a30, a31, a33, w(d3)); q2 = w(d3); a00 = a00 + w(d3) * w(d3)
      call det3(a10, a11, a21, a30, a31, a32, a20, a21, a22, w(d3)); q3 = w(d3); a00 = a00 + w(d3) * w(d3)
      a00 = ONE / SQRT(a00)
      q0 = q0 * a00
      q1 = q1 * a00
      q2 = q2 * a00
      q3 = q3 * a00
      return
    end if
!
    a00 = ONE / a00
    w(l1) = -a10 * a00
    w(l2) = -a20 * a00
    w(l3) = -a30 * a00
    a11 = a11 + a10 * w(l1)
    a21 = a21 + a20 * w(l1); a22 = a22 + a20 * w(l2)
    a31 = a31 + a30 * w(l1); a32 = a32 + a30 * w(l2); a33 = a33 + a30 * w(l3)
    q1 = a22 * a33 - a32 * a32
    q2 = a31 * a32 - a21 * a33
    q3 = a21 * a32 - a31 * a22
    a00 = ABS(q1) + ABS(q2) + ABS(q3)

    if (a00 < DEGENERACY) then
      q1 = ZERO
      q2 = a11 * a33 - a31 * a31
      q3 = a21 * a31 - a11 * a32
      a00 = ABS(q2) + ABS(q3)
      if (a00 < DEGENERACY) then
        a00 = ABS(a11 * a22 - a21 * a21)
        if (a00 < DEGENERACY) then
          ! double degeneracy
          if (ABS(a11) > DEGENERACY) then
            q0 = w(l1) * a21 - w(l2) * a11
            q1 = a21
            q2 = -a11
            q3 = ZERO
            a00 = ONE / SQRT(q0 * q0 + q1 * q1 + q2 * q2)
            q0 = q0 * a00
            q1 = q1 * a00
            q2 = q2 * a00
          elseif (ABS(a22) > DEGENERACY) then
            q0 = w(l2) * a32 - w(l3) * a22
            q1 = ZERO
            q2 = a32
            q3 = -a22
            a00 = ONE / SQRT(q0 * q0 + q2 * q2 + q3 * q3)
            q0 = q0 * a00
            q2 = q2 * a00
            q3 = q3 * a00
          elseif (ABS(a33) > DEGENERACY) then
            q0 = w(l1) * a33 - w(l3) * a31
            q1 = a33
            q2 = ZERO
            q3 = -a31
            a00 = ONE / SQRT(q0 * q0 + q1 * q1 + q3 * q3)
            q0 = q0 * a00
            q1 = q1 * a00
            q3 = q3 * a00
          else
!           ! Triple degeneracy
            q0 = ONE
            q1 = ZERO
            q2 = ZERO
            q3 = ZERO
          end if
        else
          q1 = ZERO
          q2 = ZERO
          q3 = ONE / (one + w(l3) * w(l3))
          q0 = w(l3) * q3
        end if
      else
        q0 = w(l2) * q2 + w(l3) * q3
        a00 = ONE / SQRT(q0 * q0 + q2 * q2 + q3 * q3)
        q0 = q0 * a00
        q2 = q2 * a00
        q3 = q3 * a00
      end if
    else
      q0 = w(l1) * q1 + w(l2) * q2 + w(l3) * q3
      a00 = ONE / SQRT(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
      q0 = q0 * a00
      q1 = q1 * a00
      q2 = q2 * a00
      q3 = q3 * a00
    end if
  end subroutine find_null_vector
!
  pure subroutine det3(a11, a12, a13, a21, a22, a23, a31, a32, a33, w)
    real(RK), intent(in)    :: a11, a12, a13
    real(RK), intent(in)    :: a21, a22, a23
    real(RK), intent(in)    :: a31, a32, a33
    real(RK), intent(inout) :: w(*)
    w(1) = ABS(a11)
    w(2) = ABS(a21)
    w(3) = ABS(a31)
    if (w(1) > w(2)) then
      if (w(1) > w(2)) then
! |cov(1)| > |cov(2)| .and. |cov(1)| > |cov(3)|
        if (w(1) < THRESHOLD) then ! |a31|, |a21| < |a11| < threshold
          w(1) = ZERO; return
        end if
        w(1) = ABS(a22)
        w(2) = ABS(a32)
        if (w(1) > w(2)) then
!           1 4 7
!           2 5 8
!           3 6 9, pivot = 0
          call det3p(a11, a12, a13, a21, a22, a23, a31, a32, a33, w(1)); return
        else
!           1 4 7
!           3 6 9
!           2 5 8, pivot = 1
          call det3m(a11, a12, a13, a31, a32, a33, a21, a22, a23, w(1)); return
        end if
      end if
    else
      if (w(3) < w(2)) then ! |a11| <= |a21| and |a31| < |a21|
        w(1) = ZERO
        if (w(2) < THRESHOLD) then ! |a11|, |a31| < |a21| < threshold
          w(1) = ZERO
          return
        end if
        w(1) = ABS(a12)
        w(2) = ABS(a32)
        if (w(1) > w(2)) then
!             2 5 8
!             1 4 7
!             3 6 9, pivot = 1
          call det3m(a21, a22, a23, a11, a12, a13, a31, a32, a33, w(1)); return
        else
!             2 5 8
!             3 6 9
!             1 4 7, pivot = 2
          call det3p(a21, a22, a23, a31, a32, a33, a11, a12, a13, w(1)); return
        end if
      end if
    end if

    if (w(3) < THRESHOLD) then ! |a11|, |a21| <= |a31| < threshold
      w(1) = ZERO; return
    end if
    w(1) = ABS(a12)
    w(2) = ABS(a22)
    if (w(1) < w(2)) then
!       3 6 9
!       1 4 7
!       2 5 8, pivot = 2
      call det3p(a31, a32, a33, a11, a12, a13, a21, a22, a23, w(1)); return
    else
!       3 6 9
!       2 5 8
!       1 4 7, pivot = 1
      call det3m(a31, a32, a33, a21, a22, a23, a11, a12, a13, w(1)); return
    end if
  end subroutine det3
!
  pure elemental subroutine det3p(a11, a12, a13, a21, a22, a23, a31, a32, a33, ret)
    real(RK), intent(in)    :: a11, a12, a13
    real(RK), intent(in)    :: a21, a22, a23
    real(RK), intent(in)    :: a31, a32, a33
    real(RK), intent(inout) :: ret
    ret = ((a11 * a22 - a21 * a12) * (a11 * a33 - a31 * a13) &
   &     - (a11 * a23 - a21 * a13) * (a11 * a32 - a31 * a12)) &
   &     / a11
  end subroutine det3p
!
  pure elemental subroutine det3m(a11, a12, a13, a21, a22, a23, a31, a32, a33, ret)
    real(RK), intent(in)    :: a11, a12, a13
    real(RK), intent(in)    :: a21, a22, a23
    real(RK), intent(in)    :: a31, a32, a33
    real(RK), intent(inout) :: ret
    ret = -((a11 * a22 - a21 * a12) * (a11 * a33 - a31 * a13) &
   &      - (a11 * a23 - a21 * a13) * (a11 * a32 - a31 * a12)) &
   &      / a11
  end subroutine det3m
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


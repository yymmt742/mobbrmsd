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
  real(RK), parameter    :: HALF = 0.5_RK
  real(RK), parameter    :: ONE = 1.0_RK
  real(RK), parameter    :: TWO = 2.0_RK
  real(RK), parameter    :: THREE = 3.0_RK
  real(RK), parameter    :: FOUR = 4.0_RK
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
  real(RK), parameter    :: DEGENERACY1 = -1E-4_RK
  real(RK), parameter    :: DEGENERACY2 = 1E-12_RK
#else
  real(RK), parameter    :: THRESHOLD = 1E-12_RK
  real(RK), parameter    :: DEGENERACY = 1E-4_RK
  real(RK), parameter    :: DEGENERACY1 = -1E-4_RK
  real(RK), parameter    :: DEGENERACY2 = 1E-16_RK
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
  subroutine estimate_rcmax(g, cov, w)
    !pure subroutine estimate_rcmax(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    !! work array, must be larger than worksize_sdmin().
    call find_lambda_max(g, cov, w)
    !w(1) = HALF * g * w(1)
  end subroutine estimate_rcmax
!
!| Compute the least-squares sum_i^n |x_i-Ry_i|^2 from cov = YX^T and g = tr[XX^T] + tr[YY^T].
  subroutine estimate_sdmin(g, cov, w)
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
    res = 18
  end function rotation_worksize
!
!| Compute the transpose rotation matrix for minimize tr[CR] from cov = YX^T and g = tr[XX^T] + tr[YY^T].
!  This subroutine is based on the method of Coutsias et.al. 10.1002/jcc.25802
  subroutine estimate_rotation(g, cov, rot, w)
    !pure subroutine estimate_rotation(g, cov, rot, w)
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
!   w(nrm) = TWO / g
!   w(dg4) = w(nrm) * cov(9) - w(1)
!   w(dg3) = w(nrm) * cov(9) + w(1)
!   w(dg2) = w(nrm) * (cov(1) - cov(5))
!   w(dg1) = w(nrm) * (cov(1) + cov(5))
!   w(a11) = w(dg1) + w(dg4)
!   w(a21) = w(nrm) * (cov(8) - cov(6))
!   w(a31) = w(nrm) * (cov(3) - cov(7))
!   w(a41) = w(nrm) * (cov(4) - cov(2))
!   w(a22) = w(dg2) - w(dg3)
!   w(a32) = w(nrm) * (cov(4) + cov(2))
!   w(a42) = w(nrm) * (cov(7) + cov(3))
!   w(a33) = -w(dg2) - w(dg3)
!   w(a43) = w(nrm) * (cov(8) + cov(6))
!   w(a44) = -w(dg1) + w(dg4)
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
  subroutine find_lambda_max(g, cov, w)
    !pure subroutine find_lambda_max(g, cov, w)
    real(RK), intent(in)    :: g
    !! sum of auto covariance matrix
    real(RK), intent(in)    :: cov(*)
    !! target d*n array
    real(RK), intent(inout) :: w(*)
    integer(IK), parameter  :: xk = 1, k1 = 2, k0 = 3, nrm1 = 5, nrm2 = 6
    integer(IK), parameter  :: d11 = 3, d22 = 4, d33 = 5, d21 = 6, d31 = 7, d32 = 8
    integer(IK), parameter  :: a1 = 6, a2 = 5
    integer(IK), parameter  :: r = 5, q = 6, h = 7, s = 8
    integer(IK), parameter  :: xx = 5, f = 6, df = 7
    real(RK)                :: ra2, sa2, gsa2, minima, c0, c1
    integer(IK)             :: k
!
    if (g < THRESHOLD) then
      w(1) = ZERO
      return
    end if
!
!   K1 = - 8 det|R|
    call det3(g, cov, w)
    w(k1) = -EIGHT * w(1)
!   w(k1) = w(1)
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
    w(a1) = w(d11) * (w(d22) + w(d33)) + w(d22) * w(d33) - w(d21)**2 - w(d31)**2 - w(d32)**2
    w(a1) = FOUR * w(a1)
!
!   A2 = tr[D]
!
    w(a2) = w(d11) + w(d22) + w(d33)
!
!   K0 = 16 * (A2**2 - 4*A1)
!   K2 = -8 * A2
!
    sa2 = SQRT(w(a2))
    gsa2 = g / sa2
!
    w(k1) = w(k1) / (sa2 * w(a2))
    w(k0) = ONE - w(a1) / (w(a2) * w(a2))
!
!   normalize
!
    minima = w(k1) + w(k0) - ONE
    if (minima >= ZERO) then
      print *, "minima >= ZERO"
      ! multiple root
      w(xk) = ONE
    elseif (minima > DEGENERACY1) then
      ! quasi biquadratic equations.
      c0 = w(k0) * ONEQUARTER - ONEQUARTER
      if (ABS(c0) < DEGENERACY2) then
        print *, "ABS(c0) < 1.E-14"
        ! Second order Taylor expansion around the maximum local minima (x=1).
        ! f2(x) = x**2 - (2 - 1/4*k1) * x + (k0 + 3) / 4
        w(xk) = ONE - (w(k1) - SQRT(w(k1) * w(k1) - minima * SIXTEEN)) * ONEEIGHT
      else
        ! Third order Taylor expansion around the maximum local minima (x=1).
        c1 = ONE - ONEQUARTER * w(k1)
        call find_a_cubic_root(c1, c0, w(xk), w(r), w(q), w(h), w(s))
        print *, "else ABS(c0) < 1.E-8", minima, c0, c1, w(xk)**3 - 2._RK * w(xk) + c1 * w(xk) + c0
      end if
    else
      w(xk) = gsa2
      call newton(w(k1), w(k0), w(xk), w(xx), w(f), w(df))
    end if
!
    w(xk) = w(xk) * sa2
  end subroutine find_lambda_max
!
  pure elemental subroutine newton(k1, k0, x, xx, f, df)
    real(RK), intent(in)    :: k1, k0
    real(RK), intent(inout) :: x, xx, f, df
    integer                 :: i
    do i = 1, MAXITER
      xx = x * x
      df = -TWO + xx
      f = df * xx + k1 * x + k0
      df = (k1 + (x + x) * (df + xx))
      if (ABS(f) < THRESHOLD) exit
      if (ABS(df) < THRESHOLD) exit
      f = f / df
      x = x - f
      if (ABS(f) < THRESHOLD) exit
    end do
  end subroutine newton
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
        if (w(1) < w(2)) then
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
!   369
!   147
!   258, pivot = 2
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
            !print*,3
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
! |cov(3)| > |cov(2)| => |cov(1)|
!   369
!   147
!   258, pivot = 2
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
! |cov(2)| => |cov(1)| .and. |cov(2)| => |cov(3)|
!   258
!   147
!   369, pivot = 1
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
! find a positive root of monic cubic equation, x^3 - 2 * x^2 + k1 * x + k0 = 0
  subroutine find_a_cubic_root(k1, k0, x, r, q, h, s)
    real(RK), intent(in)    :: k1, k0
    real(RK), intent(inout) :: x, r, q, h, s
    X = TWOTHIRD
    R = -FOUR * (-(EIGHTNINE - k1) * TWOTHIRD + k0)
    Q = -FOURTHIRD * (k1 - FOURTHIRD)
    if (ABS(Q) < THRESHOLD) then
      if (ABS(R) < THRESHOLD) return
      X = X + SIGN(ONE, R) * (ONEQUARTER * R)**(-ONETHIRD)
    elseif (ABS(R) < THRESHOLD) then
      if (q > ZERO) then
        X = X + HALFSQRT3 * SQRT(q)
      else
        X = X - HALFSQRT3 * SQRT(-q)
      end if
    elseif (Q > ZERO) then
      S = SQRT(Q)
      Q = S * Q
      if (ABS(R) <= Q) then
        H = R / Q
        X = X + S * cos_acos(H)
      else
        H = ABS(Q / R)
        X = X + SIGN(ONE, R) * S * cosh_acosh(H)
      end if
    else
      S = SQRT(-Q)
      H = -R / (Q * S)
      X = X + S * sinh_asinh(H)
    end if
  end subroutine find_a_cubic_root
!
  function cos_acos(x) result(res)
    real(RK), intent(in) :: x
    real(RK)             :: res, yy, df
    integer              :: i
    if (x < -1.0000000000000000E+00_RK) then
      res = 1.0000000000000000E+00_RK / 2.0000000000000000E+00_RK
      return
    elseif (x < -0.3333333333333333E+00_RK) then
      yy = SQRT(x + ONE)
      if (x < -6.6666666666666663E-01_RK) then
        if (x < -9.0236892706218241E-01_RK) then
          res = (3.9847929621330180E-04 * yy + (2.7150043732313152E-03_RK)) * x + (2.7756772944468894E-03 * yy +&
              & (8.0584699091552827E-03_RK))
          res = res * x + (9.9787382004175255E-03 * yy + (-3.0422013335204556E-04_RK))
          res = res * x + (3.2124745520119027E-02 * yy + (-6.9480186781926440E-02_RK))
          res = res * x + (4.3277149578179819E-01 * yy + (4.3616749888734951E-01_RK))
        else
          res = (3.9847929621330180E-04 * yy + (-8.1880793225013100E-04_RK)) * x + (2.7756772944468894E-03 * yy +&
              & (-5.3682021468007667E-03_RK))
          res = res * x + (9.9787382004175255E-03 * yy + (-1.9427315377984092E-02_RK))
          res = res * x + (3.2124745520119027E-02 * yy + (-8.1580601406840383E-02_RK))
          res = res * x + (4.3277149578179819E-01 * yy + (4.3329732093336737E-01_RK))
        end if
      else
        if (x < -4.3096440627115074E-01_RK) then
          res = (3.9847929621330180E-04 * yy + (-1.3962794333934664E-03_RK)) * x + (2.7756772944468894E-03 * yy +&
              & (-6.7281227183138212E-03_RK))
          res = res * x + (9.9787382004175255E-03 * yy + (-2.0640121298405912E-02_RK))
          res = res * x + (3.2124745520119027E-02 * yy + (-8.2067600037457583E-02_RK))
          res = res * x + (4.3277149578179819E-01 * yy + (4.3322280904921723E-01_RK))
        else
          res = (3.9847929621330180E-04 * yy + (1.1713621213896104E-02_RK)) * x + (2.7756772944468894E-03 * yy +&
              & (1.3560409301341475E-02_RK))
          res = res * x + (9.9787382004175255E-03 * yy + (-8.8387594483930760E-03_RK))
          res = res * x + (3.2124745520119027E-02 * yy + (-7.9004467646912643E-02_RK))
          res = res * x + (4.3277149578179819E-01 * yy + (4.3352276164963782E-01_RK))
        end if
      end if
    elseif (x < 0.3333333333333333E+00) then
      if (x < 2.0410779985789219E-17_RK) then
        if (x < -2.3570226039551581E-01_RK) then
          res = (6.3261501181342716E-01_RK) * x + (7.7758615573901924E-01_RK)
          res = res * x + (2.7435256739631902E-01_RK)
          res = res * x + (2.2752170341132330E-01_RK)
          res = res * x + (8.7030281603695292E-01_RK)
        else
          res = (-5.1407513495624758E-02_RK) * x + (1.0016161738288792E-02_RK)
          res = res * x + (-5.0149303892811525E-02_RK)
          res = res * x + (1.6657415217801974E-01_RK)
          res = res * x + (8.6602540378443871E-01_RK)
        end if
      else
        if (x < 2.3570226039551587E-01_RK) then
          res = (9.3008662532715076E-03_RK) * x + (1.4372359674536180E-02_RK)
          res = res * x + (-4.6656484877630987E-02_RK)
          res = res * x + (1.6659959026260604E-01_RK)
          res = res * x + (8.6602540378443871E-01_RK)
        else
          res = (-3.4902175486090642E-01_RK) * x + (4.1033381483828746E-01_RK)
          res = res * x + (-2.1260486589396638E-01_RK)
          res = res * x + (1.9761921280113542E-01_RK)
          res = res * x + (8.6385435215153961E-01_RK)
        end if
      end if
    elseif (x <= 1.0000000000000000E+00_RK) then
      if (x < 6.6666666666666663E-01_RK) then
        if (x < 4.3096440627115085E-01_RK) then
          res = (9.2989836420019442E-02_RK) * x + (-1.3124380721245635E-01_RK)
          res = res * x + (3.9551533881730334E-02_RK)
          res = res * x + (1.4461452858482687E-01_RK)
          res = res * x + (8.6810669969142074E-01_RK)
        else
          res = (-6.7055964504391377E-03_RK) * x + (2.2911999129408750E-02_RK)
          res = res * x + (-5.0039084385563551E-02_RK)
          res = res * x + (1.6784857275128981E-01_RK)
          res = res * x + (8.6583329929197761E-01_RK)
        end if
      else
        if (x < 9.0236892706218252E-01_RK) then
          res = (-1.0881827591406440E-03_RK) * x + (9.1295179354582787E-03_RK)
          res = res * x + (-3.7133080284270058E-02_RK)
          res = res * x + (1.6236757221824985E-01_RK)
          res = res * x + (8.6672538337508409E-01_RK)
        else
          res = (4.6017500027684044E-03_RK) * x + (-1.2962733889034270E-02_RK)
          res = res * x + (-5.1111183599694618E-03_RK)
          res = res * x + (1.4181596019335754E-01_RK)
          res = res * x + (8.7165614205287767E-01_RK)
        end if
      end if
    else
      res = cosh_acosh(ONE / x)
      return
    end if
    do i = 1, MAXITER
      yy = res * res
      df = 12.0_RK * yy - 3.0_RK
      print *, i, res, df
      if (ABS(df) < 1E-18_RK) return
      df = ((4.0_RK * yy - 3.0_RK) * res - x) / df
      res = res - df
      if (ABS(df) < THRESHOLD) return
    end do
  end function cos_acos
!
  function cosh_acosh(x) result(res)
    real(RK), intent(in) :: x
    real(RK)             :: res, yy, df
    if (x < ZERO) then
      res = ONE
      return
    elseif (x > ONE) then
      res = cos_acos(ONE / x)
      return
    else
      if (x < 5.0000000000000000E-01_RK) then
        if (x < 1.4644660940672627E-01_RK) then
          res = (2.9546700898745084E+00_RK) * x + (-6.8732996465706475E-01_RK)
          res = res * x + (-1.3519484012953734E-02_RK)
          res = res * x + (-4.6741735498024902E-03_RK)
          res = res * x + (0.0000000000000000E+00_RK)
        else
          res = (-1.0623675833993933E-01_RK) * x + (1.5654255346132895E-01_RK)
          res = res * x + (-9.6794886272540959E-02_RK)
          res = res * x + (1.8177193622178180E-03_RK)
          res = res * x + (-4.0727563438808847E-04_RK)
        end if
      else
        if (x < 8.5355339059327373E-01_RK) then
          res = (2.3919197448076808E-02_RK) * x + (-5.9903084576736369E-02_RK)
          res = res * x + (5.0794386881566130E-02_RK)
          res = res * x + (-4.7880140922667278E-02_RK)
          res = res * x + (6.4652937375348487E-03_RK)
        else
          res = (-1.3435433998222815E-01_RK) * x + (5.0040317286635916E-01_RK)
          res = res * x + (-6.9989274807980062E-01_RK)
          res = res * x + (4.0248132481808807E-01_RK)
          res = res * x + (-9.5448197561904868E-02_RK)
        end if
      end if
    end if
    !yy = invcbrt(0.5_RK * x)
    call invcbrt(0.5_RK * x, yy, res)
    res = res + 0.5_RK * (1.0_RK / yy + yy)
    yy = res * res
    df = 12.0_RK * yy - 3.0_RK
    if (ABS(df) < 1E-18_RK) return
    df = ((4.0_RK * yy - 3.0_RK) * res - ONE / x) / df
    res = res - df
  end function cosh_acosh
!
  pure elemental function sinh_asinh(x) result(res)
    real(RK), intent(in) :: x
    real(RK)             :: res, yy
    yy = x * x
    if (yy < 1.69_RK) then ! 1.3^2
      res = x * (1.0_RK / 3.0_RK) - x * yy * (4.0_RK / 81.0_RK)
    else
      !yy = invcbrt(2.0_RK * x)
      call invcbrt(2.0_RK * x, yy, res)
      res = 0.5_RK * (1.0_RK / yy - yy)
    end if
    yy = res * res
    res = res - ((4.0_RK * yy + 3.0_RK) * res - x) / (12.0_RK * yy + 3.0_RK)
    yy = res * res
    res = res - ((4.0_RK * yy + 3.0_RK) * res - x) / (12.0_RK * yy + 3.0_RK)
  end function sinh_asinh
!
! https://www.mdpi.com/1996-1073/14/4/1058
! pure elemental function invcbrt(x) result(res)
!   use mod_kinds, only: I4, R4, I8
!   implicit none
!   real(RK), intent(in) :: x
!   real(RK)             :: res, c
!   real(R4)             :: y
!   y = ABS(x)
!   res = SIGN(1.0_RK, x) * TRANSFER(INT(z"548C2B4B", I4) - TRANSFER(y, 0_I4) / 3, y)
!   c = res * res * res * x
!   res = res * (1.752319676_RK - c * (1.2509524245_RK - 0.5093818292_RK * c))
!   c = 1.0_RK - res * res * res * x
!   res = res * (1.0_RK + 0.333333333333333_RK * c)
! end function invcbrt
!
! https://www.mdpi.com/1996-1073/14/4/1058
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


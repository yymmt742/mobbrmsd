!| mobbrmsd_SLASV2 computes the singular value decomposition of a 2-by-2 triangular matrix.
!
!  mobbrmsd_SLASV2 computes the singular value decomposition of a 2-by-2
!  triangular matrix
!
!  \[
!     \left(
!       \begin{array}{}
!          F & G \\
!          0 & H \\
!       \end{array}
!     \right)
!  \]
!
!  On return, abs(SSMAX) is the larger singular value,
!  abs(SSMIN) is the smaller singular value,
!  and (CSL,SNL) and (CSR,SNR) are the left and
!  right singular vectors for abs(SSMAX),
!  giving the decomposition
!
!  \[
!     \left(
!       \begin{array}{}
!           \cos(L) & \sin(L) \\
!          -\sin(L) & \cos(L) \\
!       \end{array}
!     \right)
!     \left(
!       \begin{array}{}
!          F & G \\
!          0 & H \\
!       \end{array}
!     \right)
!     \left(
!       \begin{array}{}
!           \cos(L) & -\sin(L) \\
!           \sin(L) &  \cos(L) \\
!       \end{array}
!     \right)
!     =
!     \left(
!       \begin{array}{}
!           \text{SSMAX} &     0     \\
!                 0      &  t{SSMIN} \\
!       \end{array}
!     \right)
!  \]
!
!  Any input parameter may be aliased with any output parameter.
!
!  Barring over/underflow and assuming a guard digit in subtraction, all
!  output quantities are correct to within a few units in the last
!  place (ulps).
!
!  In IEEE arithmetic, the code works correctly if one matrix element is
!  infinite.
!
!  Overflow will not occur unless the largest singular value itself
!  overflows or is within a few ulps of overflow. (On machines with
!  partial overflow, like the Cray, overflow may occur if the largest
!  singular value is within a factor of 2 of overflow.)
!
!  Underflow is harmless if underflow is gradual. Otherwise, results
!  may correspond to a matrix modified by perturbations of size near
!  the underflow threshold.
!
!  Reference SLASV2 is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
pure elemental subroutine mobbrmsd_SLASV2(F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL)
  implicit none
  real(RK), intent(in)  :: F
!!  The (1,1) element of the 2-by-2 matrix.
!!
  real(RK), intent(in)  :: G
!!  The (1,2) element of the 2-by-2 matrix.
!!
  real(RK), intent(in)  :: H
!!  The (2,2) element of the 2-by-2 matrix.
!!
  real(RK), intent(out) :: SSMIN
!!  abs(SSMIN) is the smaller singular value.
!!
  real(RK), intent(out) :: SSMAX
!!  abs(SSMAX) is the larger singular value.
!!
  real(RK), intent(out) :: SNR
!!  The vector (CSR, SNR) is a unit right singular vector for the
!!  singular value abs(SSMAX).
!!
  real(RK), intent(out) :: CSR
!!  The vector (CSR, SNR) is a unit right singular vector for the
!!  singular value abs(SSMAX).
!!
  real(RK), intent(out) :: SNL
!!  The vector (CSL, SNL) is a unit left singular vector for the
!!  singular value abs(SSMAX).
!!
  real(RK), intent(out) :: CSL
!!  The vector (CSL, SNL) is a unit left singular vector for the
!!  singular value abs(SSMAX).
!!
  logical   :: GASMAL, SWAP
  integer   :: PMAX
  real(RK)  :: A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M
  real(RK)  :: MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT
  intrinsic :: ABS, SIGN, SQRT
! interface
!   include 'slamch.h'
! end interface
!
  FT = F
  FA = ABS(FT)
  HT = H
  HA = ABS(H)
!
! PMAX points to the maximum absolute element of matrix
! PMAX = 1 if F largest in absolute values
! PMAX = 2 if G largest in absolute values
! PMAX = 3 if H largest in absolute values
!
  PMAX = 1
  SWAP = (HA > FA)
  if (SWAP) then
    PMAX = 3
    TEMP = FT
    FT = HT
    HT = TEMP
    TEMP = FA
    FA = HA
    HA = TEMP
!
! Now FA >= HA
!
  end if
  GT = G
  GA = ABS(GT)
  if (GA == ZERO) then
!
! Diagonal matrix
!
    SSMIN = HA
    SSMAX = FA
    CLT = ONE
    CRT = ONE
    SLT = ZERO
    SRT = ZERO
  else
    GASMAL = .true.
    if (GA > FA) then
      PMAX = 2
      if ((FA / GA) < mobbrmsd_SLAMCH('EPS')) then
!
! case of very large GA
!
        GASMAL = .false.
        SSMAX = GA
        if (HA > ONE) then
          SSMIN = FA / (GA / HA)
        else
          SSMIN = (FA / GA) * HA
        end if
        CLT = ONE
        SLT = HT / GT
        SRT = ONE
        CRT = FT / GT
      end if
    end if
    if (GASMAL) then
!
! Normal case
!
      D = FA - HA
      if (D == FA) then
!
! Copes with infinite F or H
!
        L = ONE
      else
        L = D / FA
      end if
!
! Note that 0 <= L <= 1
!
      M = GT / FT
!
! Note that ABS(M) <= 1 / macheps
!
      T = TWO - L
!
! Note that T >= 1
!
      MM = M * M
      TT = T * T
      S = SQRT(TT + MM)
!
! Note that 1 <= S <= 1 + 1 / macheps
!
      if (L == ZERO) then
        R = ABS(M)
      else
        R = SQRT(L * L + MM)
      end if
!
! Note that 0 <= R <= 1 + 1 / macheps
!
      A = HALF * (S + R)
!
! Note that 1 <= A <= 1 + ABS(M)
!
      SSMIN = HA / A
      SSMAX = FA * A
      if (MM == ZERO) then
!
! Note that M is very tiny
!
        if (L == ZERO) then
          T = SIGN(TWO, FT) * SIGN(ONE, GT)
        else
          T = GT / SIGN(D, FT) + M / T
        end if
      else
        T = (M / (S + T) + M / (R + L)) * (ONE + A)
      end if
      L = SQRT(T * T + FOUR)
      CRT = TWO / L
      SRT = T / L
      CLT = (CRT + SRT * M) / A
      SLT = (HT / FT) * SRT / A
    end if
  end if
  if (SWAP) then
    CSL = SRT
    SNL = CRT
    CSR = SLT
    SNR = CLT
  else
    CSL = CLT
    SNL = SLT
    CSR = CRT
    SNR = SRT
  end if
!
! Correct signs of SSMAX and SSMIN
!
  if (PMAX == 1) TSIGN = SIGN(ONE, CSR) * SIGN(ONE, CSL) * SIGN(ONE, F)
  if (PMAX == 2) TSIGN = SIGN(ONE, SNR) * SIGN(ONE, CSL) * SIGN(ONE, G)
  if (PMAX == 3) TSIGN = SIGN(ONE, SNR) * SIGN(ONE, SNL) * SIGN(ONE, H)
  SSMAX = SIGN(SSMAX, TSIGN)
  SSMIN = SIGN(SSMIN, TSIGN * SIGN(ONE, F) * SIGN(ONE, H))
  return
!
! end of mobbrmsd_SLASV2
!
end


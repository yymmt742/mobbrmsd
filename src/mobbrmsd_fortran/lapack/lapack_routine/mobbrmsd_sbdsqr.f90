!| computes SVD using zero-shift QR algorithm.
!
!  mobbrmsd_SBDSQR computes the singular values and,
!  optionally, the right and/or left singular vectors
!  from the singular value decomposition (SVD) of
!  a real \( n \) -by- \( n \) (upper or lower)
!  bidiagonal matrix \( B \) using the implicit
!  zero-shift QR algorithm.
!
!  The SVD of \( B \) has the form
!
!  \[ B = Q S P ^ {\top} \]
!
!  where \( S \) is the diagonal matrix of singular values,
!  \( Q \) is an orthogonal matrix of left singular vectors,
!  and \( P \) is an orthogonal matrix of
!  right singular vectors.
!  If left singular vectors are requested, this
!  subroutine actually returns \( UQ \) instead of \( Q \),
!  and, if right singular vectors are requested,
!  this subroutine returns \( P ^ {\top} V ^ {\top} \)
!  instead of \( P ^ {\top} \),
!  for given real input matrices \( U \) and \( V ^ {\top} \).
!  When \( U \) and \( V ^ {\top} \) are the
!  orthogonal matrices that reduce a general matrix \( A \)
!  to bidiagonal form:  \( A = U B V ^ {\top} \),
!  as computed by mobbrmsd_DGEBRD, then
!
!  \[ A = (U Q) S (P ^ {\top} V) ^ {\top} \]
!
!  is the SVD of \( A \).
!  Optionally, the subroutine may also compute \( Q ^ {\top} C \)
!  for a given real input matrix \( C \).
!
!  See "Computing  Small Singular Values of Bidiagonal Matrices With
!  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
!  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
!  no. 5, pp. 873-912, Sept 1990) and
!  "Accurate singular values and differential qd algorithms," by
!  B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
!  Department, University of California at Berkeley, July 1992
!  for a detailed description of the algorithm.
!
!  reference SBDSQR is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!   Bug report from Cezary Dendek.
!   On March 23rd 2017, the INTEGER variable MAXIT = MAXITR*N ** 2 is
!   removed since it can overflow pretty easily (for N larger or equal
!   than 18,919). We instead use MAXITDIVN = MAXITR * N.
!
pure subroutine mobbrmsd_SBDSQR(UPLO, N, NCVT, NRU, NCC, D, E, VT, &
     &                          LDVT, U, LDU, C, LDC, WORK, INFO)
  implicit none
  character, intent(in) :: UPLO
!!  = 'U':  B is upper bidiagonal.
!!
!!  = 'L':  B is lower bidiagonal.
!!
  integer, intent(in)   :: N
!!  The order of the matrix B.  N >= 0.
!!
  integer, intent(in)   :: NCVT
!!  The number of columns of the matrix VT. NCVT >= 0.
!!
  integer, intent(in)   :: NRU
!!  The number of rows of the matrix U. NRU >= 0.
!!
  integer, intent(in)   :: NCC
!!  The number of columns of the matrix C. NCC >= 0.
!!
  real(RK), intent(inout) :: D(*)
!!  REAL array, dimension (N)
!!
!!  On entry, the n diagonal elements of the bidiagonal matrix B.
!!
!!  On exit, if INFO=0, the singular values of B in decreasing
!!  order.
!!
  real(RK), intent(inout) :: E(*)
!!  REAL array, dimension (N-1)
!!
!!  On entry, the N-1 offdiagonal elements of the bidiagonal
!!  matrix B.
!!
!!  On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E
!!  will contain the diagonal and superdiagonal elements of a
!!  bidiagonal matrix orthogonally equivalent to the one given
!!  as input.
!!
  integer, intent(in)   :: LDVT
!!  The leading dimension of the array VT.
!!  LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
!!
  real(RK), intent(inout) :: VT(LDVT, *)
!!  On entry, an N-by-NCVT matrix VT.
!!
!!  On exit, VT is overwritten by P**T * VT.
!!  Not referenced if NCVT = 0.
!!
  integer, intent(in)   :: LDU
!!  The leading dimension of the array U.  LDU >= max(1,NRU).
!!
  real(RK), intent(inout) :: U(LDU, *)
!!  On entry, an NRU-by-N matrix U.
!!
!!  On exit, U is overwritten by U * Q.
!!  Not referenced if NRU = 0.
!!
  integer, intent(in)   :: LDC
!!  The leading dimension of the array C. <br>
!!  LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
!!
  real(RK), intent(inout) :: C(LDC, *)
!!  On entry, an N-by-NCC matrix C.
!!
!!  On exit, C is overwritten by Q**T * C.
!!  Not referenced if NCC = 0.
!!
  real(RK), intent(out)   :: WORK(*)
!!  REAL array, dimension (4*(N-1))
!!
  integer, intent(out)  :: INFO
!!  = 0:  successful exit.
!!
!!  < 0:  If INFO = -i, the i-th argument had an illegal value.
!!
!!  \> 0: if NCVT = NRU = NCC = 0,
!!
!!  = 1, a split was marked by a positive value in E <br>
!!
!!  = 2, current block of Z not diagonalized after 30*N
!!       iterations (in inner while loop) <br>
!!
!!  = 3, termination criterion of outer while loop not met
!!       (program created more than N unreduced blocks) <br>
!!
!!   else NCVT = NRU = NCC = 0,
!!
!!   the algorithm did not converge; D and E contain the
!!   elements of a bidiagonal matrix which is orthogonally
!!   similar to the input matrix B;  if INFO = i, i
!!   elements of E have not converged to zero.
!!
  real(RK), parameter :: NEGONE = -1.0E0
  real(RK), parameter :: HNDRTH = 0.01E0
  real(RK), parameter :: HNDRD = 100.0E0
  real(RK), parameter :: MEIGTH = -0.125E0
  integer, parameter  :: MAXITR = 6
!..
!.. Local Scalars ..
  logical  :: LOWER, ROTATE
  integer  :: I, IDIR, ISUB, ITER, ITERDIVN, J, LL, LLL, M
  integer  :: MAXITDIVN, NM1, NM12, NM13, OLDLL, OLDM
  real(RK) :: ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU
  real(RK) :: OLDCS, OLDSN, R, SHIFT, SIGMN, SIGMX, SINL
  real(RK) :: SINR, SLL, SMAX, SMIN, SMINOA
  real(RK) :: SN, THRESH, TOL, TOLMUL, UNFL
!..
! interface
! .. External Functions ..
!   include 'lsame.h'
!   include 'slamch.h'
! .. External Subroutines ..
!   include 'slartg.h'
!   include 'slas2.h'
!   include 'slasq1.h'
!   include 'slasr.h'
!   include 'slasv2.h'
!   include 'srot.h'
!   include 'sscal.h'
!   include 'sswap.h'
! end interface
!..
!.. Intrinsic Functions ..
  intrinsic        :: ABS, MAX, MIN, real, SIGN, SQRT
!..
!.. Executable Statements ..
!
! Test the input parameters.
!
  INFO = 0
  LOWER = mobbrmsd_LSAME(UPLO, 'L')
  if (.not. mobbrmsd_LSAME(UPLO, 'U') .and. .not. LOWER) then
    INFO = -1
  else if (N < 0) then
    INFO = -2
  else if (NCVT < 0) then
    INFO = -3
  else if (NRU < 0) then
    INFO = -4
  else if (NCC < 0) then
    INFO = -5
  else if ((NCVT == 0 .and. LDVT < 1) .or. (NCVT > 0 .and. LDVT < MAX(1, N))) then
    INFO = -9
  else if (LDU < MAX(1, NRU)) then
    INFO = -11
  else if ((NCC == 0 .and. LDC < 1) .or. (NCC > 0 .and. LDC < MAX(1, N))) then
    INFO = -13
  end if
  if (INFO /= 0) then
    return
  end if
  if (N == 0) return
  if (N == 1) GO TO 160
!
!     ROTATE is true if any singular vectors desired, false otherwise
!
  ROTATE = (NCVT > 0) .or. (NRU > 0) .or. (NCC > 0)
!
!     If no singular vectors desired, use qd algorithm
!
  if (.not. ROTATE) then
    call mobbrmsd_SLASQ1(N, D, E, WORK, INFO)
!
!     If INFO equals 2, dqds didn't finish, try to finish
!
    if (INFO /= 2) return
    INFO = 0
  end if
!
  NM1 = N - 1
  NM12 = NM1 + NM1
  NM13 = NM12 + NM1
  IDIR = 0
!
!     Get machine constants
!
  EPS = mobbrmsd_SLAMCH('Epsilon')
  UNFL = mobbrmsd_SLAMCH('Safe minimum')
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left
!
  if (LOWER) then
    do I = 1, N - 1
      call mobbrmsd_SLARTG(D(I), E(I), CS, SN, R)
      D(I) = R
      E(I) = SN * D(I + 1)
      D(I + 1) = CS * D(I + 1)
      WORK(I) = CS
      WORK(NM1 + I) = SN
    end do
!
!        Update singular vectors if desired
!
    if (NRU > 0) &
&      call mobbrmsd_SLASR('R', 'V', 'F', NRU, N, WORK(1), WORK(N), &
&                  U, LDU)
    if (NCC > 0) &
&      call mobbrmsd_SLASR('L', 'V', 'F', N, NCC, WORK(1), WORK(N), &
&                  C, LDC)
  end if
!
!     Compute singular values to relative accuracy TOL
!     (By setting TOL to be negative, algorithm will compute
!     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
!
  TOLMUL = MAX(TEN, MIN(HUNDRD, EPS**MEIGTH))
  TOL = TOLMUL * EPS
!
!     Compute approximate maximum, minimum singular values
!
  SMAX = ZERO
  do I = 1, N
    SMAX = MAX(SMAX, ABS(D(I)))
  end do
  do I = 1, N - 1
    SMAX = MAX(SMAX, ABS(E(I)))
  end do
  SMIN = ZERO
  if (TOL >= ZERO) then
!
!        Relative accuracy desired
!
    SMINOA = ABS(D(1))
    if (SMINOA /= ZERO) then
      MU = SMINOA
      do I = 2, N
        MU = ABS(D(I)) * (MU / (MU + ABS(E(I - 1))))
        SMINOA = MIN(SMINOA, MU)
        if (SMINOA == ZERO) exit
      end do
    end if
    SMINOA = SMINOA / SQRT(real(N, RK))
    THRESH = MAX(TOL * SMINOA, MAXITR * (N * (N * UNFL)))
  else
!
!        Absolute accuracy desired
!
    THRESH = MAX(ABS(TOL) * SMAX, MAXITR * (N * (N * UNFL)))
!
  end if
!
!     Prepare for main iteration loop for the singular values
!     (MAXIT is the maximum number of passes through the inner
!     loop permitted before nonconvergence signalled.)
!
  MAXITDIVN = MAXITR * N
  ITERDIVN = 0
  ITER = -1
  OLDLL = -1
  OLDM = -1
!
!     M points to last element of unconverged part of matrix
!
  M = N
!
!     Begin main iteration loop
!
60 continue
!
!     Check for convergence or exceeding iteration count
!
  if (M <= 1) GO TO 160
!
  if (ITER >= N) then
    ITER = ITER - N
    ITERDIVN = ITERDIVN + 1
    if (ITERDIVN >= MAXITDIVN) GO TO 200
  end if
!
!     Find diagonal block of matrix to work on
!
  if (TOL < ZERO .and. ABS(D(M)) <= THRESH) D(M) = ZERO
  SMAX = ABS(D(M))
  do LLL = 1, M - 1
    LL = M - LLL
    ABSS = ABS(D(LL))
    ABSE = ABS(E(LL))
    if (TOL < ZERO .and. ABSS <= THRESH) D(LL) = ZERO
    if (ABSE <= THRESH) GO TO 80
    SMAX = MAX(SMAX, ABSS, ABSE)
  end do
  LL = 0
  GO TO 90
80 continue
  E(LL) = ZERO
!
! Matrix splits since E(LL) = 0
!
  if (LL == M - 1) then
!
!   Convergence of bottom singular value, return to top of loop
!
    M = M - 1
    GO TO 60
  end if
90 continue
  LL = LL + 1
!
! E(LL) through E(M-1) are nonzero, E(LL-1) is zero
!
  if (LL == M - 1) then
!
!   2 by 2 block, handle separately
!
    call mobbrmsd_SLASV2(D(M - 1), E(M - 1), D(M), SIGMN, SIGMX, SINR, &
&                        COSR, SINL, COSL)
    D(M - 1) = SIGMX
    E(M - 1) = ZERO
    D(M) = SIGMN
!
!        Compute singular vectors, if desired
!
    if (NCVT > 0) &
&      call mobbrmsd_SROT(NCVT, VT(M - 1, 1), LDVT, VT(M, 1), LDVT, COSR, SINR)
    if (NRU > 0) &
&      call mobbrmsd_SROT(NRU, U(1, M - 1), 1, U(1, M), 1, COSL, SINL)
    if (NCC > 0) &
&      call mobbrmsd_SROT(NCC, C(M - 1, 1), LDC, C(M, 1), LDC, COSL, SINL)
    M = M - 2
    GO TO 60
  end if
!
!     If working on new submatrix, choose shift direction
!     (from larger end diagonal element towards smaller)
!
  if (LL > OLDM .or. M < OLDLL) then
    if (ABS(D(LL)) >= ABS(D(M))) then
!
!           Chase bulge from top (big end) to bottom (small end)
!
      IDIR = 1
    else
!
!           Chase bulge from bottom (big end) to top (small end)
!
      IDIR = 2
    end if
  end if
!
!     Apply convergence tests
!
  if (IDIR == 1) then
!
!        Run convergence test in forward direction
!        First apply standard test to bottom of matrix
!
    if (ABS(E(M - 1)) <= ABS(TOL) * ABS(D(M)) .or. &
&       (TOL < ZERO .and. ABS(E(M - 1)) <= THRESH)) then
      E(M - 1) = ZERO
      GO TO 60
    end if
!
    if (TOL >= ZERO) then
!
!           If relative accuracy desired,
!           apply convergence criterion forward
!
      MU = ABS(D(LL))
      SMIN = MU
      do LLL = LL, M - 1
        if (ABS(E(LLL)) <= TOL * MU) then
          E(LLL) = ZERO
          GO TO 60
        end if
        MU = ABS(D(LLL + 1)) * (MU / (MU + ABS(E(LLL))))
        SMIN = MIN(SMIN, MU)
      end do
    end if
!
  else
!
!        Run convergence test in backward direction
!        First apply standard test to top of matrix
!
    if (ABS(E(LL)) <= ABS(TOL) * ABS(D(LL)) .or. &
&       (TOL < ZERO .and. ABS(E(LL)) <= THRESH)) then
      E(LL) = ZERO
      GO TO 60
    end if
!
    if (TOL >= ZERO) then
!
!           If relative accuracy desired,
!           apply convergence criterion backward
!
      MU = ABS(D(M))
      SMIN = MU
      do LLL = M - 1, LL, -1
        if (ABS(E(LLL)) <= TOL * MU) then
          E(LLL) = ZERO
          GO TO 60
        end if
        MU = ABS(D(LLL)) * (MU / (MU + ABS(E(LLL))))
        SMIN = MIN(SMIN, MU)
      end do
    end if
  end if
  OLDLL = LL
  OLDM = M
!
!     Compute shift.  First, test if shifting would ruin relative
!     accuracy, and if so set the shift to zero.
!
  if (TOL >= ZERO .and. N * TOL * (SMIN / SMAX) <= &
 &    MAX(EPS, HNDRTH * TOL)) then
!
!        Use a zero shift to avoid loss of relative accuracy
!
    SHIFT = ZERO
!
  else
!
!        Compute the shift from 2-by-2 block at end of matrix
!
    if (IDIR == 1) then
      SLL = ABS(D(LL))
      call mobbrmsd_SLAS2(D(M - 1), E(M - 1), D(M), SHIFT, R)
    else
      SLL = ABS(D(M))
      call mobbrmsd_SLAS2(D(LL), E(LL), D(LL + 1), SHIFT, R)
    end if
!
!        Test if shift negligible, and if so set to zero
!
    if (SLL > ZERO) then
      if ((SHIFT / SLL)**2 < EPS) SHIFT = ZERO
    end if
  end if
!
!     Increment iteration count
!
  ITER = ITER + M - LL
!
!     If SHIFT = 0, do simplified QR iteration
!
  if (SHIFT == ZERO) then
    if (IDIR == 1) then
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
      CS = ONE
      OLDCS = ONE
      do I = LL, M - 1
        call mobbrmsd_SLARTG(D(I) * CS, E(I), CS, SN, R)
        if (I > LL) E(I - 1) = OLDSN * R
        call mobbrmsd_SLARTG(OLDCS * R, D(I + 1) * SN, OLDCS, OLDSN, D(I))
        WORK(I - LL + 1) = CS
        WORK(I - LL + 1 + NM1) = SN
        WORK(I - LL + 1 + NM12) = OLDCS
        WORK(I - LL + 1 + NM13) = OLDSN
      end do
      H = D(M) * CS
      D(M) = H * OLDCS
      E(M - 1) = H * OLDSN
!
!           Update singular vectors
!
      if (NCVT > 0) &
&         call mobbrmsd_SLASR('L', 'V', 'F', M - LL + 1, NCVT, WORK(1), &
&                             WORK(N), VT(LL, 1), LDVT)
      if (NRU > 0) &
&         call mobbrmsd_SLASR('R', 'V', 'F', NRU, M - LL + 1, &
&                             WORK(NM12 + 1), &
&                             WORK(NM13 + 1), U(1, LL), LDU)
      if (NCC > 0) &
&         call mobbrmsd_SLASR('L', 'V', 'F', M - LL + 1, NCC, &
&                             WORK(NM12 + 1), &
&                             WORK(NM13 + 1), C(LL, 1), LDC)
!
!           Test convergence
!
      if (ABS(E(M - 1)) <= THRESH) E(M - 1) = ZERO
!
    else
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
      CS = ONE
      OLDCS = ONE
      do I = M, LL + 1, -1
        call mobbrmsd_SLARTG(D(I) * CS, E(I - 1), CS, SN, R)
        if (I < M) E(I) = OLDSN * R
        call mobbrmsd_SLARTG(OLDCS * R, D(I - 1) * SN, OLDCS, OLDSN, D(I))
        WORK(I - LL) = CS
        WORK(I - LL + NM1) = -SN
        WORK(I - LL + NM12) = OLDCS
        WORK(I - LL + NM13) = -OLDSN
      end do
      H = D(LL) * CS
      D(LL) = H * OLDCS
      E(LL) = H * OLDSN
!
!           Update singular vectors
!
      if (NCVT > 0) &
&         call mobbrmsd_SLASR('L', 'V', 'B', M - LL + 1, NCVT, &
&                             WORK(NM12 + 1), &
&                             WORK(NM13 + 1), VT(LL, 1), LDVT)
      if (NRU > 0) &
&         call mobbrmsd_SLASR('R', 'V', 'B', NRU, M - LL + 1, WORK(1), &
&                             WORK(N), U(1, LL), LDU)
      if (NCC > 0) &
&         call mobbrmsd_SLASR('L', 'V', 'B', M - LL + 1, NCC, WORK(1), &
&                             WORK(N), C(LL, 1), LDC)
!
!           Test convergence
!
      if (ABS(E(LL)) <= THRESH) E(LL) = ZERO
    end if
  else
!
!        Use nonzero shift
!
    if (IDIR == 1) then
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
      F = (ABS(D(LL)) - SHIFT) * (SIGN(ONE, D(LL)) + SHIFT / D(LL))
      G = E(LL)
      do I = LL, M - 1
        call mobbrmsd_SLARTG(F, G, COSR, SINR, R)
        if (I > LL) E(I - 1) = R
        F = COSR * D(I) + SINR * E(I)
        E(I) = COSR * E(I) - SINR * D(I)
        G = SINR * D(I + 1)
        D(I + 1) = COSR * D(I + 1)
        call mobbrmsd_SLARTG(F, G, COSL, SINL, R)
        D(I) = R
        F = COSL * E(I) + SINL * D(I + 1)
        D(I + 1) = COSL * D(I + 1) - SINL * E(I)
        if (I < M - 1) then
          G = SINL * E(I + 1)
          E(I + 1) = COSL * E(I + 1)
        end if
        WORK(I - LL + 1) = COSR
        WORK(I - LL + 1 + NM1) = SINR
        WORK(I - LL + 1 + NM12) = COSL
        WORK(I - LL + 1 + NM13) = SINL
      end do
      E(M - 1) = F
!
!           Update singular vectors
!
      if (NCVT > 0) &
&         call mobbrmsd_SLASR('L', 'V', 'F', M - LL + 1, NCVT, WORK(1), &
&                     WORK(N), VT(LL, 1), LDVT)
      if (NRU > 0) &
&         call mobbrmsd_SLASR('R', 'V', 'F', NRU, M - LL + 1, &
&                     WORK(NM12 + 1), &
&                     WORK(NM13 + 1), U(1, LL), LDU)
      if (NCC > 0) &
&         call mobbrmsd_SLASR('L', 'V', 'F', M - LL + 1, NCC, &
&                     WORK(NM12 + 1), &
&                     WORK(NM13 + 1), C(LL, 1), LDC)
!
!           Test convergence
!
      if (ABS(E(M - 1)) <= THRESH) E(M - 1) = ZERO
!
    else
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
      F = (ABS(D(M)) - SHIFT) * (SIGN(ONE, D(M)) + SHIFT / D(M))
      G = E(M - 1)
      do I = M, LL + 1, -1
        call mobbrmsd_SLARTG(F, G, COSR, SINR, R)
        if (I < M) E(I) = R
        F = COSR * D(I) + SINR * E(I - 1)
        E(I - 1) = COSR * E(I - 1) - SINR * D(I)
        G = SINR * D(I - 1)
        D(I - 1) = COSR * D(I - 1)
        call mobbrmsd_SLARTG(F, G, COSL, SINL, R)
        D(I) = R
        F = COSL * E(I - 1) + SINL * D(I - 1)
        D(I - 1) = COSL * D(I - 1) - SINL * E(I - 1)
        if (I > LL + 1) then
          G = SINL * E(I - 2)
          E(I - 2) = COSL * E(I - 2)
        end if
        WORK(I - LL) = COSR
        WORK(I - LL + NM1) = -SINR
        WORK(I - LL + NM12) = COSL
        WORK(I - LL + NM13) = -SINL
      end do
      E(LL) = F
!
!           Test convergence
!
      if (ABS(E(LL)) <= THRESH) E(LL) = ZERO
!
!           Update singular vectors if desired
!
      if (NCVT > 0) &
&         call mobbrmsd_SLASR('L', 'V', 'B', M - LL + 1, NCVT, &
&                     WORK(NM12 + 1), &
&                     WORK(NM13 + 1), VT(LL, 1), LDVT)
      if (NRU > 0) &
&         call mobbrmsd_SLASR('R', 'V', 'B', NRU, M - LL + 1, WORK(1), &
&                     WORK(N), U(1, LL), LDU)
      if (NCC > 0) &
&         call mobbrmsd_SLASR('L', 'V', 'B', M - LL + 1, NCC, WORK(1), &
&                     WORK(N), C(LL, 1), LDC)
    end if
  end if
!
!     QR iteration finished, go back and check convergence
!
  GO TO 60
!
!     All singular values converged, so make them positive
!
160 continue
  do I = 1, N
    if (D(I) == ZERO) then
!
!           Avoid -ZERO
!
      D(I) = ZERO
    end if
    if (D(I) < ZERO) then
      D(I) = -D(I)
!
!           Change sign of singular vectors, if desired
!
      if (NCVT > 0) call mobbrmsd_SSCAL(NCVT, NEGONE, VT(I, 1), LDVT)
    end if
  end do
!
!     Sort the singular values into decreasing order (insertion sort on
!     singular values, but only one transposition per singular vector)
!
  do I = 1, N - 1
!
!        Scan for smallest D(I)
!
    ISUB = 1
    SMIN = D(1)
    do J = 2, N + 1 - I
      if (D(J) <= SMIN) then
        ISUB = J
        SMIN = D(J)
      end if
    end do
    if (ISUB /= N + 1 - I) then
!
!           Swap singular values and vectors
!
      D(ISUB) = D(N + 1 - I)
      D(N + 1 - I) = SMIN
      if (NCVT > 0) &
&         call mobbrmsd_SSWAP(NCVT, VT(ISUB, 1), LDVT, VT(N + 1 - I, 1), LDVT)
      if (NRU > 0) &
&         call mobbrmsd_SSWAP(NRU, U(1, ISUB), 1, U(1, N + 1 - I), 1)
      if (NCC > 0) &
&         call mobbrmsd_SSWAP(NCC, C(ISUB, 1), LDC, C(N + 1 - I, 1), LDC)
    end if
  end do
  GO TO 220
!
!     Maximum number of iterations exceeded, failure to converge
!
200 continue
  INFO = 0
  do I = 1, N - 1
    if (E(I) /= ZERO) INFO = INFO + 1
  end do
220 continue
  return
!
!     End of SBDSQR
!
end subroutine mobbrmsd_SBDSQR


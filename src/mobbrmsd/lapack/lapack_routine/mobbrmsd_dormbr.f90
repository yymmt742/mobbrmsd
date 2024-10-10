!| multiply an orthogonal matrix to C.
!
!  If VECT = 'Q', mobbrmsd_DORMBR overwrites
!  the general real M-by-N matrix C with
!
!  |             | SIDE = 'L'       | SIDE = 'R'       |
!  | :---:       | :---:            | :---:            |
!  | TRANS = 'N' | \( Q  C \)       | \( C Q \)        |
!  | TRANS = 'T' | \( Q ^ \top C \) | \( C Q ^ \top \) |
!
!  If VECT = 'P', mobbrmsd_DORMBR overwrites
!  the general real M-by-N matrix C with
!
!  |             | SIDE = 'L'       | SIDE = 'R'       |
!  | :---:       | :---:            | :---:            |
!  | TRANS = 'N' | \( P  C \)       | \( C P \)        |
!  | TRANS = 'T' | \( P ^ \top C \) | \( C P ^ \top \) |
!
!  Here \( Q \) and \( P ^ \top \) are the orthogonal matrices determined
!  by mobbrmsd_DGEBRD when reducing
!  a real matrix \( A \) to bidiagonal form:
!  \( A = Q B P ^ \top \).
!  \( Q \) and \( P ^ \top \) are defined as products of
!  elementary reflectors \( H _ i \) and \( G _ i \) respectively.
!
!  Let \( n _ q = m \) if SIDE = 'L'
!  and \( n _ q = n \) if SIDE = 'R'.
!  Thus \( n _ q \)is the order of the orthogonal matrix
!  \( Q \) or \( P ^ \top \) that is applied.
!
!  If VECT = 'Q', \( A \) is assumed to have been
!  an \( n _ q \)-by-\( k \) matrix:
!  if \( n _ q \ge k \), \( Q = H _ 1 H _ 2 \cdots H _ k \);
!  if \( n _ q < k \), \( Q = H _ 1 H _ 2 \cdots H _ {n _ q - 1} \).
!
!  If VECT = 'P', \( A \) is assumed to have been
!  an \( k \)-by-\( n _ q \) matrix:
!  if \( n _ q < k \), \( P = G _ 1 G _ 2 \cdots G _ k \);
!  if \( n _ q \ge k \), \( P = G _ 1 G _ 2 \cdots G _ {n _ q - 1} \).
!
!  Reference DORMBR is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DORMBR(VECT, SIDE, TRANS, M, N, K, A, LDA, &
               &                TAU, C, LDC, WORK, LWORK, INFO)
  implicit none
  character, intent(in) :: VECT
!!  = 'Q': apply Q or Q**T;
!!  = 'P': apply P or P**T.
!!
  character, intent(in) :: SIDE
!!  = 'L': apply Q, Q**T, P or P**T from the Left;
!!  = 'R': apply Q, Q**T, P or P**T from the Right.
!!
  character, intent(in) :: TRANS
!!  = 'N':  No transpose, apply Q  or P;
!!
!!  = 'T':  Transpose, apply Q**T or P**T.
!!
  integer, intent(in)     :: M
!!  The number of rows of the matrix C. M >= 0.
!!
  integer, intent(in)     :: N
!!  The number of columns of the matrix C. N >= 0.
!!
  integer, intent(in)     :: K
!!  If VECT = 'Q', the number of columns in the original
!!  matrix reduced by mobbrmsd_DGEBRD.
!!
!!  If VECT = 'P', the number of rows in the original
!!  matrix reduced by mobbrmsd_DGEBRD.
!!  K >= 0.
!!
  integer, intent(in)     :: LDA
!!  The leading dimension of the array A.
!!
!!  If VECT = 'Q', LDA >= max(1,nq);
!!
!!  if VECT = 'P', LDA >= max(1,min(nq,K)).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!  DOUBLE PRECISION array, dimension
!!
!!  (LDA,min(nq,K)) if VECT = 'Q'
!!
!!  (LDA,nq)        if VECT = 'P'
!!
!!  The vectors which define the elementary reflectors
!!  H(i) and G(i), whose products determine
!!  the matrices Q and P,
!!  as returned by mobbrmsd_DGEBRD.
!!
  real(RK), intent(in)    :: TAU(*)
!!  DOUBLE PRECISION array, dimension (min(nq,K))
!!
!!  TAU(i) must contain the scalar factor of the elementary
!!  reflector H(i) or G(i) which determines Q or P, as returned
!!  by mobbrmsd_DGEBRD in the array argument TAUQ or TAUP.
!!
  integer, intent(in)     :: LDC
!!  The leading dimension of the array C. LDC >= max(1,M).
!!
  real(RK), intent(inout) :: C(LDC, *)
!!  DOUBLE PRECISION array, dimension (LDC,N)
!!
!!  On entry, the M-by-N matrix C.
!!
!!  On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
!!  or P*C or P**T*C or C*P or C*P**T.
!!
  integer, intent(in)     :: LWORK
!!  The dimension of the array WORK.
!!
!!  If SIDE = 'L', LWORK >= max(1,N);
!!
!!  if SIDE = 'R', LWORK >= max(1,M).
!!
!!  For optimum performance LWORK >= N*NB if SIDE = 'L', and
!!  LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!!  blocksize.
!!
!!  If LWORK = -1, then a workspace query is assumed; the routine
!!  only calculates the optimal size of the WORK array, returns
!!  this value as the first entry of the WORK array.
!!
  real(RK), intent(out)   :: WORK(*)
!!  DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!
!!  On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!
  integer, intent(out)    :: INFO
!!  = 0:  successful exit
!!
!!  < 0:  if INFO = -i, the i-th argument had an illegal value
!!
  logical                 :: APPLYQ, LEFT, LQUERY, NOTRAN
  character(1)            :: TRANST
  integer                 :: I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
  intrinsic               :: MAX, MIN
! interface
!   include 'dormlq.h'
!   include 'dormqr.h'
!   include 'lsame.h'
!   include 'ilaenv.h'
! end interface
!
! Test the input arguments
!
  INFO = 0
  APPLYQ = mobbrmsd_LSAME(VECT, 'Q')
  LEFT = mobbrmsd_LSAME(SIDE, 'L')
  NOTRAN = mobbrmsd_LSAME(TRANS, 'N')
  LQUERY = (LWORK == -1)
!
! NQ is the order of Q or P and NW is the minimum dimension of WORK
!
  if (LEFT) then
    NQ = M
    NW = MAX(1, N)
  else
    NQ = N
    NW = MAX(1, M)
  end if
  if (.not. APPLYQ .and. .not. mobbrmsd_LSAME(VECT, 'P')) then
    INFO = -1
  else if (.not. LEFT .and. .not. mobbrmsd_LSAME(SIDE, 'R')) then
    INFO = -2
  else if (.not. NOTRAN .and. .not. mobbrmsd_LSAME(TRANS, 'T')) then
    INFO = -3
  else if (M < 0) then
    INFO = -4
  else if (N < 0) then
    INFO = -5
  else if (K < 0) then
    INFO = -6
  else if ((APPLYQ .and. LDA < MAX(1, NQ)) .or. &
 &         (.not. APPLYQ .and. LDA < MAX(1, MIN(NQ, K)))) &
 &          then
    INFO = -8
  else if (LDC < MAX(1, M)) then
    INFO = -11
  else if (LWORK < NW .and. .not. LQUERY) then
    INFO = -13
  end if
!
  if (INFO == 0) then
    if (APPLYQ) then
      if (LEFT) then
        NB = mobbrmsd_ILAENV(1, 'DORMQR', SIDE//TRANS, M - 1, N, M - 1, -1)
      else
        NB = mobbrmsd_ILAENV(1, 'DORMQR', SIDE//TRANS, M, N - 1, N - 1, -1)
      end if
    else
      if (LEFT) then
        NB = mobbrmsd_ILAENV(1, 'DORMLQ', SIDE//TRANS, M - 1, N, M - 1, -1)
      else
        NB = mobbrmsd_ILAENV(1, 'DORMLQ', SIDE//TRANS, M, N - 1, N - 1, -1)
      end if
    end if
    LWKOPT = NW * NB
    WORK(1) = LWKOPT
  end if
!
  if (INFO /= 0) then
    return
  else if (LQUERY) then
    return
  end if
!
! Quick return if possible
!
  WORK(1) = 1
  if (M == 0 .or. N == 0) return
!
  if (APPLYQ) then
!
!   Apply Q
!
    if (NQ >= K) then
!
!     Q was determined by a call to mobbrmsd_DGEBRD with nq >= k
!
      call mobbrmsd_DORMQR(SIDE, TRANS, M, N, K, A, LDA, TAU, C, &
     &            LDC, WORK, LWORK, IINFO)
    else if (NQ > 1) then
!
!     Q was determined by a call to mobbrmsd_DGEBRD with nq < k
!
      if (LEFT) then
        MI = M - 1
        NI = N
        I1 = 2
        I2 = 1
      else
        MI = M
        NI = N - 1
        I1 = 1
        I2 = 2
      end if
      call mobbrmsd_DORMQR(SIDE, TRANS, MI, NI, NQ - 1, A(2, 1), LDA, &
     &            TAU, C(I1, I2), LDC, WORK, LWORK, IINFO)
    end if
  else
!
!   Apply P
!
    if (NOTRAN) then
      TRANST = 'T'
    else
      TRANST = 'N'
    end if
    if (NQ > K) then
!
!     P was determined by a call to mobbrmsd_DGEBRD with nq > k
!
      call mobbrmsd_DORMLQ(SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC, &
     &            WORK, LWORK, IINFO)
    else if (NQ > 1) then
!
!     P was determined by a call to mobbrmsd_DGEBRD with nq <= k
!
      if (LEFT) then
        MI = M - 1
        NI = N
        I1 = 2
        I2 = 1
      else
        MI = M
        NI = N - 1
        I1 = 1
        I2 = 2
      end if
      call mobbrmsd_DORMLQ(SIDE, TRANST, MI, NI, NQ - 1, A(1, 2), LDA, &
     &            TAU, C(I1, I2), LDC, WORK, LWORK, IINFO)
    end if
  end if
  WORK(1) = LWKOPT
  return
!
! End of mobbrmsd_DORMBR
!
end subroutine mobbrmsd_DORMBR


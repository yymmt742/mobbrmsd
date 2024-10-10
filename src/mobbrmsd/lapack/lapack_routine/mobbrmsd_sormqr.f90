!| multiply an orthogonal matrix \( Q = H _ 1 H _ 2 \cdots H _ k \).
!
!  mobbrmsd_SGELQF overwrites the general real
!  \( m \)-by-\( n \) matrix \( C \) with
!
!  | SIDE  | TRANS |                  |
!  | :---: | :---: |   :---:          |
!  |  'L'  |  'N'  | \( Q C \)        |
!  |  'R'  |  'N'  | \( C Q \)        |
!  |  'L'  |  'T'  | \( Q ^ \top C \) |
!  |  'R'  |  'T'  | \( C Q ^ \top \) |
!
!  where \( Q \) is a real orthogonal matrix defined
!  as the product of \( k \) elementary reflectors
!
!  \[
!     Q = H _ 1 H _ 2 \cdots H _ k
!  \]
!
!  as returned by mobbrmsd_SGELQF.
!  \( Q \) is of order \( m \) if SIDE = 'L'
!  and of order \( n \) if SIDE = 'R'.
!
!  Reference SORMQR is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_SORMQR(SIDE, TRANS, M, N, K, A, LDA, TAU, &
               &                C, LDC, WORK, LWORK, INFO)
  implicit none
  character, intent(in) :: SIDE
!!  = 'L': apply Q or Q**T from the Left;
!!
!!  = 'R': apply Q or Q**T from the Right.
!!
  character, intent(in) :: TRANS
!!  = 'N':  No transpose, apply Q;
!!
!!  = 'T':  Transpose, apply Q**T.
!!
  integer, intent(in)   :: M
!!  The number of rows of the matrix C. M >= 0.
!!
  integer, intent(in)   :: N
!!  The number of columns of the matrix C. N >= 0.
!!
  integer, intent(in)   :: K
!!  The number of elementary reflectors whose product defines
!!  the matrix Q.
!!
!!  If SIDE = 'L', M >= K >= 0;
!!
!!  if SIDE = 'R', N >= K >= 0.
!!
  integer, intent(in)   :: LDA
!!  The leading dimension of the array A.
!!
!!  If SIDE = 'L', LDA >= max(1,M);
!!
!!  if SIDE = 'R', LDA >= max(1,N).
!!
  real(RK), intent(inout)  :: A(LDA, *)
!!  DOUBLE PRECISION array, dimension (LDA,K)
!!
!!  The i-th column must contain the vector which defines the
!!  elementary reflector H(i), for i = 1,2,...,k, as returned by
!!  mobbrmsd_DGEQRF in the first k columns of its array argument A.
!!
  integer, intent(in)      :: LDC
!!  The leading dimension of the array C. LDC >= max(1,M).
!!
  real(RK), intent(in)     :: TAU(*)
!!  DOUBLE PRECISION array, dimension (K)
!!  TAU(i) must contain the scalar factor of the elementary
!!  reflector H(i), as returned by mobbrmsd_DGEQRF.
!!
  real(RK), intent(inout)  :: C(LDC, *)
!!  DOUBLE PRECISION array, dimension (LDC,N)
!!  On entry, the M-by-N matrix C.
!!
!!  On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!!
  real(RK), intent(out)    :: WORK(*)
!!  DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!
!!  On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!
  integer, intent(in)      :: LWORK
!!  The dimension of the array WORK.
!!
!!  If SIDE = 'L', LWORK >= max(1,N);
!!
!!  if SIDE = 'R', LWORK >= max(1,M).
!!
!!  For good performance, LWORK should generally be larger.
!!  If LWORK = -1, then a workspace query is assumed; the routine
!!  only calculates the optimal size of the WORK array, returns
!!  this value as the first entry of the WORK array, and no error
!!  message related to LWORK is issued by XERBLA.
!!
  integer, intent(out)     :: INFO
!!  = 0:  successful exit
!!
!!  < 0:  if INFO = -i, the i-th argument had an illegal value
!!
  integer, parameter :: NBMAX = 64
  integer, parameter :: LDT = NBMAX + 1
  integer, parameter :: TSIZE = LDT * NBMAX
  logical   :: LEFT, LQUERY, NOTRAN
  integer   :: I, I1, I2, I3, IB, IC, IINFO, IWT, JC
  integer   :: LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW
  intrinsic :: MAX, MIN
! interface
!   include 'lsame.h'
!   include 'ilaenv.h'
!   include 'slarfb.h'
!   include 'slarft.h'
!   include 'sorm2r.h'
! end interface
!
! Test the input arguments
!
  INFO = 0
  LEFT = mobbrmsd_LSAME(SIDE, 'L')
  NOTRAN = mobbrmsd_LSAME(TRANS, 'N')
  LQUERY = (LWORK == -1)
!
! NQ is the order of Q and NW is the minimum dimension of WORK
!
  if (LEFT) then
    NQ = M
    NW = N
  else
    NQ = N
    NW = M
  end if
  if (.not. LEFT .and. .not. mobbrmsd_LSAME(SIDE, 'R')) then
    INFO = -1
  else if (.not. NOTRAN .and. .not. mobbrmsd_LSAME(TRANS, 'T')) then
    INFO = -2
  else if (M < 0) then
    INFO = -3
  else if (N < 0) then
    INFO = -4
  else if (K < 0 .or. K > NQ) then
    INFO = -5
  else if (LDA < MAX(1, NQ)) then
    INFO = -7
  else if (LDC < MAX(1, M)) then
    INFO = -10
  else if (LWORK < MAX(1, NW) .and. .not. LQUERY) then
    INFO = -12
  end if
!
  if (INFO == 0) then
!
! Compute the workspace requirements
!
    NB = MIN(NBMAX, mobbrmsd_ILAENV(1, 'SORMQR', SIDE//TRANS, M, N, K, -1))
    LWKOPT = MAX(1, NW) * NB + TSIZE
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
  if (M == 0 .or. N == 0 .or. K == 0) then
    WORK(1) = 1
    return
  end if
!
  NBMIN = 2
  LDWORK = NW
  if (NB > 1 .and. NB < K) then
    if (LWORK < NW * NB + TSIZE) then
      NB = (LWORK - TSIZE) / LDWORK
      NBMIN = MAX(2, mobbrmsd_ILAENV(2, 'SORMQR', SIDE//TRANS, M, N, K, -1))
    end if
  end if
!
  if (NB < NBMIN .or. NB >= K) then
!
! use unblocked code
!
    call mobbrmsd_SORM2R(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO)
  else
!
! use blocked code
!
    IWT = 1 + NW * NB
    if ((LEFT .and. .not. NOTRAN) .or. (.not. LEFT .and. NOTRAN)) then
      I1 = 1
      I2 = K
      I3 = NB
    else
      I1 = ((K - 1) / NB) * NB + 1
      I2 = 1
      I3 = -NB
    end if
!
    if (LEFT) then
      NI = N
      JC = 1
    else
      MI = M
      IC = 1
    end if
!
    do I = I1, I2, I3
      IB = MIN(NB, K - I + 1)
!
! Form the triangular factor of the block reflector
! H = H(i) H(i + 1) ...H(i + ib - 1)
!
      call mobbrmsd_SLARFT('Forward', 'Columnwise', NQ - I + 1, IB, A(I, I), LDA, TAU(I), WORK(IWT), LDT)
      if (LEFT) then
!
! H or H**T is applied to C(i:m, 1:n)
!
        MI = M - I + 1
        IC = I
      else
!
! H or H**T is applied to C(1:m, i:n)
!
        NI = N - I + 1
        JC = I
      end if
!
! Apply H or H**T
!
      call mobbrmsd_SLARFB(SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
          &       IB, A(I, I), LDA, WORK(IWT), LDT, &
          &       C(IC, JC), LDC, WORK, LDWORK)
    end do
  end if
  WORK(1) = LWKOPT
  return
!
! end of mobbrmsd_SORMQR
!
end


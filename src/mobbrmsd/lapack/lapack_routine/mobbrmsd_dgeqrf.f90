!| computes a QR factorization of a real \( m \)-by-\( n \) matrix \( A \):
!
!  mobbrmsd_DGEQRF computes a QR factorization
!  of a real \( m \)-by-\( n \) matrix \( A \):
!
!  \[
!    A = Q
!    \left(
!      \begin{array}{}
!        R \\
!        0 \\
!      \end{array}
!    \right)
!  \]
!  where:
!
!  \( Q \) is a \( m \)-by-\( m \) orthogonal matrix;
!  \( R \) is an upper-triangular \( n \)-by-\( n \) matrix;
!  \( 0 \) is a \( (m-n) \)-by-\( n \) zero matrix, if \( m > n \).
!
!  The matrix \( Q \) is represented as a product of elementary reflectors
!
!  \[
!     Q = H _ 1 H _ 2 \cdots H _ k,
!  \]
!
!   where \( k = \min(m,n) \).
!   Each  \( H _ i \) has the form
!
!  \[
!     H _ i = I - \tau v v ^ \top,
!  \]
!
!   where \( \tau \) is a real scalar,
!   and \( v \) is a real vector with
!   \( v _ j = 0, 1 \le j \le i-1) \) and \( v _ i = 1 \);
!   v(i+1:m) is stored on exit in A(i+1:m,i), and \( \tau \) in TAU(i).
!
!  Reference DGEQRF is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
  implicit none
  integer, intent(in)     :: M
!!  The number of rows of the matrix A.  M >= 0.
!!
  integer, intent(in)     :: N
!!  The number of columns of the matrix A.  N >= 0.
!!
  integer, intent(in)     :: LDA
!!  The leading dimension of the array A.  LDA >= max(1,M).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!  DOUBLE PRECISION array, dimension (LDA,N)
!!
!!  On entry, the M-by-N matrix A.
!!
!!  On exit, the elements on and above the diagonal of the array
!!  contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!!  upper triangular if m >= n); the elements below the diagonal,
!!  with the array TAU, represent the orthogonal matrix Q as a
!!  product of min(m,n) elementary reflectors (see Further
!!  Details).
!!
  real(RK), intent(out)   :: TAU(*)
!!  DOUBLE PRECISION array, dimension (min(M,N))
!!  The scalar factors of the elementary reflectors (see Further
!!  Details).
!!
  real(RK), intent(out)   :: WORK(*)
!!  DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!  On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!
  integer, intent(in)     :: LWORK
!!  The dimension of the array WORK.
!!
!!  LWORK >= 1, if MIN(M,N) = 0, and LWORK >= N, otherwise.
!!  For optimum performance LWORK >= N*NB, where NB is
!!  the optimal blocksize.
!!
!!  If LWORK = -1, then a workspace query is assumed; the routine
!!  only calculates the optimal size of the WORK array, returns
!!  this value as the first entry of the WORK array.
!!
  integer, intent(out)    :: INFO
!!  = 0:  successful exit
!!
!!  < 0:  if INFO = -i, the i-th argument had an illegal value
!!
  logical :: LQUERY
  integer :: I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX
  intrinsic :: MAX, MIN
! interface
!   include 'dgeqr2.h'
!   include 'dlarfb.h'
!   include 'dlarft.h'
!   include 'ilaenv.h'
! end interface
!
!     Test the input arguments
!
  K = MIN(M, N)
  INFO = 0
  NB = mobbrmsd_ILAENV(1, 'DGEQRF', ' ', M, N, -1, -1)
  LQUERY = (LWORK == -1)
  if (M < 0) then
    INFO = -1
  else if (N < 0) then
    INFO = -2
  else if (LDA < MAX(1, M)) then
    INFO = -4
  else if (.not. LQUERY) then
    if (LWORK <= 0 .or. (M > 0 .and. LWORK < MAX(1, N))) INFO = -7
  end if
  if (INFO /= 0) then
    !CALL XERBLA( 'DGEQRF', -INFO )
    return
  else if (LQUERY) then
    if (K == 0) then
      LWKOPT = 1
    else
      LWKOPT = N * NB
    end if
    WORK(1) = LWKOPT
    return
  end if
!
! Quick return if possible
!
  if (K == 0) then
    WORK(1) = 1
    return
  end if
!
  NBMIN = 2
  NX = 0
  IWS = N
  if (NB > 1 .and. NB < K) then
!
! Determine when to cross over from blocked to unblocked code.
!
    NX = MAX(0, mobbrmsd_ILAENV(3, 'DGEQRF', ' ', M, N, -1, -1))
    if (NX < K) then
!
! Determine if workspace is large enough for blocked code.
!
      LDWORK = N
      IWS = LDWORK * NB
      if (LWORK < IWS) then
!
! Not enough workspace to use optimal NB:  reduce NB and
! determine the minimum value of NB.
!
        NB = LWORK / LDWORK
        NBMIN = MAX(2, mobbrmsd_ILAENV(2, 'DGEQRF', ' ', M, N, -1, -1))
      end if
    end if
  end if
!
  if (NB >= NBMIN .and. NB < K .and. NX < K) then
!
!        Use blocked code initially
!
    do I = 1, K - NX, NB
      IB = MIN(K - I + 1, NB)
!
!           Compute the QR factorization of the current block
!           A(i:m,i:i+ib-1)
!
      call mobbrmsd_DGEQR2(M - I + 1, IB, A(I, I), LDA, TAU(I), WORK, IINFO)
      if (I + IB <= N) then
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
        call mobbrmsd_DLARFT('Forward', 'Columnwise', M - I + 1, IB, &
            &       A(I, I), LDA, TAU(I), WORK, LDWORK)
!
!              Apply H**T to A(i:m,i+ib:n) from the left
!
        call mobbrmsd_DLARFB('Left', 'Transpose', 'Forward', &
            &       'Columnwise', M - I + 1, N - I - IB + 1, IB, &
            &       A(I, I), LDA, WORK, LDWORK, A(I, I + IB), &
            &       LDA, WORK(IB + 1), LDWORK)
      end if
    end do
  else
    I = 1
  end if
!
!     Use unblocked code to factor the last or only block.
!
  if (I <= K) call mobbrmsd_DGEQR2(M - I + 1, N - I + 1, A(I, I), LDA, TAU(I), WORK, IINFO)
!
  WORK(1) = IWS
  return
!
!     End of mobbrmsd_DGEQRF
!
end subroutine mobbrmsd_DGEQRF


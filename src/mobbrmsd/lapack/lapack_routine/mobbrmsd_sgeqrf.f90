!| mobbrmsd_DGEQRF computes a QR factorization of a real M-by-N matrix A:
!
!     A = Q * ( R ),
!             ( 0 )
!
!  where:
!
!     Q is a M-by-M orthogonal matrix;
!     R is an upper-triangular N-by-N matrix;
!     0 is a (M-N)-by-N zero matrix, if M > N.
!
!   The matrix Q is represented as a product of elementary reflectors
!
!      Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!   Each H(i) has the form
!
!      H(i) = I - tau * v * v**T
!
!   where \(\tau\) is a real scalar, and \(\mathbf{v}\) is a real vector with
!   \(v(1:i-1) = 0\) and \(v(i) = 1\); \(v(i+1:m)\) is stored on exit
!   in A(i+1:m,i), and tau in TAU(i).
!
!  Reference SGEQRF is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,--
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_SGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
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
!!  REAL array, dimension (LDA,N)
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
!!  REAL array, dimension (min(M,N))
!!  The scalar factors of the elementary reflectors (see Further
!!  Details).
!!
  real(RK), intent(out)   :: WORK(*)
!!  REAL array, dimension (MAX(1,LWORK))
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
!..external Functions..
!   include 'ilaenv.h'
! .. External Subroutines ..
!   include 'sgeqr2.h'
!   include 'slarfb.h'
!   include 'slarft.h'
! end interface
!
! Test the input arguments
!
  INFO = 0
  NB = mobbrmsd_ILAENV(1, 'SGEQRF', ' ', M, N, -1, -1)
  LWKOPT = N * NB
  WORK(1) = LWKOPT
  LQUERY = (LWORK == -1)
  if (M < 0) then
    INFO = -1
  else if (N < 0) then
    INFO = -2
  else if (LDA < MAX(1, M)) then
    INFO = -4
  else if (LWORK < MAX(1, N) .and. .not. LQUERY) then
    INFO = -7
  end if
  if (INFO /= 0) then
!   call XERBLA('SGEQRF', -INFO)
    return
  else if (LQUERY) then
    return
  end if
!
! Quick return if possible
!
  K = MIN(M, N)
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
!Determine when to cross over from blocked to unblocked code.
!
    NX = MAX(0, mobbrmsd_ILAENV(3, 'SGEQRF', ' ', M, N, -1, -1))
    if (NX < K) then
!
!Determine if workspace is large enough for blocked code.
!
      LDWORK = N
      IWS = LDWORK * NB
      if (LWORK < IWS) then
!
!Not enough workspace to use optimal NB:reduce NB and
!determine the minimum value of NB.
!
        NB = LWORK / LDWORK
        NBMIN = MAX(2, mobbrmsd_ILAENV(2, 'SGEQRF', ' ', M, N, -1, -1))
      end if
    end if
  end if
!
  if (NB >= NBMIN .and. NB < K .and. NX < K) then
!
! use blocked code initially
!
    do I = 1, K - NX, NB
      IB = MIN(K - I + 1, NB)
!
! Compute the QR factorization of the current block
! A(i:m, i:i + ib - 1)
!
      call mobbrmsd_SGEQR2(M - I + 1, IB, A(I, I), LDA, TAU(I), WORK, IINFO)
      if (I + IB <= N) then
!
! Form the triangular factor of the block reflector
! H = H(i) H(i + 1) ...H(i + ib - 1)
!
        call mobbrmsd_SLARFT('Forward', 'Columnwise', M - I + 1, IB, A(I, I), LDA, TAU(I), WORK, LDWORK)
!
! Apply H**T to A(i:m, i + ib:n) from the left
!
        call mobbrmsd_SLARFB('Left', 'Transpose', 'Forward', &
            &                'Columnwise', M - I + 1, N - I - IB + 1, IB, &
            &                A(I, I), LDA, WORK, LDWORK, A(I, I + IB), &
            &                LDA, WORK(IB + 1), LDWORK)
      end if
    end do
  else
    I = 1
  end if
!
! use unblocked code to factor the last or only block.
!
  if (I <= K) call mobbrmsd_SGEQR2(M - I + 1, N - I + 1, A(I, I), LDA, TAU(I), WORK, IINFO)
!
  WORK(1) = IWS
  return
!
! end of mobbrmsd_SGEQRF
!
end

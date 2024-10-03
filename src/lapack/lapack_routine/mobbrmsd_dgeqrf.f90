!> \brief \b mobbrmsd_DGEQRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DGEQRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dgeqrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dgeqrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dgeqrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       pure subroutine mobbrmsd_DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       real(RK)           ::   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DGEQRF computes a QR factorization of a real M-by-N matrix A:
!>
!>    A = Q * ( R ),
!>            ( 0 )
!>
!> where:
!>
!>    Q is a M-by-M orthogonal matrix;
!>    R is an upper-triangular N-by-N matrix;
!>    0 is a (M-N)-by-N zero matrix, if M > N.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is real(RK)           :: array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!>          upper triangular if m >= n); the elements below the diagonal,
!>          with the array TAU, represent the orthogonal matrix Q as a
!>          product of min(m,n) elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is real(RK)           :: array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is real(RK)           :: array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          LWORK >= 1, if MIN(M,N) = 0, and LWORK >= N, otherwise.
!>          For optimum performance LWORK >= N*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup doubleGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!>  and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
pure subroutine mobbrmsd_DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)     :: LDA, LWORK, M, N
  integer, intent(out)    :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: A(LDA, *)
  real(RK), intent(out)   :: TAU(*), WORK(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  logical :: LQUERY
  integer :: I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
  intrinsic :: MAX, MIN
!     ..
! interface
!     .. External Subroutines ..
!   include 'dgeqr2.h'
!   include 'dlarfb.h'
!   include 'dlarft.h'
!   !include 'xerbla.h'
!     .. External Functions ..
!   include 'ilaenv.h'
! end interface
!     ..
!     .. Executable Statements ..
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

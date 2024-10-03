!> \brief \b mobbrmsd_SGELQF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_SGELQF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelqf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelqf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelqf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_SGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_SGELQF computes an LQ factorization of a real M-by-N matrix A:
!>
!>    A = ( L 0 ) *  Q
!>
!> where:
!>
!>    Q is a N-by-N orthogonal matrix;
!>    L is a lower-triangular M-by-M matrix;
!>    0 is a M-by-(N-M) zero matrix, if M < N.
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and below the diagonal of the array
!>          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
!>          lower triangular if m <= n); the elements above the diagonal,
!>          with the array TAU, represent the orthogonal matrix Q as a
!>          product of elementary reflectors (see Further Details).
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
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,M).
!>          For optimum performance LWORK >= M*NB, where NB is the
!>          optimal blocksize.
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
!> \date November 2019
!
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(k) . . . H(2) H(1), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
!>  and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
pure subroutine mobbrmsd_SGELQF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
  integer, intent(in)  :: LDA, LWORK, M, N
  integer, intent(out) :: INFO
!..
!..Array Arguments..
  real(RK), intent(inout)  :: A(LDA, *)
  real(RK), intent(out)    :: TAU(*), WORK(*)
!..
!
!  =====================================================================
!
!..Local Scalars..
  logical :: LQUERY
  integer :: I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX
!..
! interface
! .. External Functions ..
!   include 'ilaenv.h'
! .. External Subroutines ..
!   include 'sgelq2.h'
!   include 'slarfb.h'
!   include 'slarft.h'
!   include 'xerbla.h'
! end interface
!..
!..intrinsic Functions..
  intrinsic :: MAX, MIN
!..
!..Executable Statements..
!
!Test the input arguments
!
  INFO = 0
  NB = mobbrmsd_ILAENV(1, 'mobbrmsd_SGELQF', ' ', M, N, -1, -1)
  LWKOPT = M * NB
  WORK(1) = LWKOPT
  LQUERY = (LWORK == -1)
  if (M < 0) then
    INFO = -1
  else if (N < 0) then
    INFO = -2
  else if (LDA < MAX(1, M)) then
    INFO = -4
  else if (LWORK < MAX(1, M) .and. .not. LQUERY) then
    INFO = -7
  end if
  if (INFO /= 0) then
!   call XERBLA('mobbrmsd_SGELQF', -INFO)
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
  IWS = M
  if (NB > 1 .and. NB < K) then
    !
    !Determine when to cross over from blocked to unblocked code.
    !
    NX = MAX(0, mobbrmsd_ILAENV(3, 'mobbrmsd_SGELQF', ' ', M, N, -1, -1))
    if (NX < K) then
      !
      !Determine if workspace is large enough for blocked code.
      !
      LDWORK = M
      IWS = LDWORK * NB
      if (LWORK < IWS) then
        !
        !Not enough workspace to use optimal NB:reduce NB and
        !determine the minimum value of NB.
        !
        NB = LWORK / LDWORK
        NBMIN = MAX(2, mobbrmsd_ILAENV(2, 'mobbrmsd_SGELQF', ' ', M, N, -1, -1))
      end if
    end if
  end if
!
  if (NB >= NBMIN .and. NB < K .and. NX < K) then
    !
    !use blocked code initially
    !
    do I = 1, K - NX, NB
      IB = MIN(K - I + 1, NB)
      !
      !Compute the LQ factorization of the current block
      !A(i:i + ib - 1, i:n)
      !
      call mobbrmsd_SGELQ2(IB, N - I + 1, A(I, I), LDA, TAU(I), WORK, IINFO)
      if (I + IB <= M) then
        !
        !Form the triangular factor of the block reflector
        !H = H(i) H(i + 1) ...H(i + ib - 1)
        !
        call mobbrmsd_SLARFT('Forward', 'Rowwise', N - I + 1, IB, A(I, I), LDA, TAU(I), WORK, LDWORK)
        !
        !Apply H to A(i + ib:m, i:n) from the right
        !
        call mobbrmsd_SLARFB('Right', 'No transpose', 'Forward', &
            &       'Rowwise', M - I - IB + 1, N - I + 1, IB, A(I, I), &
            &       LDA, WORK, LDWORK, A(I + IB, I), LDA, &
            &       WORK(IB + 1), LDWORK)
      end if
    end do
  else
    I = 1
  end if
  !
  !use unblocked code to factor the last or only block.
  !
  if (I <= K) call mobbrmsd_SGELQ2(M - I + 1, N - I + 1, A(I, I), LDA, TAU(I), WORK, IINFO)
  !
  WORK(1) = IWS
  return
  !
  !end of mobbrmsd_SGELQF
  !
end

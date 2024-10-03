!> \brief \b SORGBR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SORGBR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgbr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgbr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgbr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          VECT
!       INTEGER            INFO, K, LDA, LWORK, M, N
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
!> SORGBR generates one of the real orthogonal matrices Q or P**T
!> determined by SGEBRD when reducing a real matrix A to bidiagonal
!> form: A = Q * B * P**T.  Q and P**T are defined as products of
!> elementary reflectors H(i) or G(i) respectively.
!>
!> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
!> is of order M:
!> if m >= k, Q = H(1) H(2) . . . H(k) and SORGBR returns the first n
!> columns of Q, where m >= n >= k;
!> if m < k, Q = H(1) H(2) . . . H(m-1) and SORGBR returns Q as an
!> M-by-M matrix.
!>
!> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
!> is of order N:
!> if k < n, P**T = G(k) . . . G(2) G(1) and SORGBR returns the first m
!> rows of P**T, where n >= m >= k;
!> if k >= n, P**T = G(n-1) . . . G(2) G(1) and SORGBR returns P**T as
!> an N-by-N matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          Specifies whether the matrix Q or the matrix P**T is
!>          required, as defined in the transformation applied by SGEBRD:
!>          = 'Q':  generate Q;
!>          = 'P':  generate P**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q or P**T to be returned.
!>          M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q or P**T to be returned.
!>          N >= 0.
!>          If VECT = 'Q', M >= N >= min(M,K);
!>          if VECT = 'P', N >= M >= min(N,K).
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          If VECT = 'Q', the number of columns in the original M-by-K
!>          matrix reduced by SGEBRD.
!>          If VECT = 'P', the number of rows in the original K-by-N
!>          matrix reduced by SGEBRD.
!>          K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by SGEBRD.
!>          On exit, the M-by-N matrix Q or P**T.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension
!>                                (min(M,K)) if VECT = 'Q'
!>                                (min(N,K)) if VECT = 'P'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i) or G(i), which determines Q or P**T, as
!>          returned by SGEBRD in its array argument TAUQ or TAUP.
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
!>          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
!>          For optimum performance LWORK >= min(M,N)*NB, where NB
!>          is the optimal blocksize.
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
!> \date April 2012
!
!> \ingroup realGBcomputational
!
!  =====================================================================
pure subroutine SORGBR(VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
  implicit none
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
  character, intent(in) :: VECT
  integer, intent(in)   :: K, LDA, LWORK, M, N
  integer, intent(out)  :: INFO
!..
!..Array Arguments..
  real, intent(in)      ::  TAU(*)
  real, intent(inout)   ::  A(LDA, *)
  real, intent(out)     ::  WORK(*)
!..
!
!  =====================================================================
!..
!..Local Scalars..
  logical :: LQUERY, WANTQ
  integer :: I, IINFO, J, LWKOPT, MN
!
!..Parameters..
  real, parameter :: ONE = 1.0E+0
  real, parameter :: ZERO = 0.0E+0
!..
!..Local Scalars..
  integer :: I, J, L
!..
  interface
!..external Functions..
    include 'lsame.h'
!..external Subroutines..
    include 'sorglq.h'
    include 'sorgqr.h'
  end interface
!..
!..intrinsic Functions..
  intrinsic :: MAX, MIN
!..
!..Executable Statements..
!
!Test the input arguments
!
  INFO = 0
  WANTQ = LSAME(VECT, 'Q')
  MN = MIN(M, N)
  LQUERY = (LWORK == -1)
  if (.not. WANTQ .and. .not. LSAME(VECT, 'P')) then
    INFO = -1
  else if (M < 0) then
    INFO = -2
    else if (N < 0 .or. (WANTQ .and. (N > M .or. N < MIN(M,
    $K))) .or. (.not. WANTQ .and. (M > N .or. M <
    $MIN(N, K)))) then
    INFO = -3
  else if (K < 0) then
    INFO = -4
  else if (LDA < MAX(1, M)) then
    INFO = -6
  else if (LWORK < MAX(1, MN) .and. .not. LQUERY) then
    INFO = -9
  end if
!
  if (INFO == 0) then
    WORK(1) = 1
    if (WANTQ) then
      if (M >= K) then
        call SORGQR(M, N, K, A, LDA, TAU, WORK, -1, IINFO)
      else
        if (M > 1) then
          call SORGQR(M - 1, M - 1, M - 1, A(2, 2), LDA, TAU, WORK, -1, IINFO)
        end if
      end if
    else
      if (K < N) then
        call SORGLQ(M, N, K, A, LDA, TAU, WORK, -1, IINFO)
      else
        if (N > 1) then
          call SORGLQ(N - 1, N - 1, N - 1, A(2, 2), LDA, TAU, WORK, -1, IINFO)
        end if
      end if
    end if
    LWKOPT = WORK(1)
    LWKOPT = MAX(LWKOPT, MN)
  end if
!
  if (INFO /= 0) then
!   call XERBLA('SORGBR', -INFO)
    return
  else if (LQUERY) then
    WORK(1) = LWKOPT
    return
  end if
!
! Quick return if possible
!
  if (M == 0 .or. N == 0) then
    WORK(1) = 1
    return
  end if
!
  if (WANTQ) then
!
! Form Q, determined by a call to SGEBRD to reduce an m - by - k matrix
!
    if (M >= K) then
!
! if m >= k, assume m >= n >= k
!
      call SORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO)
!
    else
!
! if m < k, assume m = n
!
! Shift the vectors which define the elementary reflectors one
! column to the right, and set the first row and column of Q
! to those of the unit matrix
!
      do J = M, 2, -1
        A(1, J) = ZERO
        do I = J + 1, M
          A(I, J) = A(I, J - 1)
        end do
      end do
      A(1, 1) = ONE
      do I = 2, M
        A(I, 1) = ZERO
      end do
      if (M > 1) then
!
! Form Q(2:m, 2:m)
!
        call SORGQR(M - 1, M - 1, M - 1, A(2, 2), LDA, TAU, WORK, LWORK, IINFO)
      end if
    end if
  else
!
! Form P**T, determined by a call to SGEBRD to reduce a k - by - n matrix
!
    if (K < N) then
!
!if k < n, assume k <= m <= n
!
      call SORGLQ(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO)
!
    else
!
!if k >= n, assume m = n
!
!Shift the vectors which define the elementary reflectors one
!row downward, and set the first row and column of P**T to
!those of the unit matrix
!
      A(1, 1) = ONE
      do I = 2, N
        A(I, 1) = ZERO
      end do
      do J = 2, N
        do I = J - 1, 2, -1
          A(I, J) = A(I - 1, J)
        end do
        A(1, J) = ZERO
      end do
      if (N > 1) then
!
! Form P**T(2:n, 2:n)
!
        call SORGLQ(N - 1, N - 1, N - 1, A(2, 2), LDA, TAU, WORK, LWORK, IINFO)
      end if
    end if
  end if
  WORK(1) = LWKOPT
  return
!
! end of SORGBR
!
end

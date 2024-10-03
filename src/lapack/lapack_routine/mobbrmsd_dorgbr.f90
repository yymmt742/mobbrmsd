!> \brief \b mobbrmsd_DORGBR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DORGBR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dorgbr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dorgbr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dorgbr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          VECT
!       INTEGER            INFO, K, LDA, LWORK, M, N
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
!> mobbrmsd_DORGBR generates one of the real orthogonal matrices Q or P**T
!> determined by mobbrmsd_DGEBRD when reducing a real matrix A to bidiagonal
!> form: A = Q * B * P**T.  Q and P**T are defined as products of
!> elementary reflectors H(i) or G(i) respectively.
!>
!> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
!> is of order M:
!> if m >= k, Q = H(1) H(2) . . . H(k) and mobbrmsd_DORGBR returns the first n
!> columns of Q, where m >= n >= k;
!> if m < k, Q = H(1) H(2) . . . H(m-1) and mobbrmsd_DORGBR returns Q as an
!> M-by-M matrix.
!>
!> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
!> is of order N:
!> if k < n, P**T = G(k) . . . G(2) G(1) and mobbrmsd_DORGBR returns the first m
!> rows of P**T, where n >= m >= k;
!> if k >= n, P**T = G(n-1) . . . G(2) G(1) and mobbrmsd_DORGBR returns P**T as
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
!>          required, as defined in the transformation applied by mobbrmsd_DGEBRD:
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
!>          matrix reduced by mobbrmsd_DGEBRD.
!>          If VECT = 'P', the number of rows in the original K-by-N
!>          matrix reduced by mobbrmsd_DGEBRD.
!>          K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is real(RK)           :: array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by mobbrmsd_DGEBRD.
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
!>          TAU is real(RK)           :: array, dimension
!>                                (min(M,K)) if VECT = 'Q'
!>                                (min(N,K)) if VECT = 'P'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i) or G(i), which determines Q or P**T, as
!>          returned by mobbrmsd_DGEBRD in its array argument TAUQ or TAUP.
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
!> \ingroup doubleGBcomputational
!
!  =====================================================================
pure subroutine mobbrmsd_DORGBR(VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: VECT
  integer, intent(in)      :: K, LDA, LWORK, M, N
  integer, intent(out)     :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)     :: TAU(*)
  real(RK), intent(inout)  :: A(LDA, *)
  real(RK), intent(out)    :: WORK(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
  logical                 :: LQUERY, WANTQ
  integer                 :: I, IINFO, J, LWKOPT, MN
!     .. Intrinsic Functions ..
  intrinsic                :: MAX, MIN
!     .. Parameters ..
! real(RK), parameter      :: ZERO = 0.0_RK
! real(RK), parameter      :: ONE = 1.0_RK
!     ..
! interface
!     .. External Subroutines ..
!   include 'dorglq.h'
!   include 'dorgqr.h'
!   !include 'xerbla.h'
!     .. External Functions ..
!   include 'lsame.h'
! end interface
!     ..
!  =====================================================================
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  WANTQ = mobbrmsd_LSAME(VECT, 'Q')
  MN = MIN(M, N)
  LQUERY = (LWORK == -1)
  if (.not. WANTQ .and. .not. mobbrmsd_LSAME(VECT, 'P')) then
    INFO = -1
  else if (M < 0) then
    INFO = -2
  else if (N < 0 .or. (WANTQ .and. (N > M .or. N < MIN(M, K))) &
 &         .or. (.not. WANTQ .and. (M > N .or. M < MIN(N, K)))) then
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
        call mobbrmsd_DORGQR(M, N, K, A, LDA, TAU, WORK, -1, IINFO)
      else
        if (M > 1) then
          call mobbrmsd_DORGQR(M - 1, M - 1, M - 1, A, LDA, TAU, WORK, -1, IINFO)
        end if
      end if
    else
      if (K < N) then
        call mobbrmsd_DORGLQ(M, N, K, A, LDA, TAU, WORK, -1, IINFO)
      else
        if (N > 1) then
          call mobbrmsd_DORGLQ(N - 1, N - 1, N - 1, A, LDA, TAU, WORK, -1, IINFO)
        end if
      end if
    end if
    LWKOPT = WORK(1)
    LWKOPT = MAX(LWKOPT, MN)
  end if
!
  if (INFO /= 0) then
    !CALL XERBLA( 'mobbrmsd_DORGBR', -INFO )
    return
  else if (LQUERY) then
    WORK(1) = LWKOPT
    return
  end if
!
!     Quick return if possible
!
  if (M == 0 .or. N == 0) then
    WORK(1) = 1
    return
  end if
!
  if (WANTQ) then
!
!        Form Q, determined by a call to mobbrmsd_DGEBRD to reduce an m-by-k
!        matrix
!
    if (M >= K) then
!
!           If m >= k, assume m >= n >= k
!
      call mobbrmsd_DORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO)
!
    else
!
!           If m < k, assume m = n
!
!           Shift the vectors which define the elementary reflectors one
!           column to the right, and set the first row and column of Q
!           to those of the unit matrix
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
!             Form Q(2:m,2:m)
!
        call mobbrmsd_DORGQR(M - 1, M - 1, M - 1, A(2, 2), LDA, TAU, WORK, LWORK, IINFO)
      end if
    end if
  else
!
!        Form P**T, determined by a call to mobbrmsd_DGEBRD to reduce a k-by-n
!        matrix
!
    if (K < N) then
!
!           If k < n, assume k <= m <= n
!
      call mobbrmsd_DORGLQ(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO)
!
    else
!
!           If k >= n, assume m = n
!
!           Shift the vectors which define the elementary reflectors one
!           row downward, and set the first row and column of P**T to
!           those of the unit matrix
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
!            Form P**T(2:n,2:n)
!
        call mobbrmsd_DORGLQ(N - 1, N - 1, N - 1, A(2, 2), LDA, TAU, WORK, LWORK, IINFO)
      end if
    end if
  end if
  WORK(1) = LWKOPT
  return
!
!     End of mobbrmsd_DORGBR
!
end subroutine mobbrmsd_DORGBR

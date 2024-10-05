!| mobbrmsd_DORGBR generates one of the real orthogonal matrices
!  \(\mathbf{Q}\) or \(\mathbf{P}^\top\).
!  determined by mobbrmsd_DGEBRD when reducing a real matrix A to bidiagonal
!  form: \(\mathbf{A}=\mathbf{Q}\mathbf{B}\mathbf{P}^\top\).
!  \(\mathbf{Q}\) and \(\mathbf{P}^\top\) are defined as products of
!  elementary reflectors H(i) or G(i) respectively.
!
!  If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
!  is of order M:
!
!  If m >= k, Q = H(1) H(2) . . . H(k) and mobbrmsd_DORGBR returns the first n
!  columns of Q, where m >= n >= k.
!
!  if m < k, Q = H(1) H(2) . . . H(m-1) and mobbrmsd_DORGBR returns Q as an
!  M-by-M matrix.
!
!  If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
!  is of order N:
!
!  if k < n, P**T = G(k) . . . G(2) G(1) and mobbrmsd_DORGBR returns the first m
!  rows of P**T, where n >= m >= k;
!
!  if k >= n, P**T = G(n-1) . . . G(2) G(1) and mobbrmsd_DORGBR returns P**T as
!  an N-by-N matrix.
!
!  Reference DORGBR is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- LAPACK driver routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!     \date April 2012
!
pure subroutine mobbrmsd_DORGBR(VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
  implicit none
!
  character, intent(in) :: VECT
!!          Specifies whether the matrix \(\mathbf{Q}\) or the matrix \(\mathbf{P}^\top\) is
!!          required, as defined in the transformation applied by mobbrmsd_DGEBRD:
!!
!!          = 'Q':  generate \(\mathbf{Q}\).
!!
!!          = 'P':  generate \(\mathbf{P}^\top\).
!!
!! \endverbatim
  integer, intent(in)      :: M
!!          The number of rows of the matrix \(\mathbf{Q}\) or \(\mathbf{P}^\top\) to be returned.
!!          M >= 0.
!!
  integer, intent(in)      :: N
!!          The number of columns of the matrix Q or P**T to be returned.
!!          N >= 0.
!!
!!          If VECT = 'Q', M >= N >= min(M,K);
!!
!!          if VECT = 'P', N >= M >= min(N,K).
!!
  integer, intent(in)      :: K
!!          If VECT = 'Q', the number of columns in the original M-by-K
!!          matrix reduced by mobbrmsd_DGEBRD.
!!
!!          If VECT = 'P', the number of rows in the original K-by-N
!!          matrix reduced by mobbrmsd_DGEBRD.
!!
!!          K >= 0.
!!
  real(RK), intent(in)     :: TAU(*)
!!          TAU is real(RK)           :: array, dimension
!!
!!                                (min(M,K)) if VECT = 'Q'
!!
!!                                (min(N,K)) if VECT = 'P'
!!
!!          TAU(i) must contain the scalar factor of the elementary
!!          reflector H(i) or G(i), which determines \(\mathbf{Q}\) or \(\mathbf{P}^\top\), as
!!          returned by mobbrmsd_DGEBRD in its array argument TAUQ or TAUP.
!!
  integer, intent(in)      :: LDA
!!          The leading dimension of the array A. LDA >= max(1,M).
!!
  real(RK), intent(inout)  :: A(LDA, *)
!!          A is real(RK)           :: array, dimension (LDA,N)
!!
!!          On entry, the vectors which define the elementary reflectors,
!!          as returned by mobbrmsd_DGEBRD.
!!
!!          On exit, the M-by-N matrix \(\mathbf{Q}\) or \(\mathbf{P}^\top\).
!!
  integer, intent(in)      :: LWORK
!!          LWORK is INTEGER
!!
!!          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
!!          For optimum performance LWORK >= min(M,N)*NB, where NB
!!          is the optimal blocksize.
!!
!!          If LWORK = -1, then a workspace query is assumed; the routine
!!          only calculates the optimal size of the WORK array, returns
!!          this value as the first entry of the WORK array, and no error
!!          message related to LWORK is issued by XERBLA.
!!
  real(RK), intent(out)    :: WORK(*)
!!          WORK is real(RK)           :: array, dimension (MAX(1,LWORK))
!!
!!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  integer, intent(out)     :: INFO
!!          = 0:  successful exit
!!
!!          < 0:  if INFO = -i, the i-th argument had an illegal value
!!
  logical                 :: LQUERY, WANTQ
  integer                 :: I, IINFO, J, LWKOPT, MN
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
    !CALL XERBLA( 'DORGBR', -INFO )
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

!| generates a product of \( k \) elementary reflectors
!
!  mobbrmsd_DORGQR generates an \( m \)-by-\( n \) real matrix \( Q \)
!  with orthonormal columns,
!  which is defined as the first \( n \) columns of a product of
!  \( k \) K elementary reflectors of order \( m \)
!
!  \[
!     Q  =  H _ 1 H _ 2 \cdots H _ k
!  \]
!
!  as returned by mobbrmsd_DGEQRF.
!
!  Reference DORGQR is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
  implicit none
  integer, intent(in)     :: M
!!  The number of rows of the matrix Q. M >= 0.
!!
  integer, intent(in)     :: N
!!  The number of columns of the matrix Q. M >= N >= 0.
!!
  integer, intent(in)     :: K
!!  The number of elementary reflectors whose product defines the
!!  matrix Q. N >= K >= 0.
!!
  integer, intent(in)     :: LDA
!!  The first dimension of the array A. LDA >= max(1,M).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!  DOUBLE PRECISION array, dimension (LDA,N)
!!
!!  On entry, the i-th column must contain the vector which
!!  defines the elementary reflector H(i), for i = 1,2,...,k, as
!!  returned by mobbrmsd_DGEQRF in the first k columns of its array
!!  argument A.
!!
!!  On exit, the M-by-N matrix Q.
!!
  real(RK), intent(in)    :: TAU(*)
!!  DOUBLE PRECISION array, dimension (K)
!!
!!  TAU(i) must contain the scalar factor of the elementary
!!  reflector H(i), as returned by mobbrmsd_DGEQRF.
!!
  real(RK), intent(out)   :: WORK(*)
!!  DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!  On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!
  integer, intent(in)     :: LWORK
!!  The dimension of the array WORK. LWORK >= max(1,N).
!!
!!  For optimum performance LWORK >= N*NB, where NB is the
!!  optimal blocksize.
!!
  integer, intent(out)    :: INFO
!!  If LWORK = -1, then a workspace query is assumed; the routine
!!  only calculates the optimal size of the WORK array, returns
!!  this value as the first entry of the WORK array.
!!
!!  = 0:  successful exit
!!
!!  < 0:  if INFO = -i, the i-th argument has an illegal value
!!
  logical                 :: LQUERY
  integer                 :: I, IB, IINFO, IWS, J, KI, KK, L, &
 &                           LDWORK, LWKOPT, NB, NBMIN, NX
  intrinsic               :: MAX, MIN
! interface
!   include 'dlarfb.h'
!   include 'dlarft.h'
!   include 'dorg2r.h'
!   include 'ilaenv.h'
! end interface
!
! Test the input arguments
!
  INFO = 0
  NB = mobbrmsd_ILAENV(1, 'DORGQR', ' ', M, N, K, -1)
  LWKOPT = MAX(1, N) * NB
  WORK(1) = LWKOPT
  LQUERY = (LWORK == -1)
  if (M < 0) then
    INFO = -1
  else if (N < 0 .or. N > M) then
    INFO = -2
  else if (K < 0 .or. K > N) then
    INFO = -3
  else if (LDA < MAX(1, M)) then
    INFO = -5
  else if (LWORK < MAX(1, N) .and. .not. LQUERY) then
    INFO = -8
  end if
  if (INFO /= 0) then
    !CALL XERBLA( 'DORGQR', -INFO )
    return
  else if (LQUERY) then
    return
  end if
!
!     Quick return if possible
!
  if (N <= 0) then
    WORK(1) = 1
    return
  end if
!
  NBMIN = 2
  NX = 0
  IWS = N
  if (NB > 1 .and. NB < K) then
!
!        Determine when to cross over from blocked to unblocked code.
!
    NX = MAX(0, mobbrmsd_ILAENV(3, 'DORGQR', ' ', M, N, K, -1))
    if (NX < K) then
!
!           Determine if workspace is large enough for blocked code.
!
      LDWORK = N
      IWS = LDWORK * NB
      if (LWORK < IWS) then
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
        NB = LWORK / LDWORK
        NBMIN = MAX(2, mobbrmsd_ILAENV(2, 'DORGQR', ' ', M, N, K, -1))
      end if
    end if
  end if
!
  if (NB >= NBMIN .and. NB < K .and. NX < K) then
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
    KI = ((K - NX - 1) / NB) * NB
    KK = MIN(K, KI + NB)
!
!        Set A(1:kk,kk+1:n) to zero.
!
    do J = KK + 1, N
      do I = 1, KK
        A(I, J) = ZERO
      end do
    end do
  else
    KK = 0
  end if
!
!     Use unblocked code for the last or only block.
!
  if (KK < N) &
 &   call mobbrmsd_DORG2R(M - KK, N - KK, K - KK, A(KK + 1, KK + 1), LDA, &
 &                TAU(KK + 1), WORK, IINFO)
!
  if (KK > 0) then
!
!        Use blocked code
!
    do I = KI + 1, 1, -NB
      IB = MIN(NB, K - I + 1)
      if (I + IB <= N) then
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
        call mobbrmsd_DLARFT('Forward', 'Columnwise', M - I + 1, IB, &
       &             A(I, I), LDA, TAU(I), WORK, LDWORK)
!
!              Apply H to A(i:m,i+ib:n) from the left
!
        call mobbrmsd_DLARFB('Left', 'No transpose', 'Forward', &
    &                'Columnwise', M - I + 1, N - I - IB + 1, IB, &
    &                A(I, I), LDA, WORK, LDWORK, A(I, I + IB), &
    &                LDA, WORK(IB + 1), LDWORK)
      end if
!
!           Apply H to rows i:m of current block
!
      call mobbrmsd_DORG2R(M - I + 1, IB, IB, A(I, I), LDA, TAU(I), WORK, IINFO)
!
!           Set rows 1:i-1 of current block to zero
!
      do J = I, I + IB - 1
        do L = 1, I - 1
          A(L, J) = ZERO
        end do
      end do
    end do
  end if
!
  WORK(1) = IWS
  return
!
!     End of mobbrmsd_DORGQR
!
end subroutine mobbrmsd_DORGQR


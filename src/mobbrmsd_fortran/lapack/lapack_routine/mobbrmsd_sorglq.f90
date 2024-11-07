!| generates an \( Q \in \mathbb{R} ^ {m \times n} \) with orthonormal rows.
!
!  mobbrmsd_SORGLQ generates an \( M \)-by-\( N \) real matrix \( Q \)
!  with orthonormal rows, which is defined as the first-\( M \) rows
!  of a product of \( K \) elementary reflectors of order \( N \).
!
!  \[
!     Q =  H _ k \cdots H _ 2 H _ 1
!  \]
!
!  as returned by mobbrmsd_SGELQF.
!
!  Reference SORGLQ is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_SORGLQ(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
  implicit none
  integer, intent(in)      :: M
!! The number of rows of the matrix Q. M >= 0.
!!
  integer, intent(in)      :: N
!!  The number of columns of the matrix Q. N >= M.
!!
  integer, intent(in)      :: K
!!  The number of elementary reflectors whose product defines the
!!  matrix Q. M >= K >= 0.
!!
  integer, intent(in)      :: LDA
!!  The first dimension of the array A. LDA >= max(1,M).
!!
  real(RK), intent(inout)  :: A(LDA, *)
!!  DOUBLE PRECISION array, dimension (LDA,N)
!!
!!  On entry, the i-th row must contain the vector which defines
!!  the elementary reflector H(i), for i = 1,2,...,k, as returned
!!  by mobbrmsd_DGELQF in the first k rows of its array argument A.
!!
!!  On exit, the M-by-N matrix Q.
!!
  real(RK), intent(in)     :: TAU(*)
!!  DOUBLE PRECISION array, dimension (K)
!!
!!  TAU(i) must contain the scalar factor of the elementary
!!  reflector H(i), as returned by mobbrmsd_DGELQF.
!!
  real(RK), intent(out)    :: WORK(*)
!!  DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!
!!  On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!
  integer, intent(in)      :: LWORK
!!  The dimension of the array WORK. LWORK >= max(1,M).
!!
!!  For optimum performance LWORK >= M*NB, where NB is
!!  the optimal blocksize.
!!
!!  If LWORK = -1, then a workspace query is assumed; the routine
!!  only calculates the optimal size of the WORK array, returns
!!  this value as the first entry of the WORK array, and no error
!!  message related to LWORK is issued by XERBLA.
!!
  integer, intent(out)     :: INFO
!!  = 0:  successful exit
!!
!!  < 0:  if INFO = -i, the i-th argument has an illegal value
!!
  logical :: LQUERY
  integer :: I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX
  intrinsic :: MAX, MIN
!
! Test the input arguments
!
  INFO = 0
  NB = mobbrmsd_ILAENV(1, 'SORGLQ', ' ', M, N, K, -1)
  LWKOPT = MAX(1, M) * NB
  WORK(1) = LWKOPT
  LQUERY = (LWORK == -1)
  if (M < 0) then
    INFO = -1
  else if (N < M) then
    INFO = -2
  else if (K < 0 .or. K > M) then
    INFO = -3
  else if (LDA < MAX(1, M)) then
    INFO = -5
  else if (LWORK < MAX(1, M) .and. .not. LQUERY) then
    INFO = -8
  end if
  if (INFO /= 0) then
    return
  else if (LQUERY) then
    return
  end if
!
! Quick return if possible
!
  if (M <= 0) then
    WORK(1) = 1
    return
  end if
!
  NBMIN = 2
  NX = 0
  IWS = M
  if (NB > 1 .and. NB < K) then
!
!   Determine when to cross over from blocked to unblocked code.
!
    NX = MAX(0, mobbrmsd_ILAENV(3, 'SORGLQ', ' ', M, N, K, -1))
    if (NX < K) then
!
!     Determine if workspace is large enough for blocked code.
!
      LDWORK = M
      IWS = LDWORK * NB
      if (LWORK < IWS) then
!
!       Not enough workspace to use optimal NB:reduce NB and
!       determine the minimum value of NB.
!
        NB = LWORK / LDWORK
        NBMIN = MAX(2, mobbrmsd_ILAENV(2, 'SORGLQ', ' ', M, N, K, -1))
      end if
    end if
  end if
!
  if (NB >= NBMIN .and. NB < K .and. NX < K) then
!
! use blocked code after the last block.
! The first kk rows are handled by the block method.
!
    KI = ((K - NX - 1) / NB) * NB
    KK = MIN(K, KI + NB)
!
! Set A(kk + 1:m, 1:kk) to zero.
!
    do J = 1, KK
      do I = KK + 1, M
        A(I, J) = ZERO
      end do
    end do
  else
    KK = 0
  end if
!
! use unblocked code for the last or only block.
!
  if (KK < M) call mobbrmsd_SORGL2(M - KK, N - KK, K - KK, A(KK + 1, KK + 1), LDA, TAU(KK + 1), WORK, IINFO)
  !
  if (KK > 0) then
!
! use blocked code
!
    do I = KI + 1, 1, -NB
      IB = MIN(NB, K - I + 1)
      if (I + IB <= M) then
        !
        !Form the triangular factor of the block reflector
        !H = H(i) H(i + 1) ...H(i + ib - 1)
        !
        call mobbrmsd_SLARFT('Forward', 'Rowwise', N - I + 1, IB, A(I, I), LDA, TAU(I), WORK, LDWORK)
        !
        !Apply H**T to A(i + ib:m, i:n) from the right
        !
        call mobbrmsd_SLARFB('Right', 'Transpose', 'Forward', 'Rowwise', &
            &       M - I - IB + 1, N - I + 1, IB, A(I, I), LDA, WORK, &
            &       LDWORK, A(I + IB, I), LDA, WORK(IB + 1), LDWORK)
      end if
!
! Apply H**T to columns i:n of current block
!
      call mobbrmsd_SORGL2(IB, N - I + 1, IB, A(I, I), LDA, TAU(I), WORK, IINFO)
!
! Set columns 1:i - 1 of current block to zero
!
      do J = 1, I - 1
        do L = I, I + IB - 1
          A(L, J) = ZERO
        end do
      end do
    end do
  end if
!
  WORK(1) = IWS
  return
!
! end of mobbrmsd_SORGLQ
!
end


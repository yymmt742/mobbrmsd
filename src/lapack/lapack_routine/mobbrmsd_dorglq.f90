!> \brief \b mobbrmsd_DORGLQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DORGLQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dorglq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dorglq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dorglq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
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
!> mobbrmsd_DORGLQ generates an M-by-N real matrix Q with orthonormal rows,
!> which is defined as the first M rows of a product of K elementary
!> reflectors of order N
!>
!>       Q  =  H(k) . . . H(2) H(1)
!>
!> as returned by mobbrmsd_DGELQF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. N >= M.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. M >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is real(RK)           :: array, dimension (LDA,N)
!>          On entry, the i-th row must contain the vector which defines
!>          the elementary reflector H(i), for i = 1,2,...,k, as returned
!>          by mobbrmsd_DGELQF in the first k rows of its array argument A.
!>          On exit, the M-by-N matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is real(RK)           :: array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by mobbrmsd_DGELQF.
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
!>          The dimension of the array WORK. LWORK >= max(1,M).
!>          For optimum performance LWORK >= M*NB, where NB is
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
!>          < 0:  if INFO = -i, the i-th argument has an illegal value
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
pure subroutine mobbrmsd_DORGLQ(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)      :: K, LDA, LWORK, M, N
  integer, intent(out)     :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)     :: TAU(*)
  real(RK), intent(inout)  :: A(LDA, *)
  real(RK), intent(out)    :: WORK(*)
!     ..
!  =====================================================================
!     ..
!     .. Local Scalars ..
  logical                 :: LQUERY
  integer                 :: I, IB, IINFO, IWS, J, KI, KK, L, &
 &                           LDWORK, LWKOPT, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
  intrinsic               :: MAX, MIN
!
!     .. Parameters ..
! real(RK), parameter      :: ZERO = 0.0_RK
!     ..
! interface
!     .. External Subroutines ..
!   include 'dlarfb.h'
!   include 'dlarft.h'
!   include 'dorgl2.h'
!   !include 'xerbla.h'
!     .. External Functions ..
!   include 'ilaenv.h'
! end interface
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  NB = mobbrmsd_ILAENV(1, 'mobbrmsd_DORGLQ', ' ', M, N, K, -1)
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
!   CALL XERBLA( 'mobbrmsd_DORGLQ', -INFO )
    return
  else if (LQUERY) then
    return
  end if
!
!     Quick return if possible
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
!        Determine when to cross over from blocked to unblocked code.
!
    NX = MAX(0, mobbrmsd_ILAENV(3, 'mobbrmsd_DORGLQ', ' ', M, N, K, -1))
    if (NX < K) then
!
!           Determine if workspace is large enough for blocked code.
!
      LDWORK = M
      IWS = LDWORK * NB
      if (LWORK < IWS) then
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
        NB = LWORK / LDWORK
        NBMIN = MAX(2, mobbrmsd_ILAENV(2, 'mobbrmsd_DORGLQ', ' ', M, N, K, -1))
      end if
    end if
  end if
!
  if (NB >= NBMIN .and. NB < K .and. NX < K) then
!
!        Use blocked code after the last block.
!        The first kk rows are handled by the block method.
!
    KI = ((K - NX - 1) / NB) * NB
    KK = MIN(K, KI + NB)
!
!        Set A(kk+1:m,1:kk) to zero.
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
!     Use unblocked code for the last or only block.
!
  if (KK < M) &
 &   call mobbrmsd_DORGL2(M - KK, N - KK, K - KK, A(KK + 1, KK + 1), LDA, &
 &                TAU(KK + 1), WORK, IINFO)
!
  if (KK > 0) then
!
!        Use blocked code
!
    do I = KI + 1, 1, -NB
      IB = MIN(NB, K - I + 1)
      if (I + IB <= M) then
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
        call mobbrmsd_DLARFT('Forward', 'Rowwise', N - I + 1, IB, A(I, I), &
       &             LDA, TAU(I), WORK, LDWORK)
!
!              Apply H**T to A(i+ib:m,i:n) from the right
!
        call mobbrmsd_DLARFB('Right', 'Transpose', 'Forward', 'Rowwise', &
       &             M - I - IB + 1, N - I + 1, IB, A(I, I), LDA, WORK, &
       &             LDWORK, A(I + IB, I), LDA, WORK(IB + 1), &
       &             LDWORK)
      end if
!
!           Apply H**T to columns i:n of current block
!
      call mobbrmsd_DORGL2(IB, N - I + 1, IB, A(I, I), LDA, TAU(I), WORK, IINFO)
!
!           Set columns 1:i-1 of current block to zero
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
!     End of mobbrmsd_DORGLQ
!
end subroutine mobbrmsd_DORGLQ

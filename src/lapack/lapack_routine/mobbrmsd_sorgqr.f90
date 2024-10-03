!> \brief \b mobbrmsd_SORGQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_SORGQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
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
!> mobbrmsd_SORGQR generates an M-by-N real matrix Q with orthonormal columns,
!> which is defined as the first N columns of a product of K elementary
!> reflectors of order M
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by mobbrmsd_SGEQRF.
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
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by mobbrmsd_SGEQRF in the first k columns of its array
!>          argument A.
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
!>          TAU is REAL array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by mobbrmsd_SGEQRF.
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
!>          The dimension of the array WORK. LWORK >= max(1,N).
!>          For optimum performance LWORK >= N*NB, where NB is the
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
!> \date December 2016
!
!> \ingroup realOTHERcomputational
!
!  =====================================================================
pure subroutine mobbrmsd_SORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
!     .. Scalar Arguments ..
  integer, intent(in)  :: K, LDA, LWORK, M, N
  integer, intent(out) :: INFO
!..
!..Array Arguments..
  real(RK), intent(inout)  :: A(LDA, *), TAU(*)
  real(RK), intent(out)    :: WORK(*)
!..
!
!  =====================================================================
!..
!..Local Scalars..
  logical :: LQUERY
  integer :: I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX
!..
!..intrinsic Functions..
  intrinsic :: MAX, MIN
!
!..Parameters..
! real(RK), parameter :: ZERO = 0.0E+0
!..
! interface
!..external Functions..
!   include 'ilaenv.h'
!..External Subroutines ..
!   include 'slarfb.h'
!   include 'slarft.h'
!   include 'sorg2r.h'
! end interface
!..
!..Executable Statements..
!
!Test the input arguments
!
  INFO = 0
  NB = mobbrmsd_ILAENV(1, 'mobbrmsd_SORGQR', ' ', M, N, K, -1)
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
!   call XERBLA('mobbrmsd_SORGQR', -INFO)
    return
  else if (LQUERY) then
    return
  end if
!
!Quick return if possible
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
    !Determine when to cross over from blocked to unblocked code.
    !
    NX = MAX(0, mobbrmsd_ILAENV(3, 'mobbrmsd_SORGQR', ' ', M, N, K, -1))
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
        NBMIN = MAX(2, mobbrmsd_ILAENV(2, 'mobbrmsd_SORGQR', ' ', M, N, K, -1))
      end if
    end if
  end if
!
  if (NB >= NBMIN .and. NB < K .and. NX < K) then
    !
    !use blocked code after the last block.
    !The first kk columns are handled by the block method.
    !
    KI = ((K - NX - 1) / NB) * NB
    KK = MIN(K, KI + NB)
    !
    !Set A(1:kk, kk + 1:n) to zero.
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
  !use unblocked code for the last or only block.
  !
  if (KK < N) call mobbrmsd_SORG2R(M - KK, N - KK, K - KK, A(KK + 1, KK + 1), LDA, TAU(KK + 1), WORK, IINFO)
  !
  if (KK > 0) then
    !
    !use blocked code
    !
    do I = KI + 1, 1, -NB
      IB = MIN(NB, K - I + 1)
      if (I + IB <= N) then
        !
        !Form the triangular factor of the block reflector
        !H = H(i) H(i + 1) ...H(i + ib - 1)
        !
        call mobbrmsd_SLARFT('Forward', 'Columnwise', M - I + 1, IB, &
            &       A(I, I), LDA, TAU(I), WORK, LDWORK)
        !
        !Apply H to A(i:m, i + ib:n) from the left
        !
        call mobbrmsd_SLARFB('Left', 'No transpose', 'Forward', &
            &       'Columnwise', M - I + 1, N - I - IB + 1, IB, &
            &       A(I, I), LDA, WORK, LDWORK, A(I, I + IB), &
            &       LDA, WORK(IB + 1), LDWORK)
      end if
      !
      !Apply H to rows i:m of current block
      !
      call mobbrmsd_SORG2R(M - I + 1, IB, IB, A(I, I), LDA, TAU(I), WORK, IINFO)
      !
      !Set rows 1:i - 1 of current block to zero
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
  !end of mobbrmsd_SORGQR
  !
end

!> \brief \b DORMLQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORMLQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dormlq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dormlq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dormlq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       pure subroutine DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       integer            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       real(RK)           ::   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORMLQ overwrites the general real M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'T':      Q**T * C       C * Q**T
!>
!> where Q is a real orthogonal matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(k) . . . H(2) H(1)
!>
!> as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left;
!>          = 'R': apply Q or Q**T from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'T':  Transpose, apply Q**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is real(RK)           :: array, dimension
!>                               (LDA,M) if SIDE = 'L',
!>                               (LDA,N) if SIDE = 'R'
!>          The i-th row must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          DGELQF in the first k rows of its array argument A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer
!>          The leading dimension of the array A. LDA >= max(1,K).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is real(RK)           :: array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGELQF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is real(RK)           :: array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is integer
!>          The leading dimension of the array C. LDC >= max(1,M).
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
!>          LWORK is integer
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
pure subroutine DORMLQ(SIDE, TRANS, M, N, K, A, LDA, TAU, C, &
               &       LDC, WORK, LWORK, INFO)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character, intent(in) :: SIDE, TRANS
  integer, intent(in)   :: K, LDA, LDC, LWORK, M, N
  integer, intent(out)  :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)    :: TAU(*)
  real(RK), intent(inout) :: A(LDA, *), C(LDC, *)
  real(RK), intent(out)   :: WORK(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  integer, parameter :: NBMAX = 64
  integer, parameter :: LDT = NBMAX + 1
  integer, parameter :: TSIZE = LDT * NBMAX
!     ..
!     .. Local Scalars ..
  logical            :: LEFT, LQUERY, NOTRAN
  character          :: TRANST
  integer            :: I, I1, I2, I3, IB, IC, IINFO, IWT, JC, &
 &                      LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Intrinsic Functions ..
  intrinsic          :: MAX, MIN
!     ..
! interface
!     .. External Subroutines ..
!   include 'dlarfb.h'
!   include 'dlarft.h'
!   include 'dorml2.h'
!   !include 'xerbla.h'
!     .. External Functions ..
!   include 'lsame.h'
!   include 'ilaenv.h'
! end interface
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  LEFT = LSAME(SIDE, 'L')
  NOTRAN = LSAME(TRANS, 'N')
  LQUERY = (LWORK == -1)
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
  if (LEFT) then
    NQ = M
    NW = MAX(1, N)
  else
    NQ = N
    NW = MAX(1, M)
  end if
  if (.not. LEFT .and. .not. LSAME(SIDE, 'R')) then
    INFO = -1
  else if (.not. NOTRAN .and. .not. LSAME(TRANS, 'T')) then
    INFO = -2
  else if (M < 0) then
    INFO = -3
  else if (N < 0) then
    INFO = -4
  else if (K < 0 .or. K > NQ) then
    INFO = -5
  else if (LDA < MAX(1, K)) then
    INFO = -7
  else if (LDC < MAX(1, M)) then
    INFO = -10
  else if (LWORK < NW .and. .not. LQUERY) then
    INFO = -12
  end if
!
  if (INFO == 0) then
!
!        Compute the workspace requirements
!
    NB = MIN(NBMAX, ILAENV(1, 'DORMLQ', SIDE//TRANS, M, N, K, -1))
    LWKOPT = NW * NB + TSIZE
    WORK(1) = LWKOPT
  end if
!
  if (INFO /= 0) then
    !CALL XERBLA( 'DORMLQ', -INFO )
    return
  else if (LQUERY) then
    return
  end if
!
!     Quick return if possible
!
  if (M == 0 .or. N == 0 .or. K == 0) then
    WORK(1) = 1
    return
  end if
!
  NBMIN = 2
  LDWORK = NW
  if (NB > 1 .and. NB < K) then
    if (LWORK < LWKOPT) then
      NB = (LWORK - TSIZE) / LDWORK
      NBMIN = MAX(2, ILAENV(2, 'DORMLQ', SIDE//TRANS, M, N, K, -1))
    end if
  end if
!
  if (NB < NBMIN .or. NB >= K) then
!
!        Use unblocked code
!
    call DORML2(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO)
  else
!
!        Use blocked code
!
    IWT = 1 + NW * NB
    if ((LEFT .and. NOTRAN) .or. (.not. LEFT .and. .not. NOTRAN)) then
      I1 = 1
      I2 = K
      I3 = NB
    else
      I1 = ((K - 1) / NB) * NB + 1
      I2 = 1
      I3 = -NB
    end if
!
    if (LEFT) then
      NI = N
      JC = 1
    else
      MI = M
      IC = 1
    end if
!
    if (NOTRAN) then
      TRANST = 'T'
    else
      TRANST = 'N'
    end if
!
    do I = I1, I2, I3
      IB = MIN(NB, K - I + 1)
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
      call DLARFT('Forward', 'Rowwise', NQ - I + 1, IB, A(I, I), &
     &             LDA, TAU(I), WORK(IWT), LDT)
      if (LEFT) then
!
!              H or H**T is applied to C(i:m,1:n)
!
        MI = M - I + 1
        IC = I
      else
!
!              H or H**T is applied to C(1:m,i:n)
!
        NI = N - I + 1
        JC = I
      end if
!
!           Apply H or H**T
!
      call DLARFB(SIDE, TRANST, 'Forward', 'Rowwise', MI, NI, IB, &
     &             A(I, I), LDA, WORK(IWT), LDT, &
     &             C(IC, JC), LDC, WORK, LDWORK)
    end do
  end if
  WORK(1) = LWKOPT
  return
!
!     End of DORMLQ
!
end subroutine DORMLQ


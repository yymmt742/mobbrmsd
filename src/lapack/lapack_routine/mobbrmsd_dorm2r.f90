!> \brief \b mobbrmsd_DORM2R multiplies a general matrix by the orthogonal matrix from a QR factorization determined by sgeqrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DORM2R + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dorm2r.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dorm2r.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dorm2r.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, M, N
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
!> mobbrmsd_DORM2R overwrites the general real m by n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**T* C  if SIDE = 'L' and TRANS = 'T', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**T if SIDE = 'R' and TRANS = 'T',
!>
!> where Q is a real orthogonal matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by mobbrmsd_DGEQRF. Q is of order m if SIDE = 'L' and of order n
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left
!>          = 'R': apply Q or Q**T from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'T': apply Q**T (Transpose)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is real(RK)           :: array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          mobbrmsd_DGEQRF in the first k columns of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is real(RK)           :: array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by mobbrmsd_DGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is real(RK)           :: array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is real(RK)           :: array, dimension
!>                                   (N) if SIDE = 'L',
!>                                   (M) if SIDE = 'R'
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
pure subroutine mobbrmsd_DORM2R(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: SIDE, TRANS
  integer, intent(in)      :: K, LDA, LDC, M, N
  integer, intent(out)     :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)     :: TAU(*)
  real(RK), intent(inout)  :: A(LDA, *), C(LDC, *)
  real(RK), intent(out)    :: WORK(*)
!     ..
!
!  =====================================================================
!     .. Local Scalars ..
  logical                 :: LEFT, NOTRAN
  integer                 :: I, I1, I2, I3, IC, JC, MI, NI, NQ
  real(RK)                :: AII
!     ..
!     .. Intrinsic Functions ..
  intrinsic               :: MAX
!     .. Parameters ..
! real(RK), parameter      :: ONE = 1.0_RK
!     ..
! interface
!     .. External Subroutines ..
!   include 'dlarf.h'
!   !include 'xerbla.h'
!     .. External Functions ..
!   include 'lsame.h'
! end interface
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  LEFT = mobbrmsd_LSAME(SIDE, 'L')
  NOTRAN = mobbrmsd_LSAME(TRANS, 'N')
!
!     NQ is the order of Q
!
  if (LEFT) then
    NQ = M
  else
    NQ = N
  end if
  if (.not. LEFT .and. .not. mobbrmsd_LSAME(SIDE, 'R')) then
    INFO = -1
  else if (.not. NOTRAN .and. .not. mobbrmsd_LSAME(TRANS, 'T')) then
    INFO = -2
  else if (M < 0) then
    INFO = -3
  else if (N < 0) then
    INFO = -4
  else if (K < 0 .or. K > NQ) then
    INFO = -5
  else if (LDA < MAX(1, NQ)) then
    INFO = -7
  else if (LDC < MAX(1, M)) then
    INFO = -10
  end if
  if (INFO /= 0) then
    !CALL XERBLA( 'mobbrmsd_DORM2R', -INFO )
    return
  end if
!
!     Quick return if possible
!
  if (M == 0 .or. N == 0 .or. K == 0) return
!
  if ((LEFT .and. .not. NOTRAN) .or. (.not. LEFT .and. NOTRAN)) then
    I1 = 1
    I2 = K
    I3 = 1
  else
    I1 = K
    I2 = 1
    I3 = -1
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
  do I = I1, I2, I3
    if (LEFT) then
!
!           H(i) is applied to C(i:m,1:n)
!
      MI = M - I + 1
      IC = I
    else
!
!           H(i) is applied to C(1:m,i:n)
!
      NI = N - I + 1
      JC = I
    end if
!
!        Apply H(i)
!
    AII = A(I, I)
    A(I, I) = ONE
    call mobbrmsd_DLARF(SIDE, MI, NI, A(I, I), 1, TAU(I), C(IC, JC), &
   &            LDC, WORK)
    A(I, I) = AII
  end do
  return
!
!     End of mobbrmsd_DORM2R
!
end subroutine mobbrmsd_DORM2R


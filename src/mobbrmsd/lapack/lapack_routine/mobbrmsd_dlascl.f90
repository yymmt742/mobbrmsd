!> \brief \b mobbrmsd_DLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DLASCL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlascl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlascl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlascl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       INTEGER            INFO, KL, KU, LDA, M, N
!       DOUBLE PRECISION   CFROM, CTO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DLASCL multiplies the M by N real matrix A by the real scalar
!> CTO/CFROM.  This is done without over/underflow as long as the final
!> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!> A may be full, upper triangular, lower triangular, upper Hessenberg,
!> or banded.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          TYPE indices the storage type of the input matrix.
!>          = 'G':  A is a full matrix.
!>          = 'L':  A is a lower triangular matrix.
!>          = 'U':  A is an upper triangular matrix.
!>          = 'H':  A is an upper Hessenberg matrix.
!>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the lower
!>                  half stored.
!>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the upper
!>                  half stored.
!>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!>                  bandwidth KU. See DGBTRF for storage details.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] CFROM
!> \verbatim
!>          CFROM is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] CTO
!> \verbatim
!>          CTO is DOUBLE PRECISION
!>
!>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!>          without over/underflow if the final result CTO*A(I,J)/CFROM
!>          can be represented without over/underflow.  CFROM must be
!>          nonzero.
!> \endverbatim
!>
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!>          storage type.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If TYPE = 'G', 'L', 'U', 'H', LDA >= max(1,M);
!>             TYPE = 'B', LDA >= KL+1;
!>             TYPE = 'Q', LDA >= KU+1;
!>             TYPE = 'Z', LDA >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  - successful exit
!>          <0 - if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
pure subroutine mobbrmsd_DLASCL(type, KL, KU, CFROM, CTO, M, N, A, LDA, INFO)
! use LA_CONSTANTS, only: RK => dp
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: type
  integer, intent(in)      :: KL, KU, LDA, M, N
  integer, intent(out)     :: INFO
  real(RK), intent(in)     :: CFROM, CTO
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout)  :: A(LDA, *)
!     ..
!     .. Local Scalars ..
  logical                 :: DONE
  integer                 :: I, ITYPE, J, K1, K2, K3, K4
  real(RK)                :: BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. Intrinsic Functions ..
  intrinsic               :: ABS, MAX, MIN
!     .. Parameters ..
! real(RK), parameter      :: ZERO = 0.0_RK
! real(RK), parameter      :: ONE = 1.0_RK
!     ..
! interface
!   include 'lsame.h'
!   include 'disnan.h'
!   include 'dlamch.h'
!   !include 'xerbla.h'
! end interface
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
!
  if (mobbrmsd_LSAME(type, 'G')) then
    ITYPE = 0
  else if (mobbrmsd_LSAME(type, 'L')) then
    ITYPE = 1
  else if (mobbrmsd_LSAME(type, 'U')) then
    ITYPE = 2
  else if (mobbrmsd_LSAME(type, 'H')) then
    ITYPE = 3
  else if (mobbrmsd_LSAME(type, 'B')) then
    ITYPE = 4
  else if (mobbrmsd_LSAME(type, 'Q')) then
    ITYPE = 5
  else if (mobbrmsd_LSAME(type, 'Z')) then
    ITYPE = 6
  else
    ITYPE = -1
  end if
!
  if (ITYPE == -1) then
    INFO = -1
  else if (CFROM == ZERO .or. mobbrmsd_DISNAN(CFROM)) then
    INFO = -4
  else if (mobbrmsd_DISNAN(CTO)) then
    INFO = -5
  else if (M < 0) then
    INFO = -6
  else if (N < 0 .or. (ITYPE == 4 .and. N /= M) .or.&
 &         (ITYPE == 5 .and. N /= M)) then
    INFO = -7
  else if (ITYPE <= 3 .and. LDA < MAX(1, M)) then
    INFO = -9
  else if (ITYPE >= 4) then
    if (KL < 0 .or. KL > MAX(M - 1, 0)) then
      INFO = -2
    else if (KU < 0 .or. KU > MAX(N - 1, 0) .or.               &
&            ((ITYPE == 4 .or. ITYPE == 5) .and. KL /= KU)) &
&             then
      INFO = -3
    else if ((ITYPE == 4 .and. LDA < KL + 1) .or. &
&            (ITYPE == 5 .and. LDA < KU + 1) .or. &
&            (ITYPE == 6 .and. LDA < 2 * KL + KU + 1)) then
      INFO = -9
    end if
  end if
!
  if (INFO /= 0) then
    !CALL XERBLA( 'DLASCL', -INFO )
    return
  end if
!
!     Quick return if possible
!
  if (N == 0 .or. M == 0) return
!
!     Get machine parameters
!
  SMLNUM = mobbrmsd_DLAMCH('S')
  BIGNUM = ONE / SMLNUM
  DONE = .false.
!
  CFROMC = CFROM
  CTOC = CTO
!
10 continue
  CFROM1 = CFROMC!SMLNUM
  if (CFROM1 == CFROMC) then
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
    MUL = CTOC / CFROMC
    DONE = .true.
    CTO1 = CTOC
  else
    CTO1 = CTOC / BIGNUM
    if (CTO1 == CTOC) then
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
      MUL = CTOC
      DONE = .true.
      CFROMC = ONE
    else if (ABS(CFROM1) > ABS(CTOC) .and. CTOC /= ZERO) then
      MUL = SMLNUM
      DONE = .false.
      CFROMC = CFROM1
    else if (ABS(CTO1) > ABS(CFROMC)) then
      MUL = BIGNUM
      DONE = .false.
      CTOC = CTO1
    else
      MUL = CTOC / CFROMC
      DONE = .true.
    end if
  end if
!
  if (ITYPE == 0) then
!
!        Full matrix
!
    do concurrent(J=1:N, I=1:M)
      A(I, J) = A(I, J) * MUL
    end do
!         do 30 J = 1, N
!           do 20 I = 1, M
!             A(I, J) = A(I, J) * MUL
!20         continue
!30       continue
!
  else if (ITYPE == 1) then
!
!        Lower triangular matrix
!
    do J = 1, N
      do I = J, M
        A(I, J) = A(I, J) * MUL
      end do
    end do
!
  else if (ITYPE == 2) then
!
!        Upper triangular matrix
!
    do J = 1, N
      do I = 1, MIN(J, M)
        A(I, J) = A(I, J) * MUL
      end do
    end do
!
  else if (ITYPE == 3) then
!
!        Upper Hessenberg matrix
!
    do J = 1, N
      do I = 1, MIN(J + 1, M)
        A(I, J) = A(I, J) * MUL
      end do
    end do
!
  else if (ITYPE == 4) then
!
!        Lower half of a symmetric band matrix
!
    K3 = KL + 1
    K4 = N + 1
    do J = 1, N
      do I = 1, MIN(K3, K4 - J)
        A(I, J) = A(I, J) * MUL
      end do
    end do
!
  else if (ITYPE == 5) then
!
!        Upper half of a symmetric band matrix
!
    K1 = KU + 2
    K3 = KU + 1
    do J = 1, N
      do I = MAX(K1 - J, 1), K3
        A(I, J) = A(I, J) * MUL
      end do
    end do
!
  else if (ITYPE == 6) then
!
!        Band matrix
!
    K1 = KL + KU + 2
    K2 = KL + 1
    K3 = 2 * KL + KU + 1
    K4 = KL + KU + 1 + M
    do J = 1, N
      do I = MAX(K1 - J, K2), MIN(K3, K4 - J)
        A(I, J) = A(I, J) * MUL
      end do
    end do
!
  end if
!
  if (.not. DONE) GO TO 10
!
  return
!
! End of mobbrmsd_DLASCL
!
end subroutine mobbrmsd_DLASCL

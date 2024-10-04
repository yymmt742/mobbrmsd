!> \brief \b mobbrmsd_DORGL2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DORGL2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dorgl2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dorgl2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dorgl2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
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
!> mobbrmsd_DORGL2 generates an m by n real matrix Q with orthonormal rows,
!> which is defined as the first m rows of a product of k elementary
!> reflectors of order n
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
!>          On exit, the m-by-n matrix Q.
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
!>          WORK is real(RK)           :: array, dimension (M)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
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
pure subroutine mobbrmsd_DORGL2(M, N, K, A, LDA, TAU, WORK, INFO)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)     :: K, LDA, M, N
  integer, intent(out)    :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)    :: TAU(*)
  real(RK), intent(inout) :: A(LDA, *)
  real(RK), intent(out)   :: WORK(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  integer                 :: I, J, L
!     ..
!     .. Intrinsic Functions ..
  intrinsic               :: MAX
!     ..
!     .. Parameters ..
! real(RK), parameter      :: ZERO = 0.0_RK
! real(RK), parameter      :: ONE = 1.0_RK
!     ..
! interface
!     .. External Subroutines ..
!   include 'dlarf.h'
!   include 'dscal.h'
!   !include 'xerbla.h'
! end interface
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  if (M < 0) then
    INFO = -1
  else if (N < M) then
    INFO = -2
  else if (K < 0 .or. K > M) then
    INFO = -3
  else if (LDA < MAX(1, M)) then
    INFO = -5
  end if
  if (INFO /= 0) then
    !CALL XERBLA( 'DORGL2', -INFO )
    return
  end if
!
!     Quick return if possible
!
  if (M <= 0) return
!
  if (K < M) then
!
!        Initialise rows k+1:m to rows of the unit matrix
!
    do J = 1, N
      do L = K + 1, M
        A(L, J) = ZERO
      end do
      if (J > K .and. J <= M) A(J, J) = ONE
    end do
  end if
!
  do I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the right
!
    if (I < N) then
      if (I < M) then
        A(I, I) = ONE
        call mobbrmsd_DLARF('Right', M - I, N - I + 1, A(I, I), LDA, &
       &            TAU(I), A(I + 1, I), LDA, WORK)
      end if
      call mobbrmsd_DSCAL(N - I, -TAU(I), A(I, I + 1), LDA)
    end if
    A(I, I) = ONE - TAU(I)
!
!        Set A(i,1:i-1) to zero
!
    do L = 1, I - 1
      A(I, L) = ZERO
    end do
  end do
  return
!
!     End of mobbrmsd_DORGL2
!
end subroutine mobbrmsd_DORGL2


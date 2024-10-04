!> \brief \b mobbrmsd_DORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by sgeqrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DORG2R + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dorg2r.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dorg2r.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dorg2r.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
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
!> mobbrmsd_DORG2R generates an m by n real matrix Q with orthonormal columns,
!> which is defined as the first n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by mobbrmsd_DGEQRF.
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
!>          A is real(RK)           :: array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by mobbrmsd_DGEQRF in the first k columns of its array
!>          argument A.
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
!>          reflector H(i), as returned by mobbrmsd_DGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is real(RK)           :: array, dimension (N)
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
pure subroutine mobbrmsd_DORG2R(M, N, K, A, LDA, TAU, WORK, INFO)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)  :: K, LDA, M, N
  integer, intent(out) :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)    :: TAU(*)
  real(RK), intent(inout) :: A(LDA, *)
  real(RK), intent(out)   :: WORK(*)
!     ..
!
!  =====================================================================
!     .. Local Scalars ..
  integer                 :: I, J, L
!     .. Parameters ..
! real(RK), parameter      ::   ZERO = 0.0_RK
! real(RK), parameter      ::   ONE = 1.0_RK
!     ..
! interface
!     .. External Subroutines ..
!   include 'dlarf.h'
!   include 'dscal.h'
!   !include 'xerbla.h'
! end interface
!     ..
!     .. Intrinsic Functions ..
  intrinsic               :: MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  if (M < 0) then
    INFO = -1
  else if (N < 0 .or. N > M) then
    INFO = -2
  else if (K < 0 .or. K > N) then
    INFO = -3
  else if (LDA < MAX(1, M)) then
    INFO = -5
  end if
  if (INFO /= 0) then
    !CALL XERBLA( 'DORG2R', -INFO )
    return
  end if
!
!     Quick return if possible
!
  if (N <= 0) return
!
!     Initialise columns k+1:n to columns of the unit matrix
!
  do concurrent(J=K + 1:N)
    do concurrent(L=1:M)
      A(L, J) = ZERO
    end do
    A(J, J) = ONE
  end do
!       do 20 J = K + 1, N
!         do 10 L = 1, M
!           A(L, J) = ZERO
!10       continue
!         A(J, J) = ONE
!20     continue
!
  do I = K, 1, -1
!
!         Apply H(i) to A(i:m,i:n) from the left
!
    if (I < N) then
      A(I, I) = ONE
      call mobbrmsd_DLARF('Left', M - I + 1, N - I, A(I, I), 1, TAU(I), &
     &            A(I, I + 1), LDA, WORK)
    end if
    if (I < M) call mobbrmsd_DSCAL(M - I, -TAU(I), A(I + 1, I), 1)
    A(I, I) = ONE - TAU(I)
!
!    Set A(1:i-1,i) to zero
!
    do concurrent(L=1:I - 1)
      A(L, I) = ZERO
    end do
!         do 30 L = 1, I - 1
!           A(L, I) = ZERO
!30       continue
  end do
  return
!
!     End of mobbrmsd_DORG2R
!
end subroutine mobbrmsd_DORG2R


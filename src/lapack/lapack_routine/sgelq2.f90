!> \brief \b SGELQ2 computes the LQ factorization of a general rectangular matrix using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGELQ2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelq2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelq2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelq2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGELQ2( M, N, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
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
!> SGELQ2 computes an LQ factorization of a real m-by-n matrix A:
!>
!>    A = ( L 0 ) *  Q
!>
!> where:
!>
!>    Q is a n-by-n orthogonal matrix;
!>    L is a lower-triangular m-by-m matrix;
!>    0 is a m-by-(n-m) zero matrix, if m < n.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the m by n matrix A.
!>          On exit, the elements on and below the diagonal of the array
!>          contain the m by min(m,n) lower trapezoidal matrix L (L is
!>          lower triangular if m <= n); the elements above the diagonal,
!>          with the array TAU, represent the orthogonal matrix Q as a
!>          product of elementary reflectors (see Further Details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (M)
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
!> \date November 2019
!
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(k) . . . H(2) H(1), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
!>  and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
pure subroutine SGELQ2(M, N, A, LDA, TAU, WORK, INFO)
  implicit none
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
  integer, intent(in)  :: LDA, M, N
  integer, intent(out) :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: A(LDA, *)
  real(RK), intent(out)   :: TAU(*), WORK(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
  integer :: I, K
  real(RK) :: AII
!
!     .. Parameters ..
! real(RK), parameter :: ONE = 1.0E+0
!     ..
! interface
! .. External Subroutines ..
!   include 'slarf.h'
!   include 'slarfg.h'
!   include 'xerbla.h'
! end interface
!..
!..intrinsic Functions..
  intrinsic MAX, MIN
!..
!..Executable Statements..
!
! Test the input arguments
!
  INFO = 0
  if (M < 0) then
    INFO = -1
  else if (N < 0) then
    INFO = -2
  else if (LDA < MAX(1, M)) then
    INFO = -4
  end if
  if (INFO /= 0) then
! call XERBLA('SGELQ2', -INFO)
    return
  end if
!
  K = MIN(M, N)
!
  do I = 1, K
!
! Generate elementary reflector H(i) to annihilate A(i, i + 1:n)
!
    call SLARFG(N - I + 1, A(I, I), A(I, MIN(I + 1, N)), LDA, TAU(I))
    if (I < M) then
!
!Apply H(i) to A(i + 1:m, i:n) from the right
!
      AII = A(I, I)
      A(I, I) = ONE
      call SLARF('Right', M - I, N - I + 1, A(I, I), LDA, TAU(I), A(I + 1, I), LDA, WORK)
      A(I, I) = AII
    end if
  end do
  return
!
! end of SGELQ2
!
end

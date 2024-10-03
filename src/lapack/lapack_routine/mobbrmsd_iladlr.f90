!> \brief \b mobbrmsd_ILADLR scans a matrix for its last non-zero row.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_ILADLR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION mobbrmsd_ILADLR( M, N, A, LDA )
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
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
!> mobbrmsd_ILADLR scans A for its last non-zero row.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
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
pure function mobbrmsd_ILADLR(M, N, A, LDA)
! use LA_CONSTANTS, only: DP
  implicit none
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)  :: M, N, LDA
!     ..
!     .. Array Arguments ..
  real(RK), intent(in) ::  A(LDA, *)
!     ..
  integer :: mobbrmsd_ILADLR
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
  integer :: I, J
!
!     .. Parameters ..
! real(RK), parameter  :: ZERO = 0.0_DP
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
  if (M == 0) then
    mobbrmsd_ILADLR = M
  else if (A(M, 1) /= ZERO .or. A(M, N) /= ZERO) then
    mobbrmsd_ILADLR = M
  else
!     Scan up each column tracking the last zero row seen.
    mobbrmsd_ILADLR = 0
    do J = 1, N
      I = M
      do while ((A(MAX(I, 1), J) == ZERO) .and. (I >= 1))
        I = I - 1
      end do
      mobbrmsd_ILADLR = MAX(mobbrmsd_ILADLR, I)
    end do
  end if
  return
end function mobbrmsd_ILADLR

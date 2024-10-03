!> \brief \b ILADLC scans a matrix for its last non-zero column.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
! http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILADLC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/iladlc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/iladlc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/iladlc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!   INTEGER FUNCTION ILADLC( M, N, A, LDA )
!
!   .. Scalar Arguments ..
!   INTEGER            M, N, LDA
!   ..
!   .. Array Arguments ..
!   real(DP)           ::   A( LDA, * )
!   ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILADLC scans A for its last non-zero column.
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
!>          A is real(DP)           :: array, dimension (LDA,N)
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
pure function ILADLC(M, N, A, LDA)
  use LA_CONSTANTS, only: DP
  implicit none
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
! .. Scalar Arguments ..
  integer, intent(in)  :: M, N, LDA
! ..
! .. Array Arguments ..
  real(DP), intent(in) :: A(LDA, *)
  integer              :: ILADLC
! ..
!
!  =====================================================================
!
! .. Parameters ..
  real(DP), parameter  :: ZERO = 0.0_DP
! ..
! .. Local Scalars ..
  integer             :: I
! ..
! .. Executable Statements ..
!
! Quick test for the common case where one corner is non-zero.
  if (N == 0) then
    ILADLC = N
  else if (A(1, N) /= ZERO .or. A(M, N) /= ZERO) then
    ILADLC = N
  else
! Now scan each column from the end, returning with the first non-zero.
    do ILADLC = N, 1, -1
      do I = 1, M
        if (A(I, ILADLC) /= ZERO) return
      end do
    end do
  end if
  return
end function ILADLC

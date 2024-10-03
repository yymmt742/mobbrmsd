!> \brief \b mobbrmsd_ILADLC scans a matrix for its last non-zero column.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
! http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_ILADLC + dependencies
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
!   INTEGER FUNCTION mobbrmsd_ILADLC( M, N, A, LDA )
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
!> mobbrmsd_ILADLC scans A for its last non-zero column.
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
pure function mobbrmsd_ILADLC(M, N, A, LDA)
! use LA_CONSTANTS, only: DP
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
  real(RK), intent(in) :: A(LDA, *)
  integer :: mobbrmsd_ILADLC
! ..
!
!  =====================================================================
! ..
! .. Local Scalars ..
  integer :: I
!
! .. Parameters ..
! real(RK), parameter  :: ZERO = 0.0_DP
! ..
! .. Executable Statements ..
!
! Quick test for the common case where one corner is non-zero.
  if (N == 0) then
    mobbrmsd_ILADLC = N
  else if (A(1, N) /= ZERO .or. A(M, N) /= ZERO) then
    mobbrmsd_ILADLC = N
  else
! Now scan each column from the end, returning with the first non-zero.
    do mobbrmsd_ILADLC = N, 1, -1
      do I = 1, M
        if (A(I, mobbrmsd_ILADLC) /= ZERO) return
      end do
    end do
  end if
  return
end function mobbrmsd_ILADLC

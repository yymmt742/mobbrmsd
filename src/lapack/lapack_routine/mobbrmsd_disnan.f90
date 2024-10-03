!> \brief \b mobbrmsd_DISNAN tests input for NaN.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION mobbrmsd_DISNAN( DIN )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION, INTENT(IN) :: DIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN
!> \verbatim
!>          DIN is DOUBLE PRECISION
!>          Input to test for NaN.
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
pure elemental function mobbrmsd_DISNAN(DIN)
! use LA_CONSTANTS, only: RK => dp
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  real(RK), intent(in) :: DIN
  logical              :: mobbrmsd_DISNAN
!     ..
!
!  =====================================================================
!
!  .. External Functions ..
! interface
!   include 'dlaisnan.h'
! end interface
!  ..
!  .. Executable Statements ..
  mobbrmsd_DISNAN = mobbrmsd_Dmobbrmsd_LAISNAN(DIN, DIN)
  return
end function mobbrmsd_DISNAN

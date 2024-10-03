!> \brief \b mobbrmsd_SISNAN tests input for NaN.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_SISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sisnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sisnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sisnan.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION mobbrmsd_SISNAN( SIN )
!
!       .. Scalar Arguments ..
!       REAL, INTENT(IN) :: SIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_SISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIN
!> \verbatim
!>          SIN is REAL
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
!> \date June 2017
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
pure elemental function mobbrmsd_SISNAN(SI)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
  real(RK), intent(in) :: SI
  logical              :: mobbrmsd_SISNAN
!..
!
!  =====================================================================
!
! interface
! .. External Functions ..
!   include 'slaisnan.h'
! end interface
!..
!..Executable Statements..
  mobbrmsd_SISNAN = Smobbrmsd_LAISNAN(SI, SI)
  return
end

!> \brief \b SLAPY2 returns sqrt(x2+y2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAPY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLAPY2( X, Y )
!
!       .. Scalar Arguments ..
!       REAL               X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!> overflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is REAL
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is REAL
!>          X and Y specify the values x and y.
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
pure elemental function SLAPY2(X, Y)
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
  real(RK), intent(in) :: X, Y
  real(RK) :: SLAPY2
!..
!
!  =====================================================================
!
!..Local Scalars..
  real(RK) :: W, XABS, YABS, Z
  logical :: X_IS_NAN, Y_IS_NAN
!..
!..Parameters..
! real(RK), parameter :: ZERO = 0.0E0
! real(RK), parameter :: ONE = 1.0E0
!..
! interface
! .. External Functions ..
!   include 'sisnan.h'
! end interface
!..
!..intrinsic Functions..
  intrinsic :: ABS, MAX, MIN, SQRT
!..
!..Executable Statements..
!
  X_IS_NAN = SISNAN(X)
  Y_IS_NAN = SISNAN(Y)
  if (X_IS_NAN) SLAPY2 = X
  if (Y_IS_NAN) SLAPY2 = Y
!
  if (.not. (X_IS_NAN .or. Y_IS_NAN)) then
    XABS = ABS(X)
    YABS = ABS(Y)
    W = MAX(XABS, YABS)
    Z = MIN(XABS, YABS)
    if (Z == ZERO) then
      SLAPY2 = W
    else
      SLAPY2 = W * SQRT(ONE + (Z / W)**2)
    end if
  end if
  return
!
!end of SLAPY2
!
end

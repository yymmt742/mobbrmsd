!> \brief \b DLAPY2 returns sqrt(x2+y2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAPY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dlapy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dlapy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dlapy2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       real(RK)           :: FUNCTION DLAPY2( X, Y )
!
!       .. Scalar Arguments ..
!       real(RK)           ::   X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!> overflow and unnecessary underflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is real(RK)           ::
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is real(RK)           ::
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
pure function dlapy2(X, Y)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  real(RK), intent(in) :: X, Y
  real(RK)            :: dlapy2
!     ..
!  =====================================================================
!     ..
!     .. Local Scalars ..
  real(RK)             :: W, XABS, YABS, Z, HUGEVAL
  logical              :: X_IS_NAN, Y_IS_NAN
!
!     .. Parameters ..
! real(RK), parameter   :: ZERO = 0.0_RK
! real(RK), parameter   :: ONE = 1.0_RK
!     ..
! interface
!     .. External Functions ..
!   include 'disnan.h'
!     .. External Subroutines ..
!   include 'dlamch.h'
! end interface
!     ..
!     .. Intrinsic Functions ..
  intrinsic            :: ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
  X_IS_NAN = DISNAN(X)
  Y_IS_NAN = DISNAN(Y)
  if (X_IS_NAN) DLAPY2 = X
  if (Y_IS_NAN) DLAPY2 = Y
  HUGEVAL = DLAMCH('Overflow')
!
  if (.not. (X_IS_NAN .or. Y_IS_NAN)) then
    XABS = ABS(X)
    YABS = ABS(Y)
    W = MAX(XABS, YABS)
    Z = MIN(XABS, YABS)
    if (Z == ZERO .or. W > HUGEVAL) then
      DLAPY2 = W
    else
      DLAPY2 = W * SQRT(ONE + (Z / W)**2)
    end if
  end if
  return
!
!     End of DLAPY2
!
end function DLAPY2


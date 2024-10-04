!> \brief \b mobbrmsd_SLAMCH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
! http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!  REAL             FUNCTION mobbrmsd_SLAMCH( CMACH )
!
! .. Scalar Arguments ..
!  CHARACTER          CMACH
! ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_SLAMCH determines single precision machine parameters.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CMACH
!> \verbatim
!>          CMACH is CHARACTER*1
!>          Specifies the value to be returned by mobbrmsd_SLAMCH:
!>          = 'E' or 'e',   mobbrmsd_SLAMCH := eps
!>          = 'S' or 's ,   mobbrmsd_SLAMCH := sfmin
!>          = 'B' or 'b',   mobbrmsd_SLAMCH := base
!>          = 'P' or 'p',   mobbrmsd_SLAMCH := eps*base
!>          = 'N' or 'n',   mobbrmsd_SLAMCH := t
!>          = 'R' or 'r',   mobbrmsd_SLAMCH := rnd
!>          = 'M' or 'm',   mobbrmsd_SLAMCH := emin
!>          = 'U' or 'u',   mobbrmsd_SLAMCH := rmin
!>          = 'L' or 'l',   mobbrmsd_SLAMCH := emax
!>          = 'O' or 'o',   mobbrmsd_SLAMCH := rmax
!>          where
!>          eps   = relative machine precision
!>          sfmin = safe minimum, such that 1/sfmin does not overflow
!>          base  = base of the machine
!>          prec  = eps*base
!>          t     = number of (base) digits in the mantissa
!>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!>          emin  = minimum exponent before (gradual) underflow
!>          rmin  = underflow threshold - base**(emin-1)
!>          emax  = largest exponent before overflow
!>          rmax  = overflow threshold  - (base**emax)*(1-eps)
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
!> \date December 2016
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
pure elemental function mobbrmsd_SLAMCH(CMACH)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! December 2016
!
! .. Scalar Arguments ..
  character, intent(in) :: CMACH
  real(RK)              :: mobbrmsd_SLAMCH
! ..
!
! =====================================================================
!
! .. Local Scalars ..
  real(RK) :: RND, EPS, SFMIN, SMALL, RMACH
!
! .. Parameters ..
! real(RK), parameter :: ZERO = 0.0E+0
! real(RK), parameter :: ONE = 1.0E+0
! ..
! .. External Functions ..
! interface
!   include 'lsame.h'
! end interface
! ..
! ..
! .. Intrinsic Functions ..
  intrinsic :: DIGITS, EPSILON, HUGE, MAXEXPONENT, MINEXPONENT, RADIX, TINY
! ..
! .. Executable Statements ..
!
!
! Assume rounding, not chopping. Always.
!
  RND = ONE
!
  if (ONE == RND) then
    EPS = EPSILON(ZERO) ! 0.5
  else
    EPS = EPSILON(ZERO)
  end if
!
  if (mobbrmsd_LSAME(CMACH, 'E')) then
    RMACH = EPS
  else if (mobbrmsd_LSAME(CMACH, 'S')) then
    SFMIN = TINY(ZERO)
    SMALL = ONE / HUGE(ZERO)
    if (SMALL >= SFMIN) then
!
!   Use SMALL plus a bit, to avoid the possibility of rounding
!   causing overflow when computing  1/sfmin.
!
      SFMIN = SMALL * (ONE + EPS)
    end if
    RMACH = SFMIN
  else if (mobbrmsd_LSAME(CMACH, 'B')) then
    RMACH = RADIX(ZERO)
  else if (mobbrmsd_LSAME(CMACH, 'P')) then
    RMACH = EPS ! RADIX(ZERO)
  else if (mobbrmsd_LSAME(CMACH, 'N')) then
    RMACH = DIGITS(ZERO)
  else if (mobbrmsd_LSAME(CMACH, 'R')) then
    RMACH = RND
  else if (mobbrmsd_LSAME(CMACH, 'M')) then
    RMACH = MINEXPONENT(ZERO)
  else if (mobbrmsd_LSAME(CMACH, 'U')) then
    RMACH = TINY(zero)
  else if (mobbrmsd_LSAME(CMACH, 'L')) then
    RMACH = MAXEXPONENT(ZERO)
  else if (mobbrmsd_LSAME(CMACH, 'O')) then
    RMACH = HUGE(ZERO)
  else
    RMACH = ZERO
  end if
!
  mobbrmsd_SLAMCH = RMACH
  return
!
! End of mobbrmsd_SLAMCH
!
end
!***********************************************************************
!> \brief \b SLAMC3
!> \details
!> \b Purpose:
!> \verbatim
!> SLAMC3  is intended to force  A  and  B  to be stored prior to doing
!> the addition of  A  and  B ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \date December 2016
!> \ingroup auxOTHERauxiliary
!>
!> \param[in] A
!> \verbatim
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          The values A and B.
!> \endverbatim
!>
!
pure elemental function SLAMC3(A, B)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
! Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
! November 2010
!
! .. Scalar Arguments ..
  real(RK), intent(in) :: A, B
  real(RK)             :: SLAMC3
! ..
! =====================================================================
!
! .. Executable Statements ..
!
  SLAMC3 = A + B
!
  return
!
! End of SLAMC3
!
end
!
!***********************************************************************

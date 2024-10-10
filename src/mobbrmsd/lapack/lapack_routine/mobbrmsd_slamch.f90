!| mobbrmsd_SLAMCH determines single precision machine parameters.
!  Assume rounding, not chopping. Always.
!
!  reference SLAMCH is provided by [netlib](http://www.netlib.org/lapack/explore-html/)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  \author Univ. of Tennessee
!  \author Univ. of California Berkeley
!  \author Univ. of Colorado Denver
!  \author NAG Ltd.
!  \date December 2016
! December 2016
!
pure elemental function mobbrmsd_SLAMCH(CMACH)
  implicit none
  character, intent(in) :: CMACH
!!          Specifies the value to be returned by mobbrmsd_SLAMCH:
!!
!!          = 'E' or 'e',   mobbrmsd_SLAMCH := eps
!!
!!          = 'S' or 's ,   mobbrmsd_SLAMCH := sfmin
!!
!!          = 'B' or 'b',   mobbrmsd_SLAMCH := base
!!
!!          = 'P' or 'p',   mobbrmsd_SLAMCH := eps*base
!!
!!          = 'N' or 'n',   mobbrmsd_SLAMCH := t
!!
!!          = 'R' or 'r',   mobbrmsd_SLAMCH := rnd
!!
!!          = 'M' or 'm',   mobbrmsd_SLAMCH := emin
!!
!!          = 'U' or 'u',   mobbrmsd_SLAMCH := rmin
!!
!!          = 'L' or 'l',   mobbrmsd_SLAMCH := emax
!!
!!          = 'O' or 'o',   mobbrmsd_SLAMCH := rmax
!!
!!          where
!!
!!          eps   = relative machine precision
!!
!!          sfmin = safe minimum, such that 1/sfmin does not overflow <br>
!!
!!          base  = base of the machine
!!
!!          prec  = eps*base <br>
!!
!!          t     = number of (base) digits in the mantissa
!!
!!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!!
!!          emin  = minimum exponent before (gradual) underflow
!!
!!          rmin  = underflow threshold - base**(emin-1)
!!
!!          emax  = largest exponent before overflow
!!
!!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!!
  real(RK)  :: mobbrmsd_SLAMCH
!! machine parameter.
!!
  intrinsic :: DIGITS, EPSILON, HUGE, MAXEXPONENT, MINEXPONENT, RADIX, TINY
  real(RK)  :: RND, EPS, SFMIN, SMALL, RMACH
! real(RK), parameter :: ZERO = 0.0E+0
! real(RK), parameter :: ONE = 1.0E+0
! interface
!   include 'lsame.h'
! end interface
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
end


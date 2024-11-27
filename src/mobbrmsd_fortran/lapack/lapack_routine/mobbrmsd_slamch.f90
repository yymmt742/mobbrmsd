!| determines single precision machine parameters.
!
!  mobbrmsd_SLAMCH determines single precision machine parameters.
!  Assume rounding, not chopping. Always.
!
!  reference SLAMCH is provided by [netlib](http://www.netlib.org/lapack/explore-html/)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! December 2016
!
pure elemental function mobbrmsd_SLAMCH(CMACH)
  implicit none
  character, intent(in) :: CMACH
!!   Specifies the value to be returned by mobbrmsd_SLAMCH:
!!
!!   = 'E' or 'e',   mobbrmsd_SLAMCH := EPS;
!!   relative machine precision.
!!
!!   = 'S' or 's ,   mobbrmsd_SLAMCH := SFMIN;
!!   safe minimum, such that 1/sfmin does not overflow.
!!
!!   = 'B' or 'b',   mobbrmsd_SLAMCH := BASE;
!!   base of the machine.
!!
!!   = 'P' or 'p',   mobbrmsd_SLAMCH := EPS * BASE;
!!   PREC  = EPS * BASE.
!!
!!   = 'N' or 'n',   mobbrmsd_SLAMCH := T;
!!   number of (base) digits in the mantissa.
!!
!!   = 'R' or 'r',   mobbrmsd_SLAMCH := RND;
!!   1.0 when rounding occurs in addition, 0.0 otherwise.
!!
!!   = 'M' or 'm',   mobbrmsd_SLAMCH := EMIN;
!!   minimum exponent before (gradual) underflow
!!
!!   = 'U' or 'u',   mobbrmsd_SLAMCH := RMIN;
!!   underflow threshold - base**(emin-1)
!!
!!   = 'L' or 'l',   mobbrmsd_SLAMCH := EMAX;
!!   largest exponent before overflow
!!
!!   = 'O' or 'o',   mobbrmsd_SLAMCH := RMAX;
!!   overflow threshold  - (base**emax)*(1-eps)
!!
  real(RK)  :: mobbrmsd_SLAMCH
!! machine parameter.
!!
  intrinsic :: DIGITS, EPSILON, HUGE, MAXEXPONENT, MINEXPONENT, RADIX, TINY
  real(RK), parameter :: RADIX_ZERO = RADIX(ZERO)
  real(RK), parameter :: TINY_ZERO = TINY(ZERO)
  real(RK), parameter :: DIGITS_ZERO = DIGITS(ZERO)
  real(RK), parameter :: MINEXPONENT_ZERO = MINEXPONENT(ZERO)
  real(RK), parameter :: MAXEXPONENT_ZERO = MAXEXPONENT(ZERO)
  real(RK), parameter :: HUGE_ZERO = HUGE(ZERO)
  real(RK), parameter :: SMALL = ONE / HUGE_ZERO
!
! Use SMALL plus a bit, to avoid the possibility of rounding
! causing overflow when computing  1/SFMIN.
!
  real(RK), parameter :: SFMIN = MERGE(SMALL, TINY_ZERO, SMALL >= TINY_ZERO)
!
  if (mobbrmsd_LSAME(CMACH, 'E')) then
    mobbrmsd_SLAMCH = ULP
  else if (mobbrmsd_LSAME(CMACH, 'S')) then
    mobbrmsd_SLAMCH = SFMIN
  else if (mobbrmsd_LSAME(CMACH, 'B')) then
    mobbrmsd_SLAMCH = RADIX_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'P')) then
    mobbrmsd_SLAMCH = ULP ! RADIX(ZERO)
  else if (mobbrmsd_LSAME(CMACH, 'N')) then
    mobbrmsd_SLAMCH = DIGITS_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'R')) then
    mobbrmsd_SLAMCH = ONE
  else if (mobbrmsd_LSAME(CMACH, 'M')) then
    mobbrmsd_SLAMCH = MINEXPONENT_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'U')) then
    mobbrmsd_SLAMCH = TINY_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'L')) then
    mobbrmsd_SLAMCH = MAXEXPONENT_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'O')) then
    mobbrmsd_SLAMCH = HUGE_ZERO
  else
    mobbrmsd_SLAMCH = ZERO
  end if
!
  return
!
! End of mobbrmsd_SLAMCH
!
end


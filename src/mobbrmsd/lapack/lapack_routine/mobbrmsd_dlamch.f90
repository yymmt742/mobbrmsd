!| determines double precision machine parameters.
!
!  mobbrmsd_DLAMCH determines double precision machine parameters.
!  Assume rounding, not chopping. Always.
!
!  Reference DLAMCH is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- LAPACK auxiliary routine (version 3.3.0) --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     Based on LAPACK mobbrmsd_DLAMCH but with Fortran 95 query functions
!     See [documentation](http://www.cs.utk.edu/~luszczek/lapack/lamch.html)
!     and [netlib](http://www.netlib.org/lapack-dev/lapack-coding/program-style.html#id2537289).
!     July 2010
!
pure elemental function mobbrmsd_DLAMCH(CMACH)
  character, intent(in) :: CMACH
!!  Specifies the value to be returned by mobbrmsd_DLAMCH:
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
  real(RK)            :: mobbrmsd_DLAMCH
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
!     Use SMALL plus a bit, to avoid the possibility of rounding
!     causing overflow when computing  1/SFMIN.
!
  real(RK), parameter :: SFMIN = MERGE(SMALL, TINY_ZERO, SMALL >= TINY_ZERO)
!
  if (mobbrmsd_LSAME(CMACH, 'E')) then
    mobbrmsd_DLAMCH = ULP
  else if (mobbrmsd_LSAME(CMACH, 'S')) then
    mobbrmsd_DLAMCH = SFMIN
  else if (mobbrmsd_LSAME(CMACH, 'B')) then
    mobbrmsd_DLAMCH = RADIX_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'P')) then
    mobbrmsd_DLAMCH = ULP ! RADIX(ZERO)
  else if (mobbrmsd_LSAME(CMACH, 'N')) then
    mobbrmsd_DLAMCH = DIGITS_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'R')) then
    mobbrmsd_DLAMCH = ONE
  else if (mobbrmsd_LSAME(CMACH, 'M')) then
    mobbrmsd_DLAMCH = MINEXPONENT_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'U')) then
    mobbrmsd_DLAMCH = TINY_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'L')) then
    mobbrmsd_DLAMCH = MAXEXPONENT_ZERO
  else if (mobbrmsd_LSAME(CMACH, 'O')) then
    mobbrmsd_DLAMCH = HUGE_ZERO
  else
    mobbrmsd_DLAMCH = ZERO
  end if
!
  return
!
! End of mobbrmsd_DLAMCH
!
end function mobbrmsd_DLAMCH


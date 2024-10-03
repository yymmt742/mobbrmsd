pure elemental function DLAMCH(CMACH)
! use LA_CONSTANTS, only: wp => DP, ONE=>DONE, ZERO=>DZERO, EPS => DULP
!
!  -- LAPACK auxiliary routine (version 3.3.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     Based on LAPACK DLAMCH but with Fortran 95 query functions
!     See: http://www.cs.utk.edu/~luszczek/lapack/lamch.html
!     and  http://www.netlib.org/lapack-dev/lapack-coding/program-style.html#id2537289
!     July 2010
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: CMACH
  real(RK)                 :: DLAMCH
!     ..
!
!  Purpose
!  =======
!
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER!1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps!base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps!base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base!!(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base!!emax)!(1-eps)
!
! =====================================================================
!     ..
!     .. Local Scalars ..
  real(RK) :: RND, SFMIN, SMALL, RMACH
!     ..
!     .. External Functions ..
! interface
!   include 'lsame.h'
! end interface
!     ..
!     .. Intrinsic Functions ..
  intrinsic :: DIGITS, EPSILON, HUGE, MAXEXPONENT, &
 &             MINEXPONENT, RADIX, TINY
!     ..
!     .. Executable Statements ..
!
!
!     Assume rounding, not chopping. Always.
!
  RND = ONE
!
! if (ONE == RND) then
!   EPS = EPSILON(ZERO) ! 0.5
! else
!   EPS = EPSILON(ZERO)
! end if
!
  if (LSAME(CMACH, 'E')) then
    RMACH = DULP
  else if (LSAME(CMACH, 'S')) then
    SFMIN = TINY(ZERO)
    SMALL = ONE / HUGE(ZERO)
    if (SMALL >= SFMIN) then
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
      SFMIN = SMALL!( ONE+EPS )
    end if
    RMACH = SFMIN
  else if (LSAME(CMACH, 'B')) then
    RMACH = RADIX(ZERO)
  else if (LSAME(CMACH, 'P')) then
    RMACH = DULP ! RADIX(ZERO)
  else if (LSAME(CMACH, 'N')) then
    RMACH = DIGITS(ZERO)
  else if (LSAME(CMACH, 'R')) then
    RMACH = RND
  else if (LSAME(CMACH, 'M')) then
    RMACH = MINEXPONENT(ZERO)
  else if (LSAME(CMACH, 'U')) then
    RMACH = TINY(zero)
  else if (LSAME(CMACH, 'L')) then
    RMACH = MAXEXPONENT(ZERO)
  else if (LSAME(CMACH, 'O')) then
    RMACH = HUGE(ZERO)
  else
    RMACH = ZERO
  end if
!
  DLAMCH = RMACH
  return
!
!     End of DLAMCH
!
end function DLAMCH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
pure elemental function DLAMC3(A, B)
! use LA_CONSTANTS, only: wp => DP
!
!  -- LAPACK auxiliary routine (version 3.3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
  real(RK), intent(in) :: A, B
  real(RK)             :: DLAMC3
!     ..
!
!  Purpose
!  =======
!
!  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!  B       (input) DOUBLE PRECISION
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
  DLAMC3 = A + B
!
  return
!
!     End of DLAMC3
!
end function DLAMC3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

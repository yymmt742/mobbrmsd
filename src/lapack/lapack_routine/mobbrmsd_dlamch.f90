pure elemental function mobbrmsd_DLAMCH(CMACH)
! use LA_CONSTANTS, only: wp => DP, ONE=>DONE, ZERO=>DZERO, EPS => ULP
!
!  -- LAPACK auxiliary routine (version 3.3.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     Based on LAPACK mobbrmsd_DLAMCH but with Fortran 95 query functions
!     See: http://www.cs.utk.edu/~luszczek/lapack/lamch.html
!     and  http://www.netlib.org/lapack-dev/lapack-coding/program-style.html#id2537289
!     July 2010
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: CMACH
  real(RK)                 :: mobbrmsd_DLAMCH
!     ..
!
!  Purpose
!  =======
!
!  mobbrmsd_DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER!1
!          Specifies the value to be returned by mobbrmsd_DLAMCH:
!          = 'E' or 'e',   mobbrmsd_DLAMCH := eps
!          = 'S' or 's ,   mobbrmsd_DLAMCH := sfmin
!          = 'B' or 'b',   mobbrmsd_DLAMCH := base
!          = 'P' or 'p',   mobbrmsd_DLAMCH := eps!base
!          = 'N' or 'n',   mobbrmsd_DLAMCH := t
!          = 'R' or 'r',   mobbrmsd_DLAMCH := rnd
!          = 'M' or 'm',   mobbrmsd_DLAMCH := emin
!          = 'U' or 'u',   mobbrmsd_DLAMCH := rmin
!          = 'L' or 'l',   mobbrmsd_DLAMCH := emax
!          = 'O' or 'o',   mobbrmsd_DLAMCH := rmax
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
  if (mobbrmsd_LSAME(CMACH, 'E')) then
    RMACH = ULP
  else if (mobbrmsd_LSAME(CMACH, 'S')) then
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
  else if (mobbrmsd_LSAME(CMACH, 'B')) then
    RMACH = RADIX(ZERO)
  else if (mobbrmsd_LSAME(CMACH, 'P')) then
    RMACH = ULP ! RADIX(ZERO)
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
  mobbrmsd_DLAMCH = RMACH
  return
!
!     End of mobbrmsd_DLAMCH
!
end function mobbrmsd_DLAMCH
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

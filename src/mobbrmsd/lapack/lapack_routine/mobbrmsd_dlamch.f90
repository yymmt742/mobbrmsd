!| mobbrmsd_DLAMCH determines double precision machine parameters.
!  Assume rounding, not chopping. Always.
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
!!          Specifies the value to be returned by mobbrmsd_DLAMCH:
!1
!!          = 'E' or 'e',   mobbrmsd_DLAMCH := eps
!1
!!          = 'S' or 's ,   mobbrmsd_DLAMCH := sfmin
!1
!!          = 'B' or 'b',   mobbrmsd_DLAMCH := base
!1
!!          = 'P' or 'p',   mobbrmsd_DLAMCH := eps!base
!1
!!          = 'N' or 'n',   mobbrmsd_DLAMCH := t
!1
!!          = 'R' or 'r',   mobbrmsd_DLAMCH := rnd
!1
!!          = 'M' or 'm',   mobbrmsd_DLAMCH := emin
!1
!!          = 'U' or 'u',   mobbrmsd_DLAMCH := rmin
!1
!!          = 'L' or 'l',   mobbrmsd_DLAMCH := emax
!1
!!          = 'O' or 'o',   mobbrmsd_DLAMCH := rmax
!!
!!          where
!!
!!          eps   = relative machine precision
!1
!!          sfmin = safe minimum, such that 1/sfmin does not overflow
!1
!!          base  = base of the machine
!1
!!          prec  = eps!base
!1
!!          t     = number of (base) digits in the mantissa
!1
!!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!1
!!          emin  = minimum exponent before (gradual) underflow
!1
!!          rmin  = underflow threshold - base!!(emin-1)
!1
!!          emax  = largest exponent before overflow
!1
!!          rmax  = overflow threshold  - (base!!emax)!(1-eps)
!1
  real(RK)                 :: mobbrmsd_DLAMCH
!
  real(RK) :: RND, SFMIN, SMALL, RMACH
  intrinsic :: DIGITS, EPSILON, HUGE, MAXEXPONENT, MINEXPONENT, RADIX, TINY
! interface
!   include 'lsame.h'
! end interface
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


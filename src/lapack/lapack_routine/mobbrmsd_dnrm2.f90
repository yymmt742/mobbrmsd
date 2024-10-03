!> \brief \b DNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    DNRM2 := sqrt( x'*x )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER, storage spacing between elements of X
!>          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!>          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!>          If INCX = 0, x isn't a vector so there is no need to call
!>          this subroutine.  If you call it anyway, it will count x(1)
!>          in the vector norm N times.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup single_blas_level1
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!>
!  =====================================================================
pure function DNRM2(n, x, incx)
! use LA_CONSTANTS, only: RK => dp
!  ..
!  .. Scalar Arguments ..
  integer, intent(in)  :: incx, n
!  ..
!  .. Array Arguments ..
  real(RK), intent(in) :: x(*)
!
  real(RK) :: DNRM2
!
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  .. Constants ..
! real(RK), parameter :: ZERO = 0.0_RK
! real(RK), parameter :: ONE = 1.0_RK
  real(RK), parameter :: MAXN = HUGE(0.0_RK)
!  ..
!  .. Blue's scaling constants ..
! real(RK), parameter :: TSML = real(RADIX(0._RK), RK)**CEILING( &
!                      & (MINEXPONENT(0._RK) - 1) * 0.5_RK)
! real(RK), parameter :: TBIG = real(RADIX(0._RK), RK)**FLOOR( &
!                      & (MAXEXPONENT(0._RK) - DIGITS(0._RK) + 1) * 0.5_RK)
! real(RK), parameter :: SSML = real(RADIX(0._RK), RK)**(-FLOOR( &
!                      & (MINEXPONENT(0._RK) - DIGITS(0._RK)) * 0.5_RK))
! real(RK), parameter :: SBIG = real(RADIX(0._RK), RK)**(-CEILING( &
!                      & (MAXEXPONENT(0._RK) + DIGITS(0._RK) - 1) * 0.5_RK))
!  ..
!  .. Local Scalars ..
  integer :: i, ix
  logical :: notbig
  real(RK) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
!
!  Quick return if possible
!
  DNRM2 = ZERO
  if (n <= 0) return
!
  scl = ONE
  sumsq = ZERO
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
  notbig = .true.
  asml = zero
  amed = zero
  abig = zero
  ix = 1
  if (incx < 0) ix = 1 - (n - 1) * incx
  do i = 1, n
    ax = ABS(x(ix))
    if (ax > tbig) then
      abig = abig + (ax * sbig)**2
      notbig = .false.
    else if (ax < tsml) then
      if (notbig) asml = asml + (ax * ssml)**2
    else
      amed = amed + ax**2
    end if
    ix = ix + incx
  end do
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
  if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
    if ((amed > zero) .or. (amed > maxN) .or. (amed /= amed)) then
      abig = abig + (amed * sbig) * sbig
    end if
    scl = one / sbig
    sumsq = abig
  else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
    if ((amed > zero) .or. (amed > maxN) .or. (amed /= amed)) then
      amed = SQRT(amed)
      asml = SQRT(asml) / ssml
      if (asml > amed) then
        ymin = amed
        ymax = asml
      else
        ymin = asml
        ymax = amed
      end if
      scl = one
      sumsq = ymax**2 * (one + (ymin / ymax)**2)
    else
      scl = one / ssml
      sumsq = asml
    end if
  else
!
!     Otherwise all values are mid-range
!
    scl = one
    sumsq = amed
  end if
  DNRM2 = scl * SQRT(sumsq)
  return
end function


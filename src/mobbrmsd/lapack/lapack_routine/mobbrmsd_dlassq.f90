!| mobbrmsd_DLASSQ updates a sum of squares represented in scaled form.
!
!  mobbrmsd_DLASSQ  returns the values  scl  and  smsq  such that
!
!    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative.
!
!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
!  If scale * sqrt( sumsq ) > tbig then
!
!     we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
!
!  and if 0 < scale * sqrt( sumsq ) < tsml then
!
!     we require:   scale <= sqrt( HUGE ) / ssml       on entry,
!
!  where
!
!     tbig -- upper threshold for values whose square is representable;
!
!     sbig -- scaling constant for big numbers; \see la_constants.f90
!
!     tsml -- lower threshold for values whose square is representable;
!
!     ssml -- scaling constant for small numbers; \see la_constants.f90
!
!  and
!
!     TINY*EPS -- tiniest representable number;
!
!     HUGE     -- biggest representable number.
!
!   Anderson E. (2017)
!   Algorithm 978: Safe Scaling in the Level 1 BLAS
!   [ACM Trans Math Softw 44:1--28](https://doi.org/10.1145/3061665)
!
!   Blue, James L. (1978)
!   A Portable Fortran Program to Find the Euclidean Norm of a Vector
!   [ACM Trans Math Softw 4:15--23](https://doi.org/10.1145/355769.355771)
!
!  Reference DLASSQ is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK auxiliary routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     Edward Anderson, Lockheed Martin
!     Weslley Pereira, University of Colorado Denver, USA
!     Nick Papior, Technical University of Denmark, DK
!
pure subroutine mobbrmsd_DLASSQ(n, x, incx, scl, sumsq)
  integer, intent(IN)     :: N
!!   The number of elements to be used from the vector x.
!!
  real(RK), intent(IN)    :: X(*)
!!   DOUBLE PRECISION array, dimension (1+(N-1)*abs(INCX))
!!   The vector for which a scaled sum of squares is computed.
!!      x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!!
  integer, intent(IN)     :: INCX
!!   The increment between successive values of the vector x.
!!
!!   If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!!
!!   If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!!
!!   If INCX = 0, x isn't a vector so there is no need to call
!!   this subroutine.  If you call it anyway, it will count x(1)
!!   in the vector norm N times.
!!
  real(RK), intent(INOUT) :: SCL
!!   On entry, the value  scale  in the equation above.
!!
!!   On exit, SCALE is overwritten with  scl , the scaling factor
!!   for the sum of squares.
!!
  real(RK), intent(INOUT) :: SUMSQ
!!   On entry, the value  sumsq  in the equation above.
!!
!!   On exit, SUMSQ is overwritten with  smsq , the basic sum of
!!   squares from which  scl  has been factored out.
!!
  integer                :: I, IX
  logical                :: NOTBIG
  real(RK)               :: ABIG, AMED, ASML, AX, YMAX, YMIN
!
!  Quick return if possible
!
  if (IEEE_IS_NAN(scl) .or. IEEE_IS_NAN(sumsq)) return
  if (sumsq == zero) scl = one
  if (scl == zero) then
    scl = one
    sumsq = zero
  end if
  if (n <= 0) then
    return
  end if
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
    if (ax > TBIG) then
      abig = abig + (ax * SBIG)**2
      notbig = .false.
    else if (ax < TSML) then
      if (notbig) asml = asml + (ax * SSML)**2
    else
      amed = amed + ax**2
    end if
    ix = ix + incx
  end do
!
!  Put the existing sum of squares into one of the accumulators
!
  if (sumsq > zero) then
    ax = scl * SQRT(sumsq)
    if (ax > TBIG) then
!        We assume scl >= sqrt( TINY*EPS ) / SBIG
      abig = abig + (scl * SBIG)**2 * sumsq
    else if (ax < TSML) then
!        We assume scl <= sqrt( HUGE ) / SSML
      if (notbig) asml = asml + (scl * SSML)**2 * sumsq
    else
      amed = amed + scl**2 * sumsq
    end if
  end if
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
  if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
    if (amed > zero .or. IEEE_IS_NAN(amed)) then
      abig = abig + (amed * SBIG) * SBIG
    end if
    scl = one / SBIG
    sumsq = abig
  else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
    if (amed > zero .or. IEEE_IS_NAN(amed)) then
      amed = SQRT(amed)
      asml = SQRT(asml) / SSML
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
      scl = one / SSML
      sumsq = asml
    end if
  else
!
!     Otherwise all values are mid-range or zero
!
    scl = one
    sumsq = amed
  end if
  return
end subroutine mobbrmsd_DLASSQ


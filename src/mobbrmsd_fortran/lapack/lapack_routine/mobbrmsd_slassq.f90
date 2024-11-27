!| mobbrmsd_SLASSQ updates a sum of squares represented in scaled form.
!
!  mobbrmsd_SLASSQ  returns the values  scl  and  smsq  such that
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
!  ReferenceSDLASSQ is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     Edward Anderson, Lockheed Martin
!     Weslley Pereira, University of Colorado Denver, USA
!     Nick Papior, Technical University of Denmark, DK
!
pure subroutine MOBBRMSD_SLASSQ(N, X, INCX, SCL, SUMSQ)
  integer, intent(IN)     :: N
!!   The number of elements to be used from the vector x.
!!
  real(RK), intent(IN)    :: X(*)
!!   REAL array, dimension (1+(N-1)*abs(INCX))
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
  integer   :: IX
  real      :: ABSXI
  intrinsic :: ABS
!..
!..Executable Statements..
!
  if (N > 0) then
    do IX = 1, 1 + (N - 1) * INCX, INCX
      ABSXI = ABS(X(IX))
      if (ABSXI > ZERO .or. IEEE_IS_NAN(ABSXI)) then
        if (SCL < ABSXI) then
          SUMSQ = 1 + SUMSQ * (SCL / ABSXI)**2
          SCL = ABSXI
        else
          SUMSQ = SUMSQ + (ABSXI / SCL)**2
        end if
      end if
    end do
  end if
  return
!
! end of mobbrmsd_SLASSQ
!
end


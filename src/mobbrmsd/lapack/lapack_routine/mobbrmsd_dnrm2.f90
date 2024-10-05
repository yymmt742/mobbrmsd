!| mobbrmsd_DNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!  \( \text{mobbrmsd_DNRM2} := \sqrt{X X^{\top}} \)
!
!  COMPUTE THE SUM OF SQUARES IN 3 ACCUMULATORS:
!
!      ABIG -- SUMS OF SQUARES SCALED DOWN TO AVOID OVERFLOW
!      ASML -- SUMS OF SQUARES SCALED UP TO AVOID UNDERFLOW
!      AMED -- SUMS OF SQUARES THAT DO NOT REQUIRE SCALING
!
!  THE THRESHOLDS AND MULTIPLIERS ARE
!
!      TBIG -- VALUES BIGGER THAN THIS ARE SCALED DOWN BY SBIG
!      TSML -- VALUES SMALLER THAN THIS ARE SCALED UP BY SSML
!
!  See details
!
!      Anderson E. (2017)
!      Algorithm 978: Safe Scaling in the Level 1 BLAS
!      [ACM Trans Math Softw 44:1--28](https://doi.org/10.1145/3061665),
!
!      Blue, James L. (1978)
!      A Portable Fortran Program to Find the Euclidean Norm of a Vector
!      [ACM Trans Math Softw 4:15--23](https://doi.org/10.1145/355769.355771)
!
! > Reference DNRM2 is provided by [netlib.org](http://www.netlib.org/lapack/).
! > Reference BLAS level1 routine (version 3.9.1)
! > Reference BLAS is a software package provided by Univ. of Tennessee,
! > Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
! > March 2021
! > Weslley Pereira, University of Colorado Denver, USA
! > Edward Anderson, Lockheed Martin
!
pure function mobbrmsd_DNRM2(N, X, INCX)
  integer, intent(in)  :: N
!!         number of elements in input vector(s)
!!
  real(RK), intent(in) :: X(*)
!!          DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!!
  integer, intent(in)  :: INCX
!!          INCX is INTEGER, storage spacing between elements of X
!!
!!          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!!
!!          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!!
!!          If INCX = 0, x isn't a vector so there is no need to call
!!          this subroutine.  If you call it anyway, it will count x(1)
!!          in the vector norm N times.
!!
  real(RK) :: mobbrmsd_DNRM2
!! Euclidean norm, \( \sqrt{X X^{\top}} \).
!!
  real(RK), parameter :: MAXN = HUGE(0.0_RK)
  integer  :: I, IX
  logical  :: NOTBIG
  real(RK) :: ABIG, AMED, ASML, AX, SCL, SUMSQ, YMAX, YMIN
!
!  QUICK RETURN IF POSSIBLE
!
  mobbrmsd_DNRM2 = ZERO
  if (N <= 0) return
!
  SCL = ONE
  SUMSQ = ZERO
!
  NOTBIG = .true.
  ASML = ZERO
  AMED = ZERO
  ABIG = ZERO
  IX = 1
  if (INCX < 0) IX = 1 - (N - 1) * INCX
  do I = 1, N
    AX = ABS(X(IX))
    if (AX > TBIG) then
      ABIG = ABIG + (AX * SBIG)**2
      NOTBIG = .false.
    else if (AX < TSML) then
      if (NOTBIG) ASML = ASML + (AX * SSML)**2
    else
      AMED = AMED + AX**2
    end if
    IX = IX + INCX
  end do
!
!  COMBINE ABIG AND AMED OR AMED AND ASML IF MORE THAN ONE
!  ACCUMULATOR WAS USED.
!
  if (ABIG > ZERO) then
!
!     COMBINE ABIG AND AMED IF ABIG > 0.
!
    if ((AMED > ZERO) .or. (AMED > MAXN) .or. (AMED /= AMED)) then
      ABIG = ABIG + (AMED * SBIG) * SBIG
    end if
    SCL = ONE / SBIG
    SUMSQ = ABIG
  else if (ASML > ZERO) then
!
!     COMBINE AMED AND ASML IF ASML > 0.
!
    if ((AMED > ZERO) .or. (AMED > MAXN) .or. (AMED /= AMED)) then
      AMED = SQRT(AMED)
      ASML = SQRT(ASML) / SSML
      if (ASML > AMED) then
        YMIN = AMED
        YMAX = ASML
      else
        YMIN = ASML
        YMAX = AMED
      end if
      SCL = ONE
      SUMSQ = YMAX**2 * (ONE + (YMIN / YMAX)**2)
    else
      SCL = ONE / SSML
      SUMSQ = ASML
    end if
  else
!
!     OTHERWISE ALL VALUES ARE MID-RANGE
!
    SCL = ONE
    SUMSQ = AMED
  end if
  mobbrmsd_DNRM2 = SCL * SQRT(SUMSQ)
  return
end


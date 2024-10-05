!| mobbrmsd_DLAPY2 returns \( \sqrt{x^2 + y^2} \),
!  taking care not to cause unnecessary overflow
!  and unnecessary underflow.
!
!     Reference DLAPY2 is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- LAPACK driver routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure elemental function mobbrmsd_DLAPY2(X, Y)
  implicit none
  real(RK), intent(in) :: X
!!          X specify the value x.
!!
  real(RK), intent(in) :: Y
!!          Y specify the value y.
!!
  real(RK)             :: mobbrmsd_DLAPY2
!!  \( \sqrt{x^2 + y^2} \)
!!
  real(RK)             :: W, XABS, YABS, Z, HUGEVAL
  logical              :: X_IS_NAN, Y_IS_NAN
  intrinsic            :: ABS, MAX, MIN, SQRT
! real(RK), parameter  :: ZERO = 0.0_RK
! real(RK), parameter  :: ONE = 1.0_RK
! interface
!   include 'dlamch.h'
! end interface
!
  X_IS_NAN = IEEE_IS_NAN(X)
  Y_IS_NAN = IEEE_IS_NAN(Y)
!
  if (.not. (X_IS_NAN .or. Y_IS_NAN)) then
    XABS = ABS(X)
    YABS = ABS(Y)
    W = MAX(XABS, YABS)
    Z = MIN(XABS, YABS)
    HUGEVAL = mobbrmsd_DLAMCH('Overflow')
    if (Z == ZERO .or. W > HUGEVAL) then
      mobbrmsd_DLAPY2 = W
    else
      mobbrmsd_DLAPY2 = W * SQRT(ONE + (Z / W)**2)
    end if
  else
    if (X_IS_NAN) then
      mobbrmsd_DLAPY2 = X
    else
      mobbrmsd_DLAPY2 = Y
    end if
  end if
  return
!
!     End of mobbrmsd_DLAPY2
!
end function mobbrmsd_DLAPY2


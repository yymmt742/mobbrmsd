!| mobbrmsd_SLAPY2 returns \( \sqrt{x^2 + y^2} \),
!  taking care not to cause unnecessary overflow
!  and unnecessary underflow.
!
!     Reference SLAPY2 is provided by [netlib](http://www.netlib.org/lapack/)
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
pure elemental function mobbrmsd_SLAPY2(X, Y)
  implicit none
  real(RK), intent(in) :: X
!!          X specify the value x.
!!
  real(RK), intent(in) :: Y
!!          Y specify the value y.
!!
  real(RK) :: mobbrmsd_SLAPY2
!!  \( \sqrt{x^2 + y^2} \)
!!
  real(RK) :: W, XABS, YABS, Z
  logical  :: X_IS_NAN, Y_IS_NAN
  intrinsic :: ABS, MAX, MIN, SQRT
! real(RK), parameter :: ZERO = 0.0E0
! real(RK), parameter :: ONE = 1.0E0
!
  X_IS_NAN = IEEE_IS_NAN(X)
  Y_IS_NAN = IEEE_IS_NAN(Y)
!
  if (.not. (X_IS_NAN .or. Y_IS_NAN)) then
    XABS = ABS(X)
    YABS = ABS(Y)
    W = MAX(XABS, YABS)
    Z = MIN(XABS, YABS)
    if (Z == ZERO) then
      mobbrmsd_SLAPY2 = W
    else
      mobbrmsd_SLAPY2 = W * SQRT(ONE + (Z / W)**2)
    end if
  else
    if (X_IS_NAN) then
      mobbrmsd_SLAPY2 = X
    else
      mobbrmsd_SLAPY2 = Y
    end if
  end if
  return
!
!end of mobbrmsd_SLAPY2
!
end


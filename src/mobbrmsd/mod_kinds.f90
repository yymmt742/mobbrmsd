!| variable kinds correction
module mod_kinds
!
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                R4     => REAL32,      &
    &                R8     => REAL64,      &
    &                RQ     => REAL128,     &
    &                I1     => INT8,        &
    &                I2     => INT16,       &
    &                I4     => INT32,       &
    &                I8     => INT64
!
  implicit none
  private
  public  :: I1, I2, I4, I8
  public  :: R4, R8, RQ
  public  :: RK, IK
!
  integer, parameter  :: IK = KIND(0)
  !! Selected integer kind.
  integer, parameter  :: RK = KIND(0.0_R8)
  !! Selected real kind.
!
end module mod_kinds


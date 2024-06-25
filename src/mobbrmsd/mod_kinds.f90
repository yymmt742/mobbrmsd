!| Define variable kinds correction
module mod_kinds
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                I1 => INT8,        &
    &                I2 => INT16,       &
    &                I4 => INT32,       &
    &                I8 => INT64
  implicit none
  private
  public :: I1, I2, I4, I8
  public :: R4, R8
  public :: RK, IK
  !| Selected integer kind.
#ifdef USE_INT8
  integer, parameter  :: IK = I1
#elif USE_INT16
  integer, parameter  :: IK = I2
#elif USE_INT32
  integer, parameter  :: IK = I4
#elif USE_INT64
  integer, parameter  :: IK = I8
#else
  integer, parameter  :: IK = KIND(0)
#endif
  !| Single precision kind.
  integer, parameter  :: R4 = SELECTED_REAL_KIND(6)
  !| Double precision kind.
  integer, parameter  :: R8 = SELECTED_REAL_KIND(15)
  !| Selected real kind.
#ifdef USE_REAL32
  integer, parameter  :: RK = R4
#elif USE_REAL64
  integer, parameter  :: RK = R8
#else
  integer, parameter  :: RK = R8
#endif
end module mod_kinds


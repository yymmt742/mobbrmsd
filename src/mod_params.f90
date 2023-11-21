module mod_params
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                STDIN => INPUT_UNIT,  &
    &                STDOUT => OUTPUT_UNIT, &
    &                STDERR => ERROR_UNIT,  &
    &                R4 => REAL32,  &
    &                R8 => REAL64,  &
    &                RQ => REAL128, &
    &                I1 => INT8, &
    &                I2 => INT16, &
    &                I4 => INT32, &
    &                I8 => INT64
  implicit none
  private
  public  :: R4, I8, RQ
  public  :: I1, I2, I4, I8
  public  :: RK, IK, LK
  public  :: STDIN, STDOUT, STDERR
  public  :: RZERO, RONE, RHALF, RFOUR, RHUGE
!
!&<
!
  integer, parameter          :: IK = I4
  integer, parameter          :: RK = R8
  integer, parameter          :: LK = KIND(.true.)
!
  real(RK), parameter         :: RZERO = 0.0_RK
  real(RK), parameter         :: RONE  = 1.0_RK
  real(RK), parameter         :: RHALF = 0.5_RK
  real(RK), parameter         :: RFOUR = 4.0_RK
!
  real(RK), parameter         :: RHUGE = HUGE(RZERO)
!
!&>
!
end module mod_params

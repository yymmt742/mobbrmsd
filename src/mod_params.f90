module mod_params
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                STDIN => INPUT_UNIT,  &
    &                STDOUT => OUTPUT_UNIT, &
    &                STDERR => ERROR_UNIT,  &
    &                RK => REAL64,  &
    &                IK => INT32
  implicit none
  private
  public  :: RK, IK, LK
  public  :: STDIN, STDOUT, STDERR
  public  :: RZERO, RONE, RHALF
!
!&<
!
  integer(IK), parameter      :: LK = KIND(.true.)
!
  real(RK), parameter         :: RZERO = 0.0_RK
  real(RK), parameter         :: RONE  = 1.0_RK
  real(RK), parameter         :: RHALF = 0.5_RK
!
!&>
!
end module mod_params

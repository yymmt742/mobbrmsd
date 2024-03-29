!| Module with constant collection.
module mod_params
!
  use :: mod_kinds, only : I1, I2, I4, I8, R4, R8, RQ, RK, IK
  use :: blas_lapack_interface, only : D, DD, ND, setup_dimension, gemm => DGEMM
!
  implicit none
  private
  public  :: I1, I2, I4, I8
  public  :: R4, R8, RQ
  public  :: RK, IK
  public  :: D, DD, ND
  public  :: RZERO, RONE, RHALF, RFOUR, RTEN
  public  :: RPI, RNAPIER, RHUGE, LN_TO_L10
  public  :: setup_dimension, gemm
!
!&<
!
  real(RK), parameter :: RZERO = 0.0_RK
  !! Real zero.
  real(RK), parameter :: RONE  = 1.0_RK
  !! Real one.
  real(RK), parameter :: RHALF = 0.5_RK
  !! Real 1/2.
  real(RK), parameter :: RFOUR = 4.0_RK
  !! Real four.
  real(RK), parameter :: RTEN  = 10.0_RK
  !! Real ten.

  real(RK), parameter :: RPI       = ACOS(RZERO)
  !! Real circular constant.
  real(RK), parameter :: RNAPIER   = 2.71828182846_RK
  !! Napier constant.
  real(RK), parameter :: RHUGE     = HUGE(RZERO)
  !! Real large numbe
  real(RK), parameter :: LN_TO_L10 = LOG10(RNAPIER)
  !! Scaling factor, LOG10(x) =LN_TO_L10 * LN(x)
!
!&>
!
end module mod_params


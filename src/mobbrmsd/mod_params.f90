!| Module for blas lapack interface
module blas_lapack_interface
  implicit none
  private
  public  :: SCOPY,  DCOPY
  public  :: SGEMM,  DGEMM
  public  :: SGESVD, DGESVD
  public  :: SGETRF, DGETRF
  public  :: D, DD, ND
  public  :: setup_dimension
!
  interface
!
    include 'dcopy.h'
    include 'scopy.h'
!
    include 'dgemm.h'
    include 'sgemm.h'
!
    include 'dgesvd.h'
    include 'sgesvd.h'
!
    include 'dgetrf.h'
    include 'sgetrf.h'
!
  end interface
!
  integer, save, protected :: D  = 3
  !! Spatial dimension
  integer, save, protected :: DD = 9
  !! Square spatial dimension
  integer, save, protected :: ND = 9 + 2
  !! Node size, defined by [L, G, C(D,D)]
!
contains
!
!| Sets the dimensions of the space. <br>
!  Caution, this routine affects global.
  subroutine setup_dimension(d_)
    integer, intent(in) :: d_
    D = MAX(1, d_)
    DD = D * D
    ND = DD + 2
  end subroutine setup_dimension
!
end module blas_lapack_interface
!
!| Module with constant collection.
module mod_params
!
!&<
!
  use :: mod_kinds, only : I1, I2, I4, I8, R4, R8, RQ, RK, IK
  use :: blas_lapack_interface, only : D, DD, ND, setup_dimension
!
  use :: blas_lapack_interface, only : &
    &                copy   => DCOPY,  &
    &                gemm   => DGEMM,  &
    &                gesvd  => DGESVD, &
    &                getrf  => DGETRF
!
! use :: blas_lapack_interface, only : &
!   &                copy   => SCOPY,  &
!   &                gemm   => SGEMM,  &
!   &                gesvd  => SGESVD, &
!   &                getrf  => SGETRF
!
!&>
!
  implicit none
  private
  public  :: I1, I2, I4, I8
  public  :: R4, R8, RQ
  public  :: RK, IK
  public  :: D, DD, ND
  public  :: RZERO, RONE, RHALF, RFOUR, RTEN
  public  :: RNAPIER, RHUGE, LN_TO_L10
!
  public  :: setup_dimension
  public  :: copy, gemm, gesvd, getrf
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
!
  real(RK), parameter :: RNAPIER = 2.71828182846_RK
  !! Napier constant.
  real(RK), parameter :: RHUGE = HUGE(RZERO)
  !! Real large numbe
  real(RK), parameter :: LN_TO_L10 = LOG10(RNAPIER)
  !! Scaling factor, LOG10(x) =LN_TO_L10 * LN(x)
!
!&>
!
end module mod_params


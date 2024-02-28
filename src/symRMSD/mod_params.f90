!| Module for blas lapack interface
module blas_lapack_interface
  implicit none
  private
  public  :: SCOPY,  DCOPY
  public  :: SAXPY,  DAXPY
  public  :: SDOT,   DDOT
  public  :: SSYMM,  DSYMM
  public  :: SGEMM,  DGEMM
  public  :: SGESVD, DGESVD
  public  :: SGETRF, DGETRF
!
  interface
!
    include 'daxpy.h'
    include 'saxpy.h'
!
    include 'ddot.h'
    include 'sdot.h'
!
    include 'dcopy.h'
    include 'scopy.h'
!
    include 'dsymm.h'
    include 'ssymm.h'
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
end module blas_lapack_interface
!
!| Module with constant collection.
module mod_params
!
!&<
!
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                R4     => REAL32,      &
    &                R8     => REAL64,      &
    &                RQ     => REAL128,     &
    &                I1     => INT8,        &
    &                I2     => INT16,       &
    &                I4     => INT32,       &
    &                I8     => INT64,       &
    &                STDIN  => INPUT_UNIT,  &
    &                STDOUT => OUTPUT_UNIT, &
    &                STDERR => ERROR_UNIT
!
  use :: blas_lapack_interface, only : &
    &                copy   => DCOPY,       &
    &                axpy   => DAXPY,       &
    &                dot    => DDOT,        &
    &                symm   => DSYMM,       &
    &                gemm   => DGEMM,       &
    &                gesvd  => DGESVD,      &
    &                getrf  => DGETRF
!
! use :: blas_lapack_interface, only : &
!   &                copy   => SCOPY,       &
!   &                axpy   => SAXPY,       &
!   &                dot    => SDOT,        &
!   &                symm   => SSYMM,       &
!   &                gemm   => SGEMM,       &
!   &                gesvd  => SGESVD,      &
!   &                getrf  => SGETRF
!
!&>
!
  implicit none
  private
  public  :: I1, I2, I4, I8
  public  :: R4, R8, RQ
  public  :: RK, IK, LK
  public  :: STDIN, STDOUT, STDERR
  public  :: RZERO, RONE, RHALF, RFOUR, RTEN
  public  :: RNAPIER, RHUGE, LN_TO_L10
!
  public  :: D, DD
  public  :: setup_dimension
  public  :: copy, axpy, dot, gemm, gesvd, getrf, symm
!
!&<
!
  integer, parameter           :: IK = KIND(0)
  !! Selected integer kind.
  integer, parameter           :: RK = KIND(0.0_R8)
  !! Selected real kind.
  integer, parameter           :: LK = KIND(.true.)
  !! Selected logical kind.
!
  real(RK), parameter          :: RZERO = 0.0_RK
  !! Real zero.
  real(RK), parameter          :: RONE  = 1.0_RK
  !! Real one.
  real(RK), parameter          :: RHALF = 0.5_RK
  !! Real 1/2.
  real(RK), parameter          :: RFOUR = 4.0_RK
  !! Real four.
  real(RK), parameter          :: RTEN  = 10.0_RK
  !! Real ten.
!
  real(RK), parameter          :: RNAPIER = 2.71828182846_RK
  !! Napier constant.
  real(RK), parameter          :: RHUGE = HUGE(RZERO)
  !! Real large number.
  real(RK), parameter          :: LN_TO_L10 = LOG10(RNAPIER)
  !! Scaling factor, LOG10(x) =LN_TO_L10 * LN(x)
!
  integer(IK), save, protected :: D  = 3
  !! Spatial dimension
  integer(IK), save, protected :: DD = 9
  !! Square spatial dimension
!
!&>
!
contains
!
!| Sets the dimensions of the space. <br>
!  Caution, this routine affects global.
  subroutine setup_dimension(d_)
    integer(IK), intent(in) :: d_
    D = MAX(1, d_)
    DD = D * D
  end subroutine setup_dimension
!
end module mod_params


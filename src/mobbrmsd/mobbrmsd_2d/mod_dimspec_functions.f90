!| Define spatial dimension, \(D=2\),
!  and provide optimized functions for dimension.
module mod_dimspec_functions
  use mod_kinds, only: RK, IK
  implicit none
  private
  public  :: D, DD, ND
  public  :: setup_dimension
  public  :: compute_com
  public  :: compute_cov
  public  :: covcopy
!
  !! Spatial dimension
  integer(IK), parameter :: D = 2
  !! Square spatial dimension
  integer(IK), parameter :: DD = 4
  !| Node memory size, defined by \(1 + 1 + D^2\).
  !  Let \([L, G, \mathbf{C}]\) be a node,
  !  where \(L, G\in\mathbb{R}\) and \(\mathbf{C}\in\mathbb{R}^{D\times D}\).
  integer(IK), parameter :: ND = DD + 2
!
  real(RK), parameter :: ZERO = 0.0_RK
  real(RK), parameter :: HALF = 0.5_RK
  real(RK), parameter :: ONETHIRD = 1.0_RK / 3.0_RK
  real(RK), parameter :: ONE = 1.0_RK
!
  interface
    include 'dgemm.h'
    include 'sgemm.h'
  end interface
!
contains
!| Sets the dimensions of the space. <br>
!  This is dummy interface.
  subroutine setup_dimension(d_)
    integer(IK), intent(in) :: d_
  end subroutine setup_dimension
!
!| Calculate center of mass for \(D=2\)
  pure subroutine compute_com(d, n, X, C)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *)
    real(RK), intent(inout) :: C(d)
    real(RK)                :: rn
    integer(IK)             :: i
    if (n < 0) then
      C(1) = ZERO
      C(2) = ZERO
      return
    elseif (n == 1) then
      C(1) = X(1, 1)
      C(2) = X(2, 1)
      return
    elseif (n == 2) then
      C(1) = (X(1, 1) + X(1, 2)) * HALF
      C(2) = (X(2, 1) + X(2, 2)) * HALF
      return
    elseif (n == 3) then
      C(1) = (X(1, 1) + X(1, 2) + X(1, 3)) * ONETHIRD
      C(2) = (X(2, 1) + X(2, 2) + X(2, 3)) * ONETHIRD
      return
    end if
    C(1) = ZERO
    C(2) = ZERO
    do i = 2, n, 2
      C(1) = C(1) + X(1, i - 1)
      C(2) = C(2) + X(2, i - 1)
      C(1) = C(1) + X(1, i - 0)
      C(2) = C(2) + X(2, i - 0)
    end do
    if (MODULO(n, 2) == 1) then
      C(1) = C(1) + X(1, n)
      C(2) = C(2) + X(2, n)
    end if
    rn = ONE / real(n, RK)
    C(1) = C(1) * rn
    C(2) = C(2) * rn
  end subroutine compute_com
!
!| Calculate covariance matrix.
  pure subroutine compute_cov(d, n, X, Y, C)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), Y(d, *)
    real(RK), intent(inout) :: C(d, d)
#ifdef REAL32
    call SGEMM('N', 'T', D, D, n, ONE, Y, D, X, D, ZERO, C, D)
#else
    call DGEMM('N', 'T', D, D, n, ONE, Y, D, X, D, ZERO, C, D)
#endif
  end subroutine compute_cov
!
!| Compute \(\mathbf{Y} \gets \mathbf{X}-\bar{\mathbf{X}} \)
!  for \(D=2\)
  pure subroutine covcopy(d, n, X, CX, Y)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), CX(d)
    real(RK), intent(inout) :: Y(d, *)
    integer(IK)             :: i
    do concurrent(i=1:n)
      Y(1, i) = X(1, i) - CX(1)
      Y(2, i) = X(2, i) - CX(2)
    end do
  end subroutine covcopy
!
end module mod_dimspec_functions
!
!| DGEMM for \(M=N=2\). <br>
!  \(N\) and \(M\) are provided for compatibility with BLAS and are not used here. <br>
!  @warning
!    This is not a full-featured routine for GEMM. <br>
!    Do not use this routine for anything other than calculating the covariance matrix. <br>
!    Calculate only operations with \( \mathbf{C} = \mathbf{A} \mathbf{B}^\top \).
!  @endwarning
pure subroutine DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  use, intrinsic :: ISO_FORTRAN_ENV, only: RK => REAL64
  character(1), intent(in) :: TRANSA
  !! IF TRANSA(1)='T', transpose A.
  character(1), intent(in) :: TRANSB
  !! IF TRANSB(1)='T', transpose B.
  integer, intent(in)      :: M
  !! matrix dimension, not used.
  integer, intent(in)      :: N
  !! matrix dimension, not used.
  integer, intent(in)      :: K
  !! matrix dimension, not used.
  real(RK), intent(in)     :: ALPHA
  !! A coefficient, not used.
  integer, intent(in)      :: LDA
  !! leading dimension of A, must be >1.
  real(RK), intent(in)     :: A(LDA, *)
  !! matrix A.
  integer, intent(in)      :: LDB
  !! leading dimension of B, must be >1.
  real(RK), intent(in)     :: B(LDA, *)
  !! matrix B.
  real(RK), intent(in)     :: BETA
  !! A coefficient, not used.
  integer, intent(in)      :: LDC
  !! leading dimension of C, must be >1.
  real(RK), intent(inout)  :: C(LDA, *)
  !! matrix C.
  integer                  :: i
!
  C(1, 1) = 0.0_RK
  C(2, 1) = 0.0_RK
  C(1, 2) = 0.0_RK
  C(2, 2) = 0.0_RK
!
  do i = 4, K, 4
    C(1, 1) = C(1, 1) + A(1, i - 3) * B(1, i - 3) &
   &                  + A(1, i - 2) * B(1, i - 2) &
   &                  + A(1, i - 1) * B(1, i - 1) &
   &                  + A(1, i - 0) * B(1, i - 0)
    C(2, 1) = C(2, 1) + A(2, i - 3) * B(1, i - 3) &
   &                  + A(2, i - 2) * B(1, i - 2) &
   &                  + A(2, i - 1) * B(1, i - 1) &
   &                  + A(2, i - 0) * B(1, i - 0)
    C(1, 2) = C(1, 2) + A(1, i - 3) * B(2, i - 3) &
   &                  + A(1, i - 2) * B(2, i - 2) &
   &                  + A(1, i - 1) * B(2, i - 1) &
   &                  + A(1, i - 0) * B(2, i - 0)
    C(2, 2) = C(2, 2) + A(2, i - 3) * B(2, i - 3) &
   &                  + A(2, i - 2) * B(2, i - 2) &
   &                  + A(2, i - 1) * B(2, i - 1) &
   &                  + A(2, i - 0) * B(2, i - 0)
  end do
!
  select case (MODULO(K, 4))
  case (1)
    C(1, 1) = C(1, 1) + A(1, K) * B(1, K)
    C(2, 1) = C(2, 1) + A(2, K) * B(1, K)
    C(1, 2) = C(1, 2) + A(1, K) * B(2, K)
    C(2, 2) = C(2, 2) + A(2, K) * B(2, K)
  case (2)
    C(1, 1) = C(1, 1) + A(1, K - 1) * B(1, K - 1) + A(1, K) * B(1, K)
    C(2, 1) = C(2, 1) + A(2, K - 1) * B(1, K - 1) + A(2, K) * B(1, K)
    C(1, 2) = C(1, 2) + A(1, K - 1) * B(2, K - 1) + A(1, K) * B(2, K)
    C(2, 2) = C(2, 2) + A(2, K - 1) * B(2, K - 1) + A(2, K) * B(2, K)
  case (3)
    C(1, 1) = C(1, 1) + A(1, K - 2) * B(1, K - 2) + A(1, K - 1) * B(1, K - 1) + A(1, K) * B(1, K)
    C(2, 1) = C(2, 1) + A(2, K - 2) * B(1, K - 2) + A(2, K - 1) * B(1, K - 1) + A(2, K) * B(1, K)
    C(1, 2) = C(1, 2) + A(1, K - 2) * B(2, K - 2) + A(1, K - 1) * B(2, K - 1) + A(1, K) * B(2, K)
    C(2, 2) = C(2, 2) + A(2, K - 2) * B(2, K - 2) + A(2, K - 1) * B(2, K - 1) + A(2, K) * B(2, K)
  end select
!
! do i = 2, K, 2
!   C(1, 1) = C(1, 1) + A(1, i - 1) * B(1, i - 1) + A(1, i) * B(1, i)
!   C(2, 1) = C(2, 1) + A(2, i - 1) * B(1, i - 1) + A(2, i) * B(1, i)
! enddo
! do i = 2, K, 2
!   C(1, 2) = C(1, 2) + A(1, i - 1) * B(2, i - 1) + A(1, i) * B(2, i)
!   C(2, 2) = C(2, 2) + A(2, i - 1) * B(2, i - 1) + A(2, i) * B(2, i)
! enddo
! if (MODULO(K, 2) == 1) then
!   C(1, 1) = C(1, 1) + A(1, K) * B(1, K)
!   C(2, 1) = C(2, 1) + A(2, K) * B(1, K)
!   C(1, 2) = C(1, 2) + A(1, K) * B(2, K)
!   C(2, 2) = C(2, 2) + A(2, K) * B(2, K)
! end if
!
! do i = 1, K
!   C(1, 1) = C(1, 1) + A(1, i) * B(1, i)
!   C(2, 1) = C(2, 1) + A(2, i) * B(1, i)
! enddo
! do i = 1, K
!   C(1, 2) = C(1, 2) + A(1, i) * B(2, i)
!   C(2, 2) = C(2, 2) + A(2, i) * B(2, i)
! enddo
!
! do i = 1, K
!   C(1, 1) = C(1, 1) + A(1, i) * B(1, i)
!   C(2, 1) = C(2, 1) + A(2, i) * B(1, i)
!   C(1, 2) = C(1, 2) + A(1, i) * B(2, i)
!   C(2, 2) = C(2, 2) + A(2, i) * B(2, i)
! enddo
!
end subroutine DGEMM

!| SGEMM for \(M=N=2\). <br>
!  \(N\) and \(M\) are provided for compatibility with BLAS and are not used here. <br>
!  @warning
!    This is not a full-featured routine for GEMM. <br>
!    Do not use this routine for anything other than calculating the covariance matrix. <br>
!    Calculate only operations with \( \mathbf{C} = \mathbf{A} \mathbf{B}^\top \).
!  @endwarning
pure subroutine SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  use, intrinsic :: ISO_FORTRAN_ENV, only: RK => REAL32
  character(1), intent(in) :: TRANSA
  !! IF TRANSA(1)='T', transpose A.
  character(1), intent(in) :: TRANSB
  !! IF TRANSB(1)='T', transpose B.
  integer, intent(in)      :: M
  !! matrix dimension, not used.
  integer, intent(in)      :: N
  !! matrix dimension, not used.
  integer, intent(in)      :: K
  !! matrix dimension, not used.
  real(RK), intent(in)     :: ALPHA
  !! A coefficient, not used.
  integer, intent(in)      :: LDA
  !! leading dimension of A, must be >1.
  real(RK), intent(in)     :: A(LDA, *)
  !! matrix A.
  integer, intent(in)      :: LDB
  !! leading dimension of B, must be >1.
  real(RK), intent(in)     :: B(LDB, *)
  !! matrix B.
  real(RK), intent(in)     :: BETA
  !! A coefficient, not used.
  integer, intent(in)      :: LDC
  !! leading dimension of C, must be >1.
  real(RK), intent(inout)  :: C(LDC, *)
  !! matrix C.
  integer                  :: i
!
  C(1, 1) = 0.0_RK
  C(2, 1) = 0.0_RK
  C(1, 2) = 0.0_RK
  C(2, 2) = 0.0_RK
!
  do i = 4, K, 4
    C(1, 1) = C(1, 1) + A(1, i - 3) * B(1, i - 3) + A(1, i - 2) * B(1, i - 2) &
   &                  + A(1, i - 1) * B(1, i - 1) + A(1, i - 0) * B(1, i - 0)
    C(2, 1) = C(2, 1) + A(2, i - 3) * B(1, i - 3) + A(2, i - 2) * B(1, i - 2) &
   &                  + A(2, i - 1) * B(1, i - 1) + A(2, i - 0) * B(1, i - 0)
    C(1, 2) = C(1, 2) + A(1, i - 3) * B(2, i - 3) + A(1, i - 2) * B(2, i - 2) &
   &                  + A(1, i - 1) * B(2, i - 1) + A(1, i - 0) * B(2, i - 0)
    C(2, 2) = C(2, 2) + A(2, i - 3) * B(2, i - 3) + A(2, i - 2) * B(2, i - 2) &
   &                  + A(2, i - 1) * B(2, i - 1) + A(2, i - 0) * B(2, i - 0)
  end do
!
  select case (MODULO(K, 4))
  case (1)
    C(1, 1) = C(1, 1) + A(1, K) * B(1, K)
    C(2, 1) = C(2, 1) + A(2, K) * B(1, K)
    C(1, 2) = C(1, 2) + A(1, K) * B(2, K)
    C(2, 2) = C(2, 2) + A(2, K) * B(2, K)
  case (2)
    C(1, 1) = C(1, 1) + A(1, K - 1) * B(1, K - 1) + A(1, K) * B(1, K)
    C(2, 1) = C(2, 1) + A(2, K - 1) * B(1, K - 1) + A(2, K) * B(1, K)
    C(1, 2) = C(1, 2) + A(1, K - 1) * B(2, K - 1) + A(1, K) * B(2, K)
    C(2, 2) = C(2, 2) + A(2, K - 1) * B(2, K - 1) + A(2, K) * B(2, K)
  case (3)
    C(1, 1) = C(1, 1) + A(1, K - 2) * B(1, K - 2) + A(1, K - 1) * B(1, K - 1) + A(1, K) * B(1, K)
    C(2, 1) = C(2, 1) + A(2, K - 2) * B(1, K - 2) + A(2, K - 1) * B(1, K - 1) + A(2, K) * B(1, K)
    C(1, 2) = C(1, 2) + A(1, K - 2) * B(2, K - 2) + A(1, K - 1) * B(2, K - 1) + A(1, K) * B(2, K)
    C(2, 2) = C(2, 2) + A(2, K - 2) * B(2, K - 2) + A(2, K - 1) * B(2, K - 1) + A(2, K) * B(2, K)
  end select
!
end subroutine SGEMM

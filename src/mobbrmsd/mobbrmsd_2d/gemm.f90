!| DGEMM for M=N=2. <br>
!  N and M are provided for compatibility with BLAS and are not used here. <br>
!  Caution! This is not a full-featured routine for GEMM. <br>
!  Do not use this routine for anything other than calculating the covariance matrix. <br>
!  Calculate only operations with C = A@Transpose(B).
pure subroutine DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  !$ use omp_lib
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
  enddo
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

!| SGEMM for M=N=2. <br>
!  N and M are provided for compatibility with BLAS and are not used here. <br>
!  Caution! This is not a full-featured routine for GEMM. <br>
!  Do not use this routine for anything other than calculating the covariance matrix. <br>
!  Calculate only operations with C = A@Transpose(B).
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
  enddo
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


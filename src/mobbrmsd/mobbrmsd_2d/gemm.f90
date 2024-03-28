!| SGEMM for d=2. (N, M, and K are provided for compatibility with BLAS and are not used here.) <br>
!  Caution! This is not a full-featured routine for GEMM. <br>
!  Calculate only operations with C = A@B or its transposes.
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
  real(RK), intent(in)     :: B(LDB, *)
  !! matrix B.
  real(RK), intent(in)     :: BETA
  !! A coefficient, not used.
  integer, intent(in)      :: LDC
  !! leading dimension of C, must be >1.
  real(RK), intent(inout)  :: C(LDC, *)
  !! matrix C.
  if (TRANSA=='T') then
    if (TRANSB=='T') then
      C(1, 1) = A(1, 1) * B(1, 1) + A(2, 1) * B(1, 2)
      C(2, 1) = A(1, 2) * B(1, 1) + A(2, 2) * B(1, 2)
      C(1, 2) = A(1, 1) * B(2, 1) + A(2, 1) * B(2, 2)
      C(2, 2) = A(1, 2) * B(2, 1) + A(2, 2) * B(2, 2)
    else
      C(1, 1) = A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1)
      C(2, 1) = A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1)
      C(1, 2) = A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2)
      C(2, 2) = A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2)
    end if
  else
    if (TRANSB=='T') then
      C(1, 1) = A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2)
      C(2, 1) = A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2)
      C(1, 2) = A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2)
      C(2, 2) = A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2)
    else
      C(1, 1) = A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1)
      C(2, 1) = A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1)
      C(1, 2) = A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2)
      C(2, 2) = A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2)
    end if
  end if
end subroutine DGEMM

!| SGEMM for d=2. (N, M, and K are provided for compatibility with BLAS and are not used here.) <br>
!  Caution! This is not a full-featured routine for GEMM. <br>
!  Calculate only operations with C = A@B or its transposes.
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
  if (TRANSA=='T') then
    if (TRANSB=='T') then
      C(1, 1) = A(1, 1) * B(1, 1) + A(2, 1) * B(1, 2)
      C(2, 1) = A(1, 2) * B(1, 1) + A(2, 2) * B(1, 2)
      C(1, 2) = A(1, 1) * B(2, 1) + A(2, 1) * B(2, 2)
      C(2, 2) = A(1, 2) * B(2, 1) + A(2, 2) * B(2, 2)
    else
      C(1, 1) = A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1)
      C(2, 1) = A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1)
      C(1, 2) = A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2)
      C(2, 2) = A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2)
    end if
  else
    if (TRANSB=='T') then
      C(1, 1) = A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2)
      C(2, 1) = A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2)
      C(1, 2) = A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2)
      C(2, 2) = A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2)
    else
      C(1, 1) = A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1)
      C(2, 1) = A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1)
      C(1, 2) = A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2)
      C(2, 2) = A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2)
    end if
  end if
end subroutine SGEMM


program main
  use mod_mobbrmsd_lapack
  use mod_unittest
  implicit none
  type(unittest) :: z
#ifdef USE_REAL32
  integer, parameter :: place = 3
  integer, parameter :: RK = SELECTED_REAL_KIND(6)
#else
  integer, parameter :: place = 6
  integer, parameter :: RK = SELECTED_REAL_KIND(15)
#endif
!
  call z%init('test mobbrmsd_lapack')
!
  call test1(2, 2, 2)
  call test1(3, 3, 3)
  call test1(4, 4, 4)
  call test1(2, 3, 4)
  call test1(4, 3, 2)
  call test1(4, 2, 3)
  call test1(10, 100, 1000)
  call test1(4, 12, 40)
! call test1(50, 200, 35)
! call test1(441, 1230, 340)
!
  call test2(2)
  call test2(3)
  call test2(4)
!
  call test3(2)
  call test3(3)
  call test3(4)
!
  call test4(2)
  call test4(3)
  call test4(4)
!
  call z%finish_and_terminate()
!
contains
!
!| Test GEMM
  subroutine test1(M, N, K)
    integer, intent(in) :: M, N, K
    real(RK) :: A(M, K), B(K, N), C(K, M), D(N, K)
    real(RK) :: AB(M, N), CB(M, N), AD(M, N), CD(M, N)
    real(RK) :: BA(N, M), BC(N, M), DA(N, M), DC(N, M)
    real(RK) :: AB_(M, N), CB_(M, N), AD_(M, N), CD_(M, N)
    real(RK) :: BA_(N, M), BC_(N, M), DA_(N, M), DC_(N, M)
!
    call RANDOM_NUMBER(A)
    call RANDOM_NUMBER(B)
    C = TRANSPOSE(A)
    D = TRANSPOSE(B)
#ifdef USE_REAL32
    call SGEMM("N", "N", M, N, K, 1.0_RK, A, M, B, K, 0.0_RK, AB, M)
    call SGEMM("T", "N", M, N, K, 1.0_RK, C, K, B, K, 0.0_RK, CB, M)
    call SGEMM("N", "T", M, N, K, 1.0_RK, A, M, D, N, 0.0_RK, AD, M)
    call SGEMM("T", "T", M, N, K, 1.0_RK, C, K, D, N, 0.0_RK, CD, M)
    call SGEMM("T", "T", N, M, K, 1.0_RK, B, K, A, M, 0.0_RK, BA, N)
    call SGEMM("T", "N", N, M, K, 1.0_RK, B, K, C, K, 0.0_RK, BC, N)
    call SGEMM("N", "T", N, M, K, 1.0_RK, D, N, A, M, 0.0_RK, DA, N)
    call SGEMM("N", "N", N, M, K, 1.0_RK, D, N, C, K, 0.0_RK, DC, N)
#else
    call DGEMM("N", "N", M, N, K, 1.0_RK, A, M, B, K, 0.0_RK, AB, M)
    call DGEMM("T", "N", M, N, K, 1.0_RK, C, K, B, K, 0.0_RK, CB, M)
    call DGEMM("N", "T", M, N, K, 1.0_RK, A, M, D, N, 0.0_RK, AD, M)
    call DGEMM("T", "T", M, N, K, 1.0_RK, C, K, D, N, 0.0_RK, CD, M)
    call DGEMM("T", "T", N, M, K, 1.0_RK, B, K, A, M, 0.0_RK, BA, N)
    call DGEMM("T", "N", N, M, K, 1.0_RK, B, K, C, K, 0.0_RK, BC, N)
    call DGEMM("N", "T", N, M, K, 1.0_RK, D, N, A, M, 0.0_RK, DA, N)
    call DGEMM("N", "N", N, M, K, 1.0_RK, D, N, C, K, 0.0_RK, DC, N)
#endif
    AB_ = MATMUL(A, B)
    CB_ = MATMUL(TRANSPOSE(C), B)
    AD_ = MATMUL(A, TRANSPOSE(D))
    CD_ = MATMUL(TRANSPOSE(C), TRANSPOSE(D))
    BA_ = MATMUL(TRANSPOSE(B), TRANSPOSE(A))
    BC_ = MATMUL(TRANSPOSE(B), C)
    DA_ = MATMUL(D, TRANSPOSE(A))
    DC_ = MATMUL(D, C)
    call z%assert_almost_equal(AB, AB_, 'gemm : A @ B', place=place)
    call z%assert_almost_equal(CB, CB_, 'gemm : C^T @ B', place=place)
    call z%assert_almost_equal(AD, AD_, 'gemm : A @ D^T', place=place)
    call z%assert_almost_equal(CD, CD_, 'gemm : C^T @ D^T', place=place)
    call z%assert_almost_equal(BA, BA_, 'gemm : B^T @ A^T', place=place)
    call z%assert_almost_equal(BC, BC_, 'gemm : B^T @ C', place=place)
    call z%assert_almost_equal(DA, DA_, 'gemm : D @ A^T', place=place)
    call z%assert_almost_equal(DC, DC_, 'gemm : D @ C', place=place)
    call z%assert_almost_equal(AB, CB, 'gemm : AB = CB', place=place)
    call z%assert_almost_equal(AB, AD, 'gemm : AB = AD', place=place)
    call z%assert_almost_equal(AB, CD, 'gemm : AB = CD', place=place)
    call z%assert_almost_equal(BA, BC, 'gemm : BA = BC', place=place)
    call z%assert_almost_equal(BA, DA, 'gemm : BA = DA', place=place)
    call z%assert_almost_equal(BA, DC, 'gemm : BA = DC', place=place)
!
    CD = AB
    AD = AB
    CD = AB
    CD_ = AB
    AD_ = AB
    CD_ = AB
    BA = BC
    BA = DA
    BA = DC
    BA_ = BC
    BA_ = DA
    BA_ = DC
!
#ifdef USE_REAL32
    call SGEMM("N", "N", M, N, K, 1.0_RK, A, M, B, K, 1.0_RK, AB, M)
    call SGEMM("T", "N", M, N, K, 1.0_RK, C, K, B, K, 1.0_RK, CB, M)
    call SGEMM("N", "T", M, N, K, 1.0_RK, A, M, D, N, 1.0_RK, AD, M)
    call SGEMM("T", "T", M, N, K, 1.0_RK, C, K, D, N, 1.0_RK, CD, M)
    call SGEMM("T", "T", N, M, K, 1.0_RK, B, K, A, M, 1.0_RK, BA, N)
    call SGEMM("T", "N", N, M, K, 1.0_RK, B, K, C, K, 1.0_RK, BC, N)
    call SGEMM("N", "T", N, M, K, 1.0_RK, D, N, A, M, 1.0_RK, DA, N)
    call SGEMM("N", "N", N, M, K, 1.0_RK, D, N, C, K, 1.0_RK, DC, N)
#else
    call DGEMM("N", "N", M, N, K, 1.0_RK, A, M, B, K, 1.0_RK, AB, M)
    call DGEMM("T", "N", M, N, K, 1.0_RK, C, K, B, K, 1.0_RK, CB, M)
    call DGEMM("N", "T", M, N, K, 1.0_RK, A, M, D, N, 1.0_RK, AD, M)
    call DGEMM("T", "T", M, N, K, 1.0_RK, C, K, D, N, 1.0_RK, CD, M)
    call DGEMM("T", "T", N, M, K, 1.0_RK, B, K, A, M, 1.0_RK, BA, N)
    call DGEMM("T", "N", N, M, K, 1.0_RK, B, K, C, K, 1.0_RK, BC, N)
    call DGEMM("N", "T", N, M, K, 1.0_RK, D, N, A, M, 1.0_RK, DA, N)
    call DGEMM("N", "N", N, M, K, 1.0_RK, D, N, C, K, 1.0_RK, DC, N)
#endif
    AB_ = AB_ + MATMUL(A, B)
    CB_ = CB_ + MATMUL(TRANSPOSE(C), B)
    AD_ = AD_ + MATMUL(A, TRANSPOSE(D))
    CD_ = CD_ + MATMUL(TRANSPOSE(C), TRANSPOSE(D))
    BA_ = BA_ + MATMUL(TRANSPOSE(B), TRANSPOSE(A))
    BC_ = BC_ + MATMUL(TRANSPOSE(B), C)
    DA_ = DA_ + MATMUL(D, TRANSPOSE(A))
    DC_ = DC_ + MATMUL(D, C)
    call z%assert_almost_equal(AB, AB_, 'gemm : A @ B + E', place=place)
    call z%assert_almost_equal(CB, CB_, 'gemm : C^T @ B + E', place=place)
    call z%assert_almost_equal(AD, AD_, 'gemm : A @ D^T + E', place=place)
    call z%assert_almost_equal(CD, CD_, 'gemm : C^T @ D^T + E', place=place)
    call z%assert_almost_equal(BA, BA_, 'gemm : B^T @ A^T + E', place=place)
    call z%assert_almost_equal(BC, BC_, 'gemm : B^T @ C + E', place=place)
    call z%assert_almost_equal(DA, DA_, 'gemm : D @ A^T + E', place=place)
    call z%assert_almost_equal(DC, DC_, 'gemm : D @ C + E', place=place)
    call z%assert_almost_equal(AB, CB, 'gemm : AB = CB', place=place)
    call z%assert_almost_equal(AB, AD, 'gemm : AB = AD', place=place)
    call z%assert_almost_equal(AB, CD, 'gemm : AB = CD', place=place)
    call z%assert_almost_equal(BA, BC, 'gemm : BA = BC', place=place)
    call z%assert_almost_equal(BA, DA, 'gemm : BA = DA', place=place)
    call z%assert_almost_equal(BA, DC, 'gemm : BA = DC', place=place)
!
    CD = AB
    AD = AB
    CD = AB
    CD_ = AB
    AD_ = AB
    CD_ = AB
    BA = BC
    BA = DA
    BA = DC
    BA_ = BC
    BA_ = DA
    BA_ = DC
!
#ifdef USE_REAL32
    call SGEMM("N", "N", M, N, K, 0.5_RK, A, M, B, K, 1.0_RK, AB, M)
    call SGEMM("T", "N", M, N, K, 0.5_RK, C, K, B, K, 1.0_RK, CB, M)
    call SGEMM("N", "T", M, N, K, 0.5_RK, A, M, D, N, 1.0_RK, AD, M)
    call SGEMM("T", "T", M, N, K, 0.5_RK, C, K, D, N, 1.0_RK, CD, M)
    call SGEMM("T", "T", N, M, K, 0.5_RK, B, K, A, M, 1.0_RK, BA, N)
    call SGEMM("T", "N", N, M, K, 0.5_RK, B, K, C, K, 1.0_RK, BC, N)
    call SGEMM("N", "T", N, M, K, 0.5_RK, D, N, A, M, 1.0_RK, DA, N)
    call SGEMM("N", "N", N, M, K, 0.5_RK, D, N, C, K, 1.0_RK, DC, N)
#else
    call DGEMM("N", "N", M, N, K, 0.5_RK, A, M, B, K, 1.0_RK, AB, M)
    call DGEMM("T", "N", M, N, K, 0.5_RK, C, K, B, K, 1.0_RK, CB, M)
    call DGEMM("N", "T", M, N, K, 0.5_RK, A, M, D, N, 1.0_RK, AD, M)
    call DGEMM("T", "T", M, N, K, 0.5_RK, C, K, D, N, 1.0_RK, CD, M)
    call DGEMM("T", "T", N, M, K, 0.5_RK, B, K, A, M, 1.0_RK, BA, N)
    call DGEMM("T", "N", N, M, K, 0.5_RK, B, K, C, K, 1.0_RK, BC, N)
    call DGEMM("N", "T", N, M, K, 0.5_RK, D, N, A, M, 1.0_RK, DA, N)
    call DGEMM("N", "N", N, M, K, 0.5_RK, D, N, C, K, 1.0_RK, DC, N)
#endif
    AB_ = AB_ + 0.5_RK * MATMUL(A, B)
    CB_ = CB_ + 0.5_RK * MATMUL(TRANSPOSE(C), B)
    AD_ = AD_ + 0.5_RK * MATMUL(A, TRANSPOSE(D))
    CD_ = CD_ + 0.5_RK * MATMUL(TRANSPOSE(C), TRANSPOSE(D))
    BA_ = BA_ + 0.5_RK * MATMUL(TRANSPOSE(B), TRANSPOSE(A))
    BC_ = BC_ + 0.5_RK * MATMUL(TRANSPOSE(B), C)
    DA_ = DA_ + 0.5_RK * MATMUL(D, TRANSPOSE(A))
    DC_ = DC_ + 0.5_RK * MATMUL(D, C)
    call z%assert_almost_equal(AB, AB_, 'gemm : 0.5 A @ B + E', place=place)
    call z%assert_almost_equal(CB, CB_, 'gemm : 0.5 C^T @ B + E', place=place)
    call z%assert_almost_equal(AD, AD_, 'gemm : 0.5 A @ D^T + E', place=place)
    call z%assert_almost_equal(CD, CD_, 'gemm : 0.5 C^T @ D^T + E', place=place)
    call z%assert_almost_equal(BA, BA_, 'gemm : 0.5 B^T @ A^T + E', place=place)
    call z%assert_almost_equal(BC, BC_, 'gemm : 0.5 B^T @ C + E', place=place)
    call z%assert_almost_equal(DA, DA_, 'gemm : 0.5 D @ A^T + E', place=place)
    call z%assert_almost_equal(DC, DC_, 'gemm : 0.5 D @ C + E', place=place)
    call z%assert_almost_equal(AB, CB, 'gemm : AB = CB', place=place)
    call z%assert_almost_equal(AB, AD, 'gemm : AB = AD', place=place)
    call z%assert_almost_equal(AB, CD, 'gemm : AB = CD', place=place)
    call z%assert_almost_equal(BA, BC, 'gemm : BA = BC', place=place)
    call z%assert_almost_equal(BA, DA, 'gemm : BA = DA', place=place)
    call z%assert_almost_equal(BA, DC, 'gemm : BA = DC', place=place)
  end subroutine test1
!
!| Test GEQRF and ORMQR
  subroutine test2(D)
    integer, intent(in) :: D
    real(RK)            :: X(D * D + 1), A(D * D), SO(D, D), tau(D), w(D * 3)
    integer             :: i, j, info
    call RANDOM_NUMBER(X)
    do i = 1, D * D
      A(i) = SQRT(-2._RK * LOG(X(i))) * COS(6.2831853070_RK * X(i + 1))
    end do
    do concurrent(i=1:D, j=1:D)
      SO(i, j) = MERGE(1.0_RK, 0.0_RK, i == j)
    end do
#ifdef USE_REAL32
    call SGEQRF(D, D, A, D, tau, w, SIZE(w), info)
    call SORMQR("L", "N", D, D, D, A, D, tau, SO, D, w, SIZE(w), info)
#else
    call DGEQRF(D, D, A, D, tau, w, SIZE(w), info)
    call DORMQR("L", "N", D, D, D, A, D, tau, SO, D, w, SIZE(w), info)
#endif
    call z%assert_is_eye(MATMUL(SO, TRANSPOSE(SO)), 'geqrf and ormqr : SO @ SO^T = I', place=place)
    call z%assert_is_eye(MATMUL(TRANSPOSE(SO), SO), 'geqrf and ormqr : SO^T @ SO = I', place=place)
  end subroutine test2
!
!| Test SGETRF
  subroutine test3(D)
    integer, intent(in)  :: D
    real(RK)             :: X(D, D), Y(D, D), det, ref
    real(RK)             :: L(D, D), U(D, D), LU(D, D), P(D, D), S(D)
    integer              :: i, j, info, ipiv(D)
!
    call RANDOM_NUMBER(X)
    Y = X
    if (D == 2) then
      ref = X(1, 1) * X(2, 2) - X(1, 2) * X(2, 1)
    else if (D == 3) then
      ref = X(1, 1) * (X(2, 2) * X(3, 3) - X(2, 3) * X(3, 2)) &
     &    - X(1, 2) * (X(2, 1) * X(3, 3) - X(2, 3) * X(3, 1)) &
     &    + X(1, 3) * (X(2, 1) * X(3, 2) - X(2, 2) * X(3, 1))
    else if (D == 4) then
      ref = X(1, 1) * (X(2, 2) * (X(3, 3) * X(4, 4) - X(3, 4) * X(4, 3)) &
     &               - X(2, 3) * (X(3, 2) * X(4, 4) - X(3, 4) * X(4, 2)) &
     &               + X(2, 4) * (X(3, 2) * X(4, 3) - X(3, 3) * X(4, 2))) &
     &    - X(2, 1) * (X(1, 2) * (X(3, 3) * X(4, 4) - X(3, 4) * X(4, 3)) &
     &               - X(1, 3) * (X(3, 2) * X(4, 4) - X(3, 4) * X(4, 2)) &
     &               + X(1, 4) * (X(3, 2) * X(4, 3) - X(3, 3) * X(4, 2))) &
     &    + X(3, 1) * (X(1, 2) * (X(2, 3) * X(4, 4) - X(2, 4) * X(4, 3)) &
     &               - X(1, 3) * (X(2, 2) * X(4, 4) - X(2, 4) * X(4, 2)) &
     &               + X(1, 4) * (X(2, 2) * X(4, 3) - X(2, 3) * X(4, 2))) &
     &    - X(4, 1) * (X(1, 2) * (X(2, 3) * X(3, 4) - X(2, 4) * X(3, 3)) &
     &               - X(1, 3) * (X(2, 2) * X(3, 4) - X(2, 4) * X(3, 2)) &
     &               + X(1, 4) * (X(2, 2) * X(3, 3) - X(2, 3) * X(3, 2)))
    else
      ref = 0.0_RK
    end if
!
#ifdef USE_REAL32
    call SGETRF(D, D, X, D, ipiv, info)
#else
    call DGETRF(D, D, X, D, ipiv, info)
#endif
!
    U = 0.0_RK
    do concurrent(i=1:D, j=1:D)
      P(i, j) = MERGE(1.0_RK, 0.0_RK, i == j)
      L(i, j) = MERGE(1.0_RK, 0.0_RK, i == j)
    end do
    do i = 1, D
      S = P(:, i)
      P(:, i) = P(:, ipiv(i))
      P(:, ipiv(i)) = S
    end do
    do j = 1, D
      do i = 1, j
        U(i, j) = X(i, j)
      end do
      do i = j + 1, D
        L(i, j) = X(i, j)
      end do
    end do
    LU = MATMUL(P, MATMUL(L, U))
!
    det = 1.0_RK
    do i = 1, D
      det = det * X(i, i)
      if (ipiv(i) /= i) det = -det
    end do
!
    call z%assert_almost_equal(Y, LU, 'getrf : X = L @ U', place=place)
    call z%assert_almost_equal(det, ref, 'getrf : |X| = |U|', place=place)
!
  end subroutine test3
!
!| Test SGESVD
  subroutine test4(D)
    integer, intent(in)   :: D
    real(RK)              :: X(D, D), Y(D, D), S(D), U(D, D), VT(D, D), UVT(D, D), Q(1)
    integer               :: i, info
    call RANDOM_NUMBER(X)
    Y = X
#ifdef USE_REAL32
    call SGESVD('A', 'A', D, D, X, D, S, U, D, VT, D, Q, -1, info)
    block
      real(RK) :: W(NINT(Q(1)))
      call SGESVD('A', 'A', D, D, X, D, S, U, D, VT, D, W, SIZE(W), info)
    end block
#else
    call DGESVD('A', 'A', D, D, X, D, S, U, D, VT, D, Q, -1, info)
    block
      real(RK) :: W(NINT(Q(1)))
      call DGESVD('A', 'A', D, D, X, D, S, U, D, VT, D, W, SIZE(W), info)
    end block
#endif
    X = 0.0_RK
    do i = 1, D
      X(i, i) = S(i)
    end do
    UVT = MATMUL(U, VT)
    call z%assert_almost_equal(MATMUL(MATMUL(U, X), VT), Y, 'gesvd : X = U @ S @ V^T', place=place)
    call z%assert_is_eye(MATMUL(UVT, TRANSPOSE(UVT)), 'gesvd : U @ V^T @ V @ U^T = I', place=place)
    call z%assert_is_eye(MATMUL(TRANSPOSE(UVT), UVT), 'gesvd : V @ U^T @ U @ V^T = I', place=place)
  end subroutine test4
!
end program main


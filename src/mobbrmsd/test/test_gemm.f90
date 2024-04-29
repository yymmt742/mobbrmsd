program main
  use mod_dimspec_functions, only: D, DD
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer        :: i
!
  interface
    include 'dgemm.h'
    include 'sgemm.h'
  end interface
!
  call u%init('test gemm')
  do i = 1, 10
    call test1(500 * (i - 1) + 1)
  end do
!
  call u%init('test gemm time')
  print'(A)', '       K     mygemm (ms)     matmul (ms)        error'
  call test2(1)
  call test2(2)
  call test2(5)
  call test2(7)
  call test2(32)
  call test2(49)
  call test2(64)
  call test2(75)
  call test2(128)
  call test2(151)
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(K)
    integer, intent(in) :: K
    real(RK) :: A(D, K), B(D, K), C1(D, D), C2(D, D)
!
    call RANDOM_NUMBER(A)
    call RANDOM_NUMBER(B)
!
#ifdef REAL32
    call sgemm('N', 'T', D, D, K, ONE, A, D, B, D, ZERO, C1, D)
#else
    call dgemm('N', 'T', D, D, K, ONE, A, D, B, D, ZERO, C1, D)
#endif
    C2 = MATMUL(A, TRANSPOSE(B))
    call u%assert_almost_equal([C1 - C2], ZERO, ' A @ TB = C')
!
  end subroutine test1
!
  subroutine test2(K)
    integer(IK), intent(in) :: K
    integer(IK), parameter  :: N = 5000000 / DD
    real(RK)                :: A0(D, K), A(D, K), B(D, K), C1(D, D), C2(D, D)
    real(RK)                :: time_begin_s, time_end_s, time_gemm, time_matmul
    integer(IK)             :: i
!
    call RANDOM_NUMBER(A)
    call RANDOM_NUMBER(B)
!
    A0 = A
!
    call CPU_TIME(time_begin_s)
    do i = 1, N
#ifdef REAL32
      call sgemm('N', 'T', D, D, K, ONE, A, D, B, D, ZERO, C1, D)
#else
      call dgemm('N', 'T', D, D, K, ONE, A, D, B, D, ZERO, C1, D)
#endif
      A = SUM(C1) / N + B
    end do
    call CPU_TIME(time_end_s)
!
    time_gemm = 1000 * (time_end_s - time_begin_s)
!
    A = A0
!
    call CPU_TIME(time_begin_s)
    do i = 1, N
      C2 = MATMUL(A, TRANSPOSE(B))
      A = SUM(C2) / N + B
    end do
    call CPU_TIME(time_end_s)
!
    time_matmul = 1000 * (time_end_s - time_begin_s)
    print'(I8,2f16.3,F16.9)', K, time_gemm, time_matmul, SUM(C1 - C2)
!
  end subroutine test2
!
end program main

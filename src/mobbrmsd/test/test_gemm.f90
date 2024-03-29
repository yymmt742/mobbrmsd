program main
  use mod_params, only: D, DD, RK, IK, gemm, ONE => RONE, ZERO => RZERO
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer(IK)    :: i
!
  call u%init('test gemm')
  do i = 1, 10
    call test1()
  end do
!
  call u%init('test gemm time')
  call test2()
  call test3()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    real(RK) :: A(D, D), B(D, D), C1(D, D), C2(D, D)
!
    call RANDOM_NUMBER(A)
    call RANDOM_NUMBER(B)
    call gemm('A', 'A', D, D, D, ONE, A, D, B, D, ZERO, C1, D)
    C2 = MATMUL(A, B)
    call u%assert_almost_equal([C1 - C2], ZERO, ' A @  B = C')
!
    call gemm('T', 'A', D, D, D, ONE, A, D, B, D, ZERO, C1, D)
    C2 = MATMUL(TRANSPOSE(A), B)
    call u%assert_almost_equal([C1 - C2], ZERO, 'TA @  B = C')
!
    call gemm('A', 'T', D, D, D, ONE, A, D, B, D, ZERO, C1, D)
    C2 = MATMUL(A, TRANSPOSE(B))
    call u%assert_almost_equal([C1 - C2], ZERO, ' A @ TB = C')
    return
!
    call gemm('T', 'T', D, D, D, ONE, A, D, B, D, ZERO, C1, D)
    C2 = MATMUL(TRANSPOSE(A), TRANSPOSE(B))
    call u%assert_almost_equal([C1 - C2], ZERO, 'TA @ TB = C')
!
  end subroutine test1
!
  subroutine test2()
    integer(IK), parameter :: N = 1000000 / DD
    real(RK)    :: A0(D, D), A(D, D), B(D, D), C1(D, D), C2(D, D)
    real(RK)    :: time_begin_s, time_end_s, time_gemm, time_matmul
    integer(IK) :: i
!
    call RANDOM_NUMBER(A)
    call RANDOM_NUMBER(B)
!
    A0 = A
!
    call CPU_TIME(time_begin_s)
    do i = 1, N
      call gemm('A', 'A', D, D, D, ONE, A, D, B, D, ZERO, C1, D)
      A = C1 / N + B
    end do
    call CPU_TIME(time_end_s)
!
    time_gemm = 1000 * (time_end_s - time_begin_s)
!
    A = A0
!
    call CPU_TIME(time_begin_s)
    do i = 1, N
      C2 = MATMUL(A, B)
      A = C2 / N + B
    end do
    call CPU_TIME(time_end_s)
!
    time_matmul = 1000 * (time_end_s - time_begin_s)
    print'(3f16.9)', time_gemm, time_matmul, SUM(C1 - C2)
!
  end subroutine test2
!
  subroutine test3()
    integer(IK), parameter :: N = 1000000 / DD
    real(RK)    :: A0(D, D), A(D, D), B(D, D), C1(D, D), C2(D, D)
    real(RK)    :: time_begin_s, time_end_s, time_gemm, time_matmul
    integer(IK) :: i
!
    call RANDOM_NUMBER(A)
    call RANDOM_NUMBER(B)
!
    A0 = A
!
    call CPU_TIME(time_begin_s)
    do i = 1, N
      call gemm('T', 'A', D, D, D, ONE, A, D, B, D, ZERO, C1, D)
      A = C1 / N + B
    end do
    call CPU_TIME(time_end_s)
!
    time_gemm = 1000 * (time_end_s - time_begin_s)
!
    A = A0
!
    call CPU_TIME(time_begin_s)
    do i = 1, N
      C2 = MATMUL(TRANSPOSE(A), B)
      A = C2 / N + B
    end do
    call CPU_TIME(time_end_s)
!
    time_matmul = 1000 * (time_end_s - time_begin_s)
    print'(3f16.9)', time_gemm, time_matmul, SUM(C1 - C2)
!
  end subroutine test3
!
end program main

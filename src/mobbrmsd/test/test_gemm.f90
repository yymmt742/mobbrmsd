program main
  use mod_params, only: D, RK, IK, gemm, ONE => RONE, ZERO => RZERO
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer(IK)    :: i
!
  call u%init('test det d=1')
  do i = 1, 10
    call test1()
  end do
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    real(RK) :: A(D, D), B(D, D), C(D, D)
!
    call random_number(A)
    call random_number(B)
    call gemm('A', 'A', D, D, D, ONE, A, D, B, D, ZERO, C, D)
    call u%assert_almost_equal([C - MATMUL(A, B)], ZERO, ' A @  B = C')
!
    call gemm('T', 'A', D, D, D, ONE, A, D, B, D, ZERO, C, D)
    call u%assert_almost_equal([C - MATMUL(TRANSPOSE(A), B)], ZERO, 'TA @  B = C')
!
    call gemm('A', 'T', D, D, D, ONE, A, D, B, D, ZERO, C, D)
    call u%assert_almost_equal([C - MATMUL(A, TRANSPOSE(B))], ZERO, ' A @ TB = C')
!
    call gemm('T', 'T', D, D, D, ONE, A, D, B, D, ZERO, C, D)
    call u%assert_almost_equal([C - MATMUL(TRANSPOSE(A), TRANSPOSE(B))], ZERO, 'TA @ TB = C')
!
  end subroutine test1
!
end program main

program main
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use mod_dimspec_functions, only: D, compute_com
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test compute_com')
  call test1(1)
  call test1(2)
  call test1(3)
  call test1(4)
  call test1(5)
  call test1(7)
  call test1(32)
  call test1(49)
  call test1(64)
  call test1(75)
  call test1(128)
  call test1(1051)
!
  call u%init('test compute_com time')
  print'(A)', '       K  mycom (ms)    SUM (ms)'
  call test2(1)
  call test2(2)
  call test2(3)
  call test2(4)
  call test2(5)
  call test2(7)
  call test2(32)
  call test2(49)
  call test2(64)
  call test2(75)
  call test2(128)
  call test2(1051)
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(K)
    integer, intent(in) :: K
    real(RK) :: X(D, K), CX(D)
!
    call RANDOM_NUMBER(X)
    call compute_com(D, K, X, CX)
    call u%assert_almost_equal(CX, SUM(X, 2) / K, 'compute_com')
!
  end subroutine test1
!
  subroutine test2(K)
    integer(IK), intent(in) :: K
    integer(IK)             :: N
    real(RK)                :: X(D, K), CX(D)
    real(RK)                :: res
    real(RK)                :: time_begin_s, time_end_s, time_covdot, time_ref
    integer(IK)             :: i
!
    call RANDOM_NUMBER(X)
    CX = SUM(X, 2) / K
    res = 0.0_RK
    N = 1000000 / (INT(LOG10(real(K, RK))) + 1)
!
    call CPU_TIME(time_begin_s)
    do i = 1, N
      call compute_com(D, K, X, CX)
      res = res + SUM(CX)
      call RANDOM_NUMBER(X(:, 1))
    end do
    call CPU_TIME(time_end_s)
!
    time_covdot = 1000 * (time_end_s - time_begin_s)
!
    call CPU_TIME(time_begin_s)
    do i = 1, N
      CX = SUM(X, 2) / K
      res = res + SUM(CX)
      call RANDOM_NUMBER(X(:, 1))
    end do
    call CPU_TIME(time_end_s)
!
    time_ref = 1000 * (time_end_s - time_begin_s)
    print'(I8,3f12.3)', K, time_covdot, time_ref, time_covdot / time_ref
    FLUSH (OUTPUT_UNIT)
!
  end subroutine test2
!
end program main


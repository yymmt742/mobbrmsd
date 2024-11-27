program main
  use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
  use mod_kinds
  use mod_params, only: RK, IK
  use mod_unittest
  implicit none
  type(unittest) :: u
#ifdef USE_REAL32
  integer, parameter    :: place = 4
#else
  integer, parameter    :: place = 8
#endif
!
  call u%init('test invsqrt')
  call test1()
  call u%init('test invsqrt time')
  call test2()
  call test3()
  call u%finish_and_terminate()
!
contains
!
#include "invsqrt.f90"
!
  subroutine test1()
    real(RK) :: x, y, e1
    integer  :: i
    do i = -100, 100
      x = EXP(i * 0.1)
      y = invsqrt(x)
      e1 = (x - 1 / y**2) / x
      call u%assert_is_zero(e1, 'x = f(x)**(-2)', place=place)
    end do
  end subroutine test1
!
  subroutine test2()
    integer(IK), parameter :: N = 5000000, M = 10
    real(RK)               :: A(N), B(N), C(N)
    real(RK)               :: time_begin_s, time_end_s, time1, time2, tot1, tot2
    real(RK)               :: merr1, merr2, nerr1, nerr2
    integer                :: i
    tot1 = 0.0_RK
    tot2 = 0.0_RK
    merr1 = 0.0_RK
    merr2 = 0.0_RK
    nerr1 = 0.0_RK
    nerr2 = 0.0_RK
    print'(A)', "                   implement       intrinsic"
    print'(A)', "            --------------------------------"
    do i = 1, M
      call RANDOM_NUMBER(A)
      A = EXP((A - 0.5_RK) * 10.0_RK)
!
      call CPU_TIME(time_begin_s)
      B = invsqrt(A)
      call CPU_TIME(time_end_s)
      time1 = 1000 * (time_end_s - time_begin_s)
      tot1 = tot1 + time1
      merr1 = merr1 + mean_error(N, A, B)
      nerr1 = nerr1 + norm_error(N, A, B)
!
      call CPU_TIME(time_begin_s)
      C = 1.0_RK / SQRT(A)
      call CPU_TIME(time_end_s)
      time2 = 1000 * (time_end_s - time_begin_s)
      tot2 = tot2 + time2
      merr2 = merr2 + mean_error(N, A, C)
      nerr2 = nerr2 + norm_error(N, A, C)
      print '(I12,2f16.3)', i, time1, time2
    end do
    print'(A)', "            --------------------------------"
    print'(12X,2f16.3)', tot1, tot2
    print'(A,2f16.9)', "  mean error", merr1 / M, merr2 / M
    print'(A,2f16.9)', "  norm error", nerr1 / M, nerr2 / M
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
  end subroutine test2

  subroutine test3()
    integer(IK), parameter :: N = 5000000, M = 10
    real(RK)               :: A(N), B(N)
    real(RK)               :: time_begin_s, time_end_s, time1, time2, tot1, tot2
    integer                :: i
    tot1 = 0.0_RK
    tot2 = 0.0_RK
    print'(A)', " div vs inv   A*invcbrt(A*A)    1/invcbrt(A)"
    print'(A)', "            --------------------------------"
    do i = 1, M
      call RANDOM_NUMBER(A)
      A = EXP((A - 0.5_RK) * 10.0_RK)
!
      call CPU_TIME(time_begin_s)
      B = A * invsqrt(A)
      call CPU_TIME(time_end_s)
      time1 = 1000 * (time_end_s - time_begin_s)
      tot1 = tot1 + time1
!
      call CPU_TIME(time_begin_s)
      B = 1.0_RK / invsqrt(A)
      call CPU_TIME(time_end_s)
      time2 = 1000 * (time_end_s - time_begin_s)
      tot2 = tot2 + time2
      print '(I12,2f16.3)', i, time1, time2
    end do
    print'(A)', "            --------------------------------"
    print'(12X,2f16.3)', tot1, tot2
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
  end subroutine test3
!
  pure function mean_error(n, X, Y) result(res)
    integer, intent(in)  :: n
    real(RK), intent(in) :: X(n), Y(n)
    real(RK)             :: res
    res = SQRT(SUM((X - 1.0_RK / Y**2)**2) / SIZE(X))
  end function mean_error
!
  pure function norm_error(n, X, Y) result(res)
    integer, intent(in)  :: n
    real(RK), intent(in) :: X(n), Y(n)
    real(RK)             :: res
    res = MAXVAL(ABS(X - 1.0_RK / Y**2))
  end function norm_error
!
end program main

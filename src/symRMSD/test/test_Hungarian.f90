program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_Hungarian
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer(IK), parameter :: NTEST = 25
  integer(IK)            :: itest
!
  call u%init('test Hungarian')
  do itest =1,NTEST
    call test1()
  enddo
  do itest =1,NTEST
    call test2()
  enddo
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer(IK), parameter :: N_TEST=10
    integer(IK), parameter :: d=3
    integer(IK), parameter :: n=3
    real(RK)           :: X(d, n), Y(d, n), C(n, n), P(n, n)
    real(RK)           :: minsp
    integer            :: piv(n), i, j
!
    call random_number(X)
    call random_number(Y)
    C = MATMUL(TRANSPOSE(X), Y)
    call Hungarian(n, C, P, piv)
    do concurrent(i=1:n, j=1:n)
      P(i, j) = MERGE(ONE, ZERO, piv(i) == j)
    end do
    minsp = 999
    minsp = MIN(minsp, SP(n, [1_IK, 2_IK, 3_IK], C))
    minsp = MIN(minsp, SP(n, [1_IK, 3_IK, 2_IK], C))
    minsp = MIN(minsp, SP(n, [2_IK, 1_IK, 3_IK], C))
    minsp = MIN(minsp, SP(n, [2_IK, 3_IK, 1_IK], C))
    minsp = MIN(minsp, SP(n, [3_IK, 1_IK, 2_IK], C))
    minsp = MIN(minsp, SP(n, [3_IK, 2_IK, 1_IK], C))
    call u%assert_almost_equal(SUM(C * P), Hungarian_value(n, C), 'CP = X    ')
    call u%assert_almost_equal(Hungarian_value(n, C), minsp,      'X  = minsp')
!
  end subroutine test1
!
  subroutine test2()
    integer, parameter :: N_TEST = 10
    integer, parameter :: n = 1000
    real(RK)           :: C(n, n), P(n, n)
    integer            :: i, j, piv(n)
!
    call RANDOM_NUMBER(C)
    call Hungarian(n, C, P, piv)
    do concurrent(i=1:n, j=1:n)
      P(i, j) = MERGE(ONE, ZERO, piv(i) == j)
    end do
    call u%assert_almost_equal(SUM(C * P), Hungarian_value(n, C), 'CP = X    ')
!
  end subroutine test2
!
  pure function SP(n, ix, C) result(res)
    integer, intent(in)  :: n, ix(n)
    real(RK), intent(in) :: C(n, n)
    integer              :: i
    real(RK)             :: res
    res = 0D0
    do i = 1, n
      res = res + C(i, ix(i))
    end do
  end function SP
!
end program main

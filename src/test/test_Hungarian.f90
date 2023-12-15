program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_Hungarian
  use mod_unittest
  implicit none
  type(unittest) :: z
!
  call z%init('test pca d=3')
  call test1()
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter :: N_TEST=10
    integer, parameter :: d=3
    integer, parameter :: n=10
    real(RK)           :: X(d, n), Y(d, n), C(n, n), P(n, n)
!
    call random_number(X)
    call random_number(Y)
    C = MATMUL(TRANSPOSE(X), Y)
    call Hungarian(n, C, P)
    print'(10f9.3)', P
!
  end subroutine test1
!
end program main

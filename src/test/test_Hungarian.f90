program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_Hungarian
  use mod_unittest
  implicit none
  type(unittest) :: z
!
  call z%init('test pca d=3')
  call test1()
  call test2()
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter :: N_TEST=10
    integer, parameter :: d=3
    integer, parameter :: n=3
    integer, parameter :: i1(3) = [1,2,3]
    integer, parameter :: i2(3) = [1,3,2]
    integer, parameter :: i3(3) = [2,1,3]
    integer, parameter :: i4(3) = [2,3,1]
    integer, parameter :: i5(3) = [3,1,2]
    integer, parameter :: i6(3) = [3,2,1]
    real(RK)           :: X(d, n), Y(d, n), C(n, n), P(n, n)
    integer            :: piv(n), i, j
!
    call random_number(X)
    call random_number(Y)
    C = MATMUL(TRANSPOSE(X), Y)
    call Hungarian(n, C, piv, P)
    do concurrent(i=1:n, j=1:n)
      P(i, j) = MERGE(ONE, ZERO, i == piv(j))
    end do
    print'(3f9.3)', C,P
    print'(3f9.3)', SUM(C * TRANSPOSE(P))
    print'(3f9.3)', Hungarian_value(n, C)
    print*
    print'(3f9.3)', matmul(C, P)
    print*
    print'(3f9.3)', matmul(P, C)
    print*
    print'(3i4,f9.6)',i1, SP(n, i1, C)
    print'(3i4,f9.6)',i2, SP(n, i2, C)
    print'(3i4,f9.6)',i3, SP(n, i3, C)
    print'(3i4,f9.6)',i4, SP(n, i4, C)
    print'(3i4,f9.6)',i5, SP(n, i5, C)
    print'(3i4,f9.6)',i6, SP(n, i6, C)
!
  end subroutine test1
!
  subroutine test2()
    integer, parameter :: N_TEST=10
    integer, parameter :: n=1000
    real(RK)           :: C(n, n), P(n, n)
    integer            :: i, j, piv(n)
!
    call random_number(C)
    call Hungarian(n, C, piv, P)
    do concurrent(i=1:n, j=1:n)
    P(i, j) = MERGE(ONE, ZERO, i == piv(j))
  end do
    print'(f9.3)', SUM(C * TRANSPOSE(P))
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

program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_Kabsch
  use mod_unittest
  implicit none
  type(unittest) :: z
!
  call test1()
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter :: N_TEST=10
    real(RK)           :: Y(18), X(18), rot(3,3), krot(3,3)
    integer            :: i
!
    call z%init('test Kabsch d=3')
!
    call random_number(X)
!
    do i=1,N_TEST
      rot = SO3()
      Y = [MATMUL(rot, RESHAPE(X, [3, 6]))]
      call Kabsch(3, 6, X, Y, krot)
      print'(3f9.3)',rot, transpose(krot)
      print*
!     call z%assert_almost_equal(res, X(1), 'det d=1')
!     call z%assert_almost_equal(sig, SIGN(ONE, X(1)), 'det_sign d=1')
    enddo
!
  end subroutine test1
!
  function SO2() result(res)
    real(RK) :: a(1), res(2, 2)
    call RANDOM_NUMBER(a)
    res(:, 1) = [ COS(a(1)),-SIN(a(1))]
    res(:, 2) = [ SIN(a(1)), COS(a(1))]
  end function SO2
!
  function SO3() result(res)
    real(RK) :: a(3), res(3, 3)
    call RANDOM_NUMBER(a)
    a = a / SQRT(DOT_PRODUCT(a, a))
    res(:,1) = [a(1)**2,        a(1)*a(2)-a(3), a(1)*a(3)+a(2)]
    res(:,2) = [a(1)*a(2)+a(3), a(2)*a(2),      a(2)*a(3)-a(1)]
    res(:,3) = [a(1)*a(3)-a(2), a(2)*a(3)+a(1), a(3)*a(3)     ]
  end function SO3
!
end program main

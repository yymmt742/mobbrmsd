program main
  use mod_params, only: RK, IK
  use mod_rot
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
!
    real(RK), parameter    :: X(12) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    real(RK)               :: Y(12)
    real(RK)               :: U2(2,2), U3(3,3)
    integer                :: i
!
    call u%init('test rot')
!
    do i=1,10
      U2 = SO2()
      call rot(2, 6, U2, X, Y)
      call u%assert_almost_equal(Y, [MATMUL(U2, RESHAPE(X, [2, 6]))], 'rot  2d  ')
      call u%assert_almost_equal(DOT_PRODUCT(Y( 1: 2),Y( 1: 2)),  5D0, 'norm 2d 1')
      call u%assert_almost_equal(DOT_PRODUCT(Y( 3: 4),Y( 3: 4)), 25D0, 'norm 2d 2')
      call u%assert_almost_equal(DOT_PRODUCT(Y( 5: 6),Y( 5: 6)), 61D0, 'norm 2d 3')
      call u%assert_almost_equal(DOT_PRODUCT(Y( 7: 8),Y( 7: 8)),113D0, 'norm 2d 4')
      call u%assert_almost_equal(DOT_PRODUCT(Y( 9:10),Y( 9:10)),181D0, 'norm 2d 5')
      call u%assert_almost_equal(DOT_PRODUCT(Y(11:12),Y(11:12)),265D0, 'norm 2d 6')
    enddo
!
    do i=1,10
      U3 = SO3()
      call rot(3, 4, U3, X, Y)
      call u%assert_almost_equal(Y, [MATMUL(U3, RESHAPE(X, [3, 4]))], 'rot  3d  ')
      call u%assert_almost_equal(DOT_PRODUCT(Y( 1: 3),Y( 1: 3)), 14D0, 'norm 3d 1')
      call u%assert_almost_equal(DOT_PRODUCT(Y( 4: 6),Y( 4: 6)), 77D0, 'norm 3d 2')
      call u%assert_almost_equal(DOT_PRODUCT(Y( 7: 9),Y( 7: 9)),194D0, 'norm 3d 3')
      call u%assert_almost_equal(DOT_PRODUCT(Y(10:12),Y(10:12)),365D0, 'norm 3d 4')
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

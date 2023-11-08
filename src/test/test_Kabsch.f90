program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_Kabsch
  use mod_unittest
  implicit none
  type(unittest) :: z
!
! call test1()
  call test2()
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter    :: N_TEST=10
    real(RK)              :: Y(18), X(18)
    real(RK)              :: rot(3, 3), krot(3, 3)
    real(RK)              :: tor(6, 6), ktor(6, 6)
    real(RK)              :: E3(3,3), E6(6,6)
    real(RK), allocatable :: w(:)
    integer               :: i, j, lw
!
    call z%init('test Kabsch d=3')
!
    call random_number(X)
!
    lw = Kabsch_worksize(6)
!
    allocate(w(lw))
!
    do concurrent(j=1:3, i=1:3)
      E3(i, j) = MERGE(1D0, 0D0, i == j)
    enddo
    do concurrent(j=1:6, i=1:6)
      E6(i, j) = MERGE(1D0, 0D0, i == j)
    enddo
!
    do i=1,N_TEST
!
      rot = SO3()
      Y = [MATMUL(rot, RESHAPE(X, [3, 6]))]
      call Kabsch(3, 6, X, Y, krot, w)
      print'(3f9.3)',krot
      call z%assert_almost_equal([rot - TRANSPOSE(krot)], 0D0, 'Krot = rot^{-1}')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - E3], 0D0, 'Krot@Krot^{T}=I')
!
      tor = SO6()
      tor = E6
      print'(6f9.3)',matmul(transpose(tor),tor)
      Y = [MATMUL(RESHAPE(X, [3, 6]), TRANSPOSE(tor))]
      call Kabsch(3, 5, X, Y, ktor, w, row_major=.TRUE.)
!
      print'(5f9.3)',ktor
      print*
!     print'(6f9.3)',tor
!     print*
!
      !call z%assert_almost_equal([rot - TRANSPOSE(krot)], 0D0, 'Krot = rot^{-1}')
    enddo
!
  end subroutine test1
!
  subroutine test2()
    integer, parameter    :: N_TEST=10
    real(RK)              :: Y(18), X(18)
    real(RK)              :: rot(6, 6), krot(6, 6)
    real(RK)              :: E6(6,6)
    real(RK), allocatable :: w(:)
    integer               :: i, j, lw
!
    call z%init('test Kabsch d=6')
!
    call random_number(X)
!
    lw = Kabsch_worksize(6)
!
    allocate(w(lw))
!
    do concurrent(j=1:6, i=1:6)
      E6(i, j) = MERGE(1D0, 0D0, i == j)
    enddo
!
    do i=1,N_TEST
!
      rot = SO6()
      rot = E6
      Y = [MATMUL(rot, RESHAPE(X, [6, 3]))]
      call Kabsch(6, 3, X, Y, krot, w)
      print'(6f9.3)',rot
      print*
      print'(6f9.3)',krot
      print*
!      call z%assert_almost_equal([rot - TRANSPOSE(krot)], 0D0, 'Krot = rot^{-1}')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - E6], 0D0, 'Krot@Krot^{T}=I')
!
    enddo
!
  end subroutine test2
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
  function SO6() result(res)
    real(RK) :: res(6, 6), tmp(6, 6)
    res = 0D0
    res(1:2,3:4) = SO2()
    res(3:4,1:2) = SO2()
    res(5:6,5:6) = SO2()
!
    tmp = 0D0
    tmp(4:6,1:3) = SO3()
    tmp(1:3,4:6) = SO3()
!
    res = MATMUL(res, tmp)
!
  end function SO6
!
end program main

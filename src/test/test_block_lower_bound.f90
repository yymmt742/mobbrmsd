program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_rmsd
  use mod_rmsd_brute
  use mod_block_lower_bound
  use mod_unittest
  implicit none
  type(unittest)    :: u
  integer,parameter :: NTEST=1000
  integer           :: fail
  integer           :: i
!
  call u%init('test block_lower_bound')
  fail = 0
  do i=1,NTEST
    call test1(fail)
  enddo
  print*,fail,'/',NTEST
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(fail)
    integer, intent(inout) :: fail
    integer, parameter     :: d = 3
    integer, parameter     :: m = 5
    integer, parameter     :: n = 3
    integer, parameter     :: mlist(3) = [2, 3, 4]
    integer, parameter     :: nlist(2) = [1, 3]
    real(RK)               :: X(d * m * n), Y(d * m * n)
    real(RK)               :: w(block_lower_bound_worksize(d, m, n, mlist, nlist))
    real(RK)               :: b
!
    X = [sample(d, m * n)]
    Y = [MATMUL(MATMUL(SO3(), RESHAPE(X, [d, m * n])), SO15())]
!
!   call block_lower_bound(d, m, n, mlist, nlist, X, Y, w)
    b = w(1)
!   call block_lower_bound(d, m, n, mlist, nlist, X, Y, w, R0=SO3()) ; if(w(1)<b) b = w(1)
!   call block_lower_bound(d, m, n, mlist, nlist, X, Y, w, R0=SO3()) ; if(w(1)<b) b = w(1)
!   call block_lower_bound(d, m, n, mlist, nlist, X, Y, w, R0=SO3()) ; if(w(1)<b) b = w(1)
!   call block_lower_bound(d, m, n, mlist, nlist, X, Y, w, R0=SO3()) ; if(w(1)<b) b = w(1)
!
    if (b > 0.0001D0) then
      fail = fail + 1
      print'(I8,F9.6)', fail, w(1)
    endif
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
    res(:, 1) = [a(1) * a(1),        a(1) * a(2) - a(3), a(1) * a(3) + a(2)]
    res(:, 2) = [a(1) * a(2) + a(3), a(2) * a(2),        a(2) * a(3) - a(1)]
    res(:, 3) = [a(1) * a(3) - a(2), a(2) * a(3) + a(1), a(3) * a(3)]
  end function SO3
!
  function SO6() result(res)
    real(RK) :: res(6, 6), tmp(6, 6)
!
    res = 0D0
    res(1:2,3:4) = SO2()
    res(3:4,1:2) = SO2()
    res(5:6,5:6) = SO2()
!
    tmp = 0D0
    tmp(4:6,1:3) = SO3()
    tmp(1:3,4:6) = - SO3()
!
    res = MATMUL(res, tmp)
!
  end function SO6
!
  function SO15() result(res)
    real(RK) :: res(15, 15)
!
    res = 0D0
!
!   res(1,11) = 1D0
!   res(2,12) = 1D0
!   res(3,13) = 1D0
!   res(4,14) = 1D0
!   res(5,15) = 1D0
!
!   res(6,6) = 1D0
!   res(7,7) = 1D0
!   res(8,8) = 1D0
!   res(9,9) = 1D0
!   res(10,10) = 1D0
!
!   res(11,1) =-1D0
!   res(12,2) =-1D0
!   res(13,3) =-1D0
!   res(14,4) =-1D0
!   res(15,5) =-1D0
!
    res(6:10,6:10) = eye(5)
!
    res(1,11) = 1D0
    res(2:4,12:14) = SO3()
    res(5,15) = 1D0
!
    res(11,1) = 1D0
    res(12:14,2:4) = SO3()
    res(15,5) = 1D0
!
  end function SO15
!
  pure function eye(d) result(res)
    integer,intent(in) :: d
    real(RK)           :: res(d, d)
    integer            :: i, j
    do concurrent(j=1:d, i=1:d)
      res(i, j) = MERGE(1D0, 0D0, i == j)
    enddo
  end function eye
!
  function sample(d, n) result(res)
    integer, intent(in)  :: d, n
    real(RK)             :: cnt(d)
    real(RK)             :: res(d, n)
    integer              :: i
    call RANDOM_NUMBER(res)
    cnt = centroid(d, n, res)
    do concurrent(i=1:n)
      res(:, i) = res(:, i) - cnt
    enddo
  end function sample
!
  pure function centroid(d, n, X) result(res)
    integer, intent(in)  :: d, n
    real(RK), intent(in) :: X(d, n)
    real(RK)             :: res(d)
    res = SUM(X, 2) / n
  end function centroid
!
end program main

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
! fail = 0
! do i = 1, NTEST
!   call test1(fail)
! end do
! print*,fail,'/',NTEST
!
  fail = 0
  do i = 1, NTEST
    call test2(fail)
  enddo
  print*,fail,'/',NTEST
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(fail)
    integer, intent(inout)     :: fail
    integer, parameter         :: d = 3
    integer, parameter         :: m = 5
    integer, parameter         :: n = 3
    integer, parameter         :: f = 3
    integer, parameter         :: g = 2
    type(mol_block), parameter :: b(1) = [mol_block(m, n, f, g)]
    real(RK)                   :: X(d * m * n), Y(d * m * n)
    real(RK)                   :: w(block_lower_bound_worksize(d, 1, b))
    real(RK)                   :: score
!
    X = [sample(d, m * n)]
    Y = [MATMUL(MATMUL(SO3(), RESHAPE(X, [d, m * n])), PER15())]
!
    call block_lower_bound(d, 1, [b], X, Y, w)
    score = w(1)
!
    if (score > 0.0001D0) then
      fail = fail + 1
      print'(I8,F9.6)', fail, w(1)
    endif
!
  end subroutine test1
!
  subroutine test2(fail)
    integer, intent(inout)     :: fail
    integer, parameter         :: d = 3
    integer, parameter         :: m = 5
    integer, parameter         :: n = 3
    integer, parameter         :: f = 3
    integer, parameter         :: g = 2
    real(RK), parameter        :: lambda = 5.0_RK
    type(mol_block), parameter :: b(1) = [mol_block(m, n, f, g)]
    real(RK)                   :: X(d * m * n), Y(d * m * n)
    real(RK)                   :: w(block_lower_bound_worksize(d, 1, b))
    real(RK)                   :: score
    integer                    :: i
!
    X = [([lambda * i + 0.5 * sample(d, m)], i=1, n)]
    call centering(d, m * n, X)
    Y = [MATMUL(MATMUL(SO3(), RESHAPE(X, [d, m * n])), PER15())]
!
    call block_lower_bound(d, 1, b, X, Y, w, threshold=0.000001D0, maxiter=100)
    !call block_lower_bound(d, 1, b, X, Y, w)
    score = w(1)
!
    if (score > 0.01D0) then
      fail = fail + 1
      print'(I8,F9.6)', fail, w(1)
    endif
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
  function PER3() result(res)
    real(RK) :: res(3, 3)
    integer  :: r(3)
    integer  :: i, j
    call RANDOM_NUMBER(res(:, 1))
    r(1) = MAXLOC(res(:, 1), 1)
    r(2) = MAXLOC([res(:r(1)-1, 1),-100D0,res(r(1)+1:, 1)], 1)
    r(3) = MINLOC(res(:, 1), 1)
    do concurrent(i=1:3, j=1:3)
      res(i, j) = MERGE(1, 0, r(i)==j)
    end do
  end function PER3
!
  function PER15() result(res)
    real(RK) :: res(15, 15)
!
    res = 0D0
!
!   res(6:10,1:5) = eye(5)
!   res(1:5,6:10) = eye(5)
!
    res(1:3,6:8) = PER3()
    res(4:5,9:10) = eye(2)
!
!   res(1:2,6:7) = SO2()
!   res(3:5,8:10) = eye(3)
!
    res(6:8,1:3) = PER3()
    res(9:10,4:5) = eye(2)
!
!   res(6:7,1:2) = SO2()
!   res(8:10,3:5) = eye(3)
!
    res(11:15,11:15) = eye(5)
!
  end function PER15
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
  pure subroutine centering(d, n, X)
    integer, intent(in)     :: d, n
    real(RK), intent(inout) :: X(d, n)
    real(RK)                :: c(d)
    integer(IK)             :: i
    c = centroid(d, n, X)
    do concurrent(i=1:n)
      X(:,i) = X(:,i) - c
    end do
  end subroutine centering
!
  pure function centroid(d, n, X) result(res)
    integer, intent(in)  :: d, n
    real(RK), intent(in) :: X(d, n)
    real(RK)             :: res(d)
    res = SUM(X, 2) / n
  end function centroid
!
end program main

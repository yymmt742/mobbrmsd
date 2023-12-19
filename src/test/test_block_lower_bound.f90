program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_rmsd
  use mod_rmsd_brute
  use mod_mol_block
  use mod_block_lower_bound
  use mod_unittest
  implicit none
  type(unittest)    :: u
  integer, parameter :: NTEST = 1000
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
    call test2(i, fail)
  end do
  print *, fail, '/', NTEST
!
! fail = 0
! do i = 1, NTEST
!   call test3(fail)
! enddo
! print*,fail,'/',NTEST
!
! call test4()
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
    integer, parameter     :: f = 3
    integer, parameter     :: g = 2
    type(mol_block)        :: b = mol_block(0, m, n, f, g)
    type(mol_block_list)   :: blk
    real(RK)               :: X(d * m * n), Y(d * m * n)
    real(RK), allocatable  :: w(:)
!
    blk = mol_block_list(d, 1, [b]) ! g = n
    blk%b(1)%g = g
    allocate (w(block_lower_bound_worksize(blk)))
    X = [sample(d, m * n)]
    call centering(d, m * n, X)
    Y = [MATMUL(MATMUL(SO3(), RESHAPE(X, [d, m * n])), PER15())]
    call centering(d, m * n, Y)
!
    call block_lower_bound(blk, X, Y, w, nrand=10)
!
    if (w(1) > 0.0001D0) then
      fail = fail + 1
      print'(I8,F9.6)', fail, w(1)
    end if
!
  end subroutine test1
!
  subroutine test2(itest, fail)
    integer, intent(in)    :: itest
    integer, intent(inout) :: fail
    integer, parameter     :: d = 3
    integer, parameter     :: m = 15
    integer, parameter     :: n = 15
    integer, parameter     :: f = 3
    integer, parameter     :: g = 8
    real(RK), parameter    :: lambda = 2.0_RK
    type(mol_block)        :: b = mol_block(0, m, n, f, g)
    type(mol_block_list)   :: blk
    real(RK)               :: X(d * m * n), Y(d * m * n)
    real(RK), allocatable  :: w(:)
    integer                :: i
!
    blk = mol_block_list(d, 1, [b]) ! g = n
    allocate (w(block_lower_bound_worksize(blk)))
    X = [([lambda * i + sample(d, m)], i=1, n)]
    call centering(d, m * n, X)
    Y = [MATMUL(MATMUL(SO3(), RESHAPE(X, [d, m * n])), PER25(m, n, g))]
    call centering(d, m * n, Y)
!
    call block_lower_bound(blk, X, Y, w, nrand=5)
!
    if (w(1) > 0.01D0) then
      fail = fail + 1
      print'(I8,A,I8,2F9.4)', fail, '/', itest, w(1), real(fail, RK) / real(itest, RK)
    end if
!
  end subroutine test2
!
! subroutine test3(fail)
!   integer, intent(inout) :: fail
!   integer, parameter     :: d = 3
!   integer, parameter     :: m = 5
!   integer, parameter     :: n = 12
!   integer, parameter     :: f = 3
!   integer, parameter     :: g = 2
!   real(RK), parameter    :: lambda = 2.0_RK
!   real(RK)               :: X(d * m * n * 2), Y(d * m * n * 2)
!   type(mol_block_list)   :: b
!   real(RK), allocatable  :: w(:)
!   integer                :: i, j, k
!
!   b = mol_block_list(d, 2, [m, m], [n, n], [f, f], [n, n]) ! g = n
!   b%b(1)%g = g
!   b%b(2)%g = g
!
!   X = [([([lambda * RESHAPE([([i, j, 0], k=1, m)], [d, m]) + sample(d, m)], i=1, n)], j=0, 1)]
!   call centering(d, 2 * m * n, X)
!   Y = [MATMUL(SO3(), RESHAPE([(MATMUL(RESHAPE(X(d * m * n * j + 1:d * m * n * (j + 1)), [d, m * n]), &
!     & PER25(m, n, g)), j=0, 1)], [d, m * n * 2]))]
!   call centering(d, 2 * m * n, Y)
!
!   allocate (w(block_lower_bound_worksize(b)))
!   call block_lower_bound(b, X, Y, w, nrand=5)
!
!   if (w(1) > 0.01D0) then
!     fail = fail + 1
!     print'(I8,F9.4)', fail, w(1)
!   end if
!
! end subroutine test3
!
! subroutine test4()
!   integer, parameter     :: d = 3
!   integer, parameter     :: m = 5
!   integer, parameter     :: n = 5
!   integer, parameter     :: f = 3
!   integer, parameter     :: g = 2
!   real(RK), parameter    :: lambda = 5.0_RK
!   real(RK)               :: X(d, m * n * 2), Y(d, m * n * 2)
!   type(mol_block_list)   :: b
!   real(RK), allocatable  :: w(:)
!   integer                :: i, j, k
!
!   b = mol_block_list(d, 2, [m, m], [n, n], [f, f], [n, n]) ! g = n
!   b%b(1)%g = g
!   b%b(2)%g = 0
!
!   X = RESHAPE([([([lambda * RESHAPE([([i, j, 0], k=1, m)], [d, m]) + sample(d, m)], i=1, n)], j=0, 1)], [d, m * n * 2])
!   call centering(d, 2 * m * n, X)
!   Y = MATMUL(SO3(), X)
!   Y(:, :m * n) = MATMUL(Y(:, :m * n), PER25(m, n, g))
!   call centering(d, 2 * m * n, Y)
!
!   allocate (w(block_lower_bound_worksize(b)))
!   call block_lower_bound(b, X, Y, w, nrand=5)
!   print *, w(1)
!
! end subroutine test4
!
  function SO3() result(res)
    real(RK) :: a(3), res(3, 3)
    call RANDOM_NUMBER(a)
    a = a / SQRT(DOT_PRODUCT(a, a))
    res(:, 1) = [a(1) * a(1), a(1) * a(2) - a(3), a(1) * a(3) + a(2)]
    res(:, 2) = [a(1) * a(2) + a(3), a(2) * a(2), a(2) * a(3) - a(1)]
    res(:, 3) = [a(1) * a(3) - a(2), a(2) * a(3) + a(1), a(3) * a(3)]
  end function SO3
!
  function PER3() result(res)
    real(RK) :: res(3, 3)
    integer  :: r(3)
    integer  :: i, j
    call RANDOM_NUMBER(res(:, 1))
    r(1) = MAXLOC(res(:, 1), 1)
    r(2) = MAXLOC([res(:r(1) - 1, 1), -100D0, res(r(1) + 1:, 1)], 1)
    r(3) = MINLOC(res(:, 1), 1)
    do concurrent(i=1:3, j=1:3)
      res(i, j) = MERGE(1, 0, r(i) == j)
    end do
  end function PER3
!
  function PER15() result(res)
    real(RK) :: res(15, 15), a
!
    call RANDOM_NUMBER(a)
    res = 0D0
!
    if (a < 0.5D0) then
      res(1:3, 1:3) = PER3()
      res(4:5, 4:5) = eye(2)
      res(6:8, 6:8) = PER3()
      res(9:10, 9:10) = eye(2)
    else
      res(1:3, 6:8) = PER3()
      res(4:5, 9:10) = eye(2)
      res(6:8, 1:3) = PER3()
      res(9:10, 4:5) = eye(2)
    end if
!
    res(11:15, 11:15) = eye(5)
!
  end function PER15
!
  function PER25(m, n, g) result(res)
    integer(IK), intent(in) :: m, n, g
    real(RK) :: res(m * n, m * n)
    integer(IK) :: i
!
    res = eye(m * n)
!
    do i = 0, g - 1
      res(i * m + 1:i * m + 3, i * m + 1:i * m + 3) = PER3()
    end do
!
  end function PER25
!
  pure function eye(d) result(res)
    integer, intent(in) :: d
    real(RK)           :: res(d, d)
    integer            :: i, j
    do concurrent(j=1:d, i=1:d)
      res(i, j) = MERGE(1D0, 0D0, i == j)
    end do
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
    end do
  end function sample
!
  pure subroutine centering(d, n, X)
    integer, intent(in)     :: d, n
    real(RK), intent(inout) :: X(d, n)
    real(RK)                :: c(d)
    integer(IK)             :: i
    c = centroid(d, n, X)
    do concurrent(i=1:n)
      X(:, i) = X(:, i) - c
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

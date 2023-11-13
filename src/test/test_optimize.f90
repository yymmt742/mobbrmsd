program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_cov
  use mod_rmsd
  use mod_Kabsch
  use mod_procrustes
  use mod_unittest
  implicit none
  type(unittest)    :: u
  integer,parameter :: NTEST=5000
  integer           :: fail
  integer           :: i
!
  call u%init('test optimimze')
  fail = 0
  do i=1,NTEST
    call test1(fail)
  enddo
  print*,fail,'/',NTEST
  call u%init('test partial optimimze')
  fail = 0
  do i=1,NTEST
    call test2(fail)
  enddo
  print*,fail,'/',NTEST
  call u%init('test partial molecular optimimze')
  fail = 0
  do i=1,NTEST
    call test3(fail)
  enddo
  print*,fail,'/',NTEST
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(fail)
    integer, intent(inout) :: fail
    integer, parameter     :: n_iter = 500
    integer, parameter     :: d = 3
    integer, parameter     :: n = 6
    real(RK)               :: X(d * n), Y(d * n), Z(d * n)
    real(RK)               :: Xbar(d)
    real(RK)               :: dcov(d * d), ncov(n * n)
    real(RK)               :: drot(d * d), nrot(n * n)
    real(RK)               :: w(Kabsch_worksize(n))
    real(RK)               :: conv, prev
    integer                :: i
!
    call RANDOM_NUMBER(X); Xbar = centroid(d, n, X)
    do i = 1, n
      X((i - 1) * d + 1:i * d) = X((i - 1) * d + 1:i * d) - Xbar
    end do
    Y = [MATMUL(MATMUL(SO3(), RESHAPE(X, [d, n])), SO6())]
    Z = Y
!
    prev = rmsd(d, n, X, Y)
!
    do i=1,n_iter
!
      call cov(d, n, X, Z, dcov, reset=.true.)
      call Kabsch(3, n, dcov, drot, w)
!
      Z = [MATMUL(RESHAPE(drot, [d, d]), RESHAPE(Y, [d, n]))]
      call cov_row_major(d, n, X, Z, ncov, reset=.true.)
      call procrustes(n, ncov, nrot, w)
      Z = [MATMUL(RESHAPE(Z, [d, n]), TRANSPOSE(RESHAPE(nrot, [n, n])))]
      conv = rmsd(d, n, X, Z)
      if (ABS(prev - conv) < 1D-16) exit
      prev = conv
      Z = [MATMUL(RESHAPE(Y, [d, n]), TRANSPOSE(RESHAPE(nrot, [n, n])))]
!
    enddo
!
    Z = [MATMUL(MATMUL(RESHAPE(drot, [d, d]), RESHAPE(Y, [d, n])), TRANSPOSE(RESHAPE(nrot, [n, n])))]
    if(rmsd(d, n, X, Z)>0.0001D0) fail = fail + 1
!   call u%assert_almost_equal([X - Z], 0D0, 'R = RY')
!
  end subroutine test1
!
  subroutine test2(fail)
    integer, intent(inout) :: fail
    integer, parameter  :: n_iter = 500
    integer, parameter  :: d = 3
    integer, parameter  :: n = 12
    integer, parameter  :: nlist(6) = [1,2,3,4,5,6]
    real(RK)            :: X(d * n), Y(d * n), Z(d * n)
    real(RK)            :: Xbar(d)
    real(RK)            :: dcov(d * d), ncov(6 * 6)
    real(RK)            :: drot(d * d), nrot(6 * 6)
    real(RK)            :: w(Kabsch_worksize(n))
    real(RK)            :: conv, prev
    integer             :: i
!
    call RANDOM_NUMBER(X); Xbar = centroid(d, n, X)
    do i = 1, n
      X((i - 1) * d + 1:i * d) = X((i - 1) * d + 1:i * d) - Xbar
    end do
    Y = [MATMUL(MATMUL(SO3(), RESHAPE(X, [d, n])), SO12())]
    Z = Y
!
    prev = rmsd(d, n, X, Y)
!
    do i=1,n_iter
!
      call cov(d, n, X, Z, dcov, reset=.true.)
      call Kabsch(3, n, dcov, drot, w)
!
      Z = [MATMUL(RESHAPE(drot, [d, d]), RESHAPE(Y, [d, n]))]
      call cov_row_major(d, nlist, X, Z, ncov, reset=.true.)
      call procrustes(6, ncov, nrot, w)
      Z = [MATMUL(RESHAPE(Z, [d, n]), RI(nrot))]
      conv = rmsd(d, n, X, Z)
      !print*,i, prev, prev-conv
      if (ABS(prev - conv) < 1D-16) exit
      prev = conv
      Z = [MATMUL(RESHAPE(Y, [d, n]), RI(nrot))]
!
    enddo
!
    Z = [MATMUL(MATMUL(RESHAPE(drot, [d, d]), RESHAPE(Y, [d, n])), RI(nrot))]
    if(rmsd(d, n, X, Z)>0.0001D0) fail = fail + 1
!   call u%assert_almost_equal([X - Z], 0D0, 'R = RY')
!
  end subroutine test2
!
  subroutine test3(fail)
    integer, intent(inout) :: fail
    integer, parameter  :: d = 3
    integer, parameter  :: m = 3
    integer, parameter  :: n = 26
    integer, parameter  :: p = 6
    integer, parameter  :: n_iter = 500
    integer, parameter  :: nlist(p) = [21,22,23,24,25,26]
    real(RK)            :: X(d * m * n), Y(d * m * n), Z(d * m * n)
    real(RK)            :: Xbar(d)
    real(RK)            :: dcov(d * d), ncov(p * p)
    real(RK)            :: drot(d * d), nrot(p * p)
    real(RK)            :: w(Kabsch_worksize(n))
    real(RK)            :: conv, prev
    integer             :: i
!
    call RANDOM_NUMBER(X); Xbar = centroid(d, n, X)
    do i = 1, m * n
      X((i - 1) * d + 1:i * d) = X((i - 1) * d + 1:i * d) - Xbar
    end do
    call mol_rot(d, m, n, MATMUL(SO3(), RESHAPE(X, [d, m * n])), Y, SO26())
    Z = Y
!
    prev = rmsd(d, m * n, X, Y)
!
    do i=1,n_iter
!
      call cov(d, n, X, Z, dcov, reset=.true.)
      call Kabsch(d, m * n, dcov, drot, w)
      Z = [MATMUL(RESHAPE(drot, [d, d]), RESHAPE(Y, [d, m * n]))]
      call cov_row_major(d * m, nlist, X, Z, ncov, reset=.true.)
      call procrustes(p, ncov, nrot, w)
      call mol_rot(d, m, n, [MATMUL(RESHAPE(drot, [d, d]), RESHAPE(Y, [d, m * n]))], Z, RI2(nrot))
      conv = rmsd(d, m * n, X, Z)
      !print*,i, prev, prev-conv
      if (ABS(prev - conv) < 1D-16) exit
      prev = conv
      call mol_rot(d, m, n, Y, Z, RI2(nrot))
!
    enddo
!
    call mol_rot(d, m, n, MATMUL(RESHAPE(drot, [d, d]), RESHAPE(Y, [d, m * n])), Z, RI2(nrot))
!   call u%assert_almost_equal([X - Z], 0D0, 'R = RY')
    if(rmsd(d, n, X, Z)>0.0001D0) fail = fail + 1
!
  end subroutine test3
!
  pure subroutine mol_rot(d, m, n, X, Y, R)
    integer, intent(in)     :: d, m, n
    real(RK), intent(in)    :: X(d, m, n), R(n, n)
    real(RK), intent(inout) :: Y(d, m, n)
    integer                 :: i
    do i=1,m
      Y(:, i, :) = MATMUL(X(:, i, :), R)
    enddo
  end subroutine mol_rot
!
  pure function RI(R) result(res)
    real(RK), intent(in) :: R(6, 6)
    real(RK)             :: res(12, 12)
    res(1:6, 1:6) = TRANSPOSE(R)
    res(7:12,1:6) = 0D0
    res(1:6,7:12) = 0D0
    res(7:12,7:12) = eye(6)
  end function RI
!
  pure function RI2(R) result(res)
    real(RK), intent(in) :: R(6, 6)
    real(RK)             :: res(26, 26)
    res = 0D0
    res(1:20,1:20) = eye(20)
    res(21:26, 21:26) = TRANSPOSE(R)
  end function RI2
!
  pure function centroid(d, n, X) result(res)
    integer, intent(in)  :: d, n
    real(RK), intent(in) :: X(d, n)
    real(RK)             :: res(d)
    res = SUM(X, 2) / n
  end function centroid
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
  function SO12() result(res)
    real(RK) :: res(12, 12)
!
    res = 0D0
    res(1:6,1:6) = SO6()
    res(7:12,7:12) = eye(6)
!
  end function SO12
!
  function SO26() result(res)
    real(RK) :: res(26, 26)
!
    res = 0D0
    res(1:20,1:20) = eye(20)
    res(21:26,21:26) = SO6()
!
  end function SO26
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
end program main

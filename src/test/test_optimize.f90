program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_rmsd
  use mod_Kabsch
  use mod_procrustes
  use mod_unittest
  implicit none
  type(unittest) :: z
!
  call z%init('test optimimze')
  call test1(20)
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1(n_test)
    integer, intent(in) :: n_test
    integer, parameter  :: d = 3
    integer, parameter  :: n = 12
    real(RK)            :: Y(d, n), X(d, n)
    real(RK)            :: Ybar(d), Xbar(d)
    real(RK)            :: dcov(d, d), ncov(n, n)
    real(RK)            :: drot(d, d), nrot(n, n)
    real(RK)            :: w(Kabsch_worksize(n))
    integer             :: i
!
    call RANDOM_NUMBER(X); Xbar = centroid(d, n, X)
    do i = 1, n
      X(:, i) = X(:, i) - Xbar
    end do
!   call RANDOM_NUMBER(Y); Ybar = centroid(d, n, Y)
!   do i = 1, n
!     Y(:, i) = Y(:, i) - Ybar
!   end do
    Y = MATMUL(MATMUL(SO3(), X), SO12())

    print *, rmsd(d, n, X, Y)
!
    do i=1,N_TEST
!
      dcov = MATMUL(X, TRANSPOSE(Y))
      call Kabsch(3, dcov, drot, w)
      Y = MATMUL(drot, Y)
      print'(3f7.2)',drot
      print *, rmsd(d, n, X, Y)
      ncov = MATMUL(TRANSPOSE(X), Y)
      call procrustes(n, ncov, nrot, w)
      Y = MATMUL(Y, TRANSPOSE(nrot))
      print'(12f7.2)',nrot
      print *, rmsd(d, n, X, Y)
!
    enddo
!
  end subroutine test1
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
    res(1:6,7:12) = SO6()
    res(7:12,1:6) = - SO6()
!
  end function SO12
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

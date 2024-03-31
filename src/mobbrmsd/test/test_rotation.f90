program main
  use blas_lapack_interface, only : D, DD
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_rotation
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: z
!
  interface
    include 'dgemm.h'
    include 'dgesvd.h'
    include 'dgetrf.h'
  end interface
!
  call z%init('test rotation')
!
  call test1(1, 10)
  call test1(2, 10)
  call test1(3, 10)
  call test1(10, 10)
  call test1(20, 10)
  call test1(100, 10)
!
  call z%finish_and_terminate()
!
contains
!
  subroutine test1(n, n_test)
    integer, intent(in)   :: n, n_test
    real(RK)              :: Y(D, n), X(D, n), cov(D, D), g
    real(RK)              :: rot(D, D), krot(D, D), sd, kd, sm
    real(RK), allocatable :: w(:)
    integer               :: i
!
    allocate (w(MAX(rotation_worksize(), sdmin_worksize())))
!
    X = sample(n)
!
    do i = 1, N_TEST
      rot = SO()
      Y = MATMUL(rot, X)
      g = SUM(X * X) + SUM(Y * Y)
      cov = MATMUL(X, TRANSPOSE(Y))
!
      call estimate_rotation(g, cov, krot, w)
      call z%assert_almost_equal([X - MATMUL(krot, Y)], ZERO, 'X = YR   ')
!
      if (D <= n) call z%assert_almost_equal([MATMUL(rot, krot) - eye()], ZERO, 'S@RT = I ')
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye()], ZERO, 'R@RT = I ')
!
      call estimate_sdmin(g, cov, w)
      call z%assert_almost_equal(w(1), ZERO, 'sdmin=0  ')
    end do
!
    do i = 1, N_TEST
      Y = 0.8 * X + 0.2 * sample(n)
      cov = MATMUL(X, TRANSPOSE(Y))
      g = SUM(X**2) + SUM(Y**2)
      call estimate_sdmin(g, cov, w)
      call Kabsch(cov, krot)
!
      sm = w(1)
      sd = SUM((X - MATMUL(krot, Y))**2)
      kd = SUM(cov * krot)
      kd = g - kd - kd
!
      call z%assert_almost_equal(sm / sd, ONE, 'vs Kabsch', place=4)
!
      call estimate_rotation(g, cov, krot, w)
      call z%assert_greater_equal(SUM(cov * krot), SUM(cov * SO()), 'CR >= CQ ')
!
      sd = SUM((X - MATMUL(krot, Y))**2)
      kd = SUM(cov * krot)
      kd = g - kd - kd
      call z%assert_almost_equal(sm / sd, ONE, 'sdmin-sd ', place=4)
      call z%assert_almost_equal(sm / kd, ONE, 'sdmin-kd ', place=4)
      call z%assert_almost_equal([MATMUL(krot, TRANSPOSE(krot)) - eye()], ZERO, 'R@RT = I ')
    end do
!
  end subroutine test1
!
!| Calculate the rotation matrix R^T from covariance matrix.
  subroutine Kabsch(cov, rot)
    real(RK), intent(in)    :: cov(*)
    real(RK), intent(inout) :: rot(*)
    real(RK)                :: M(D, D), U(D, D), VT(D, D), UVT(D, D), S(D)
    integer(IK)             :: nw, info
!
    nw = worksize_Kabsch()
!
    block
      real(RK) :: w(nw)
      call copy(DD, cov, M)
      call DGESVD('A', 'A', D, D, M, D, S, U, D, VT, D, w, nw, info)
      UVT = MATMUL(U, VT)
      call det_sign(UVT)
      if (UVT(1,1) < ZERO) call neg(d, U(1, D))
      rot(:DD) = [MATMUL(U, VT)]
    end block
!
  end subroutine Kabsch
!
  !| work array size for Kabsch algorithm.
  pure elemental function worksize_Kabsch() result(res)
    real(RK)    :: w(1)
    integer(IK) :: res, info
    call DGESVD('A', 'A', D, D, w, D, w, w, D, w, D, w, -1, info)
    res = NINT(w(1)) + DD * 3 + D
  end function worksize_Kabsch
!
!| calculate determinant sign of square matrix x, with leading dimension.
   pure subroutine det_sign(x)
     real(RK), intent(inout) :: x(*)
     !! square matrix, on exit, x(1) is assigned the determinant sign of x, <br>
     !! and the other elements are undefined.
!
     if (D < 1) then
       return
     elseif (D == 1) then
       x(1) = SIGN(ONE, x(1))
     elseif (D == 2) then
       x(1) = SIGN(ONE, x(1) * x(4) - x(2) * x(3))
     elseif (D == 3) then
       x(1) = SIGN(ONE, x(1) * (x(5) * x(9) - x(8) * x(6)) +&
         &              x(4) * (x(8) * x(3) - x(2) * x(9)) +&
         &              x(7) * (x(2) * x(6) - x(5) * x(3)))
     else
       block
         integer(IK) :: i, j, k, ipiv(D)
         call DGETRF(D, D, x, D, ipiv, j)
         ipiv(1) = COUNT([(ipiv(i) == i, i=1, D)])
         j = 1
         k = D + 1
         do i = 1, D
           if (x(j) <= ZERO) ipiv(1) = ipiv(1) + 1
           j = j + k
         end do
         if (MODULO(ipiv(1), 2) == 0) then
           x(1) = ONE
         else
           x(1) = -ONE
         end if
       end block
     end if
!
   end subroutine det_sign
!
  pure subroutine neg(N, X)
    integer(IK), intent(in) :: N
    real(RK), intent(inout) :: X(*)
    integer(IK)             :: i
    do concurrent(i=1:N)
      X(i) = -X(i)
    end do
  end subroutine neg
!
  pure subroutine copy(N, X, Y)
    integer(IK), intent(in) :: N
    real(RK), intent(in)    :: X(*)
    real(RK), intent(inout) :: Y(*)
    integer(IK)             :: i
    do concurrent(i=1:N)
      Y(i) = X(i)
    end do
  end subroutine copy
end program main

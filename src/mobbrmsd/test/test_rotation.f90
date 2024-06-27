program main
  use mod_dimspec_functions, only: D, DD, setup_dimension
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_rotation
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: z
#ifdef USE_REAL32
  integer, parameter :: place = 1
#else
  integer, parameter :: place = 6
#endif
!
  interface
#ifdef USE_REAL32
    include 'sgesvd.h'
    include 'sgetrf.h'
#elif USE_REAL64
    include 'dgesvd.h'
    include 'dgetrf.h'
#else
    include 'dgesvd.h'
    include 'dgetrf.h'
#endif
  end interface
!
  call setup_dimension(4) ! for gd
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
  subroutine test1(n, n_test)
    integer, intent(in)   :: n, n_test
    real(RK)              :: Y(D, n), X(D, n)
    real(RK)              :: cov(D, D), g
    real(RK)              :: rot(D, D), krot(D, D), nrm, sd, kd, sm
    real(RK), allocatable :: w(:)
    character(4)          :: cn, ct
    integer               :: i
!
    allocate (w(MAX(rotation_worksize(), sdmin_worksize())))
!
    nrm = ONE / n
    X = sample(n)
    write (cn, '(I0)') n
!
    do i = 1, N_TEST
      write (ct, '(I0)') i
      rot = SO()
      Y = MATMUL(rot, X)
      g = SUM(X * X) + SUM(Y * Y)
      cov = MATMUL(X, TRANSPOSE(Y))
!
      call estimate_rotation(g, cov, krot, w)
      call z%assert_almost_equal([X], [MATMUL(krot, Y)], 'X = YR '//cn//ct, place=place)
!
      if (D <= n) call z%assert_is_eye(MATMUL(rot, krot), 'S@RT = I '//cn//ct, place=place)
      call z%assert_is_eye(MATMUL(krot, TRANSPOSE(krot)), 'R@RT = I '//cn//ct, place=place)
!
      call estimate_sdmin(g, cov, w)
      call z%assert_is_zero(nrm * w(1), 'sdmin = 0 '//cn//ct, place=place)
    end do
!
    do i = 1, N_TEST
      Y = 0.8 * X + 0.2 * sample(n)
      cov = MATMUL(X, TRANSPOSE(Y))
      g = SUM(X**2) + SUM(Y**2)
      call estimate_sdmin(g, cov, w)
      call Kabsch(cov, krot)
      sm = nrm * w(1)
      sd = nrm * SUM((X - MATMUL(krot, Y))**2)
      call z%assert_almost_equal(sm, sd, 'vs Kabsch', place=place)
      call estimate_rotation(g, cov, krot, w)
      call z%assert_greater_equal(SUM(cov * krot), SUM(cov * SO()), 'CR >= CQ ')
      sd = nrm * SUM((X - MATMUL(krot, Y))**2)
      kd = SUM(cov * krot)
      kd = nrm * (g - kd - kd)
      call z%assert_almost_equal(sm, sd, 'sdmin-sd ', place=place)
      call z%assert_almost_equal(sm, kd, 'sdmin-kd ', place=place)
      call z%assert_is_eye(MATMUL(krot, TRANSPOSE(krot)), 'R@RT = I ', place=place)
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
#ifdef USE_REAL32
      call SGESVD('A', 'A', D, D, M, D, S, U, D, VT, D, w, nw, info)
#else
      call DGESVD('A', 'A', D, D, M, D, S, U, D, VT, D, w, nw, info)
#endif
      UVT = MATMUL(U, VT)
      call det_sign(UVT)
      if (UVT(1, 1) < ZERO) call neg(d, U(1, D))
      rot(:DD) = [MATMUL(U, VT)]
    end block
!
  end subroutine Kabsch
!
  !| work array size for Kabsch algorithm.
  pure elemental function worksize_Kabsch() result(res)
    real(RK)    :: w(1)
    integer(IK) :: res, info
#ifdef USE_REAL32
    call SGESVD('A', 'A', D, D, w, D, w, w, D, w, D, w, -1, info)
#elif USE_REAL64
    call DGESVD('A', 'A', D, D, w, D, w, w, D, w, D, w, -1, info)
#else
    call DGESVD('A', 'A', D, D, w, D, w, w, D, w, D, w, -1, info)
#endif
    res = NINT(w(1)) + DD * 3 + D
  end function worksize_Kabsch
!
!| calculate determinant sign of square matrix x, with leading dimension.
  subroutine det_sign(x)
    real(RK), intent(inout) :: x(*)
     !! square matrix, on exit, x(1) is assigned the determinant sign of x, <br>
     !! and the other elements are undefined.
    if (D < 1) then
      return
    elseif (D == 1) then
      x(1) = SIGN(ONE, x(1))
    elseif (D == 2) then
      x(1) = SIGN(ONE, x(1) * x(4) - x(2) * x(3))
    elseif (D == 3) then
      x(1) = SIGN(ONE, x(1) * (x(5) * x(9) - x(8) * x(6))&
        &            - x(4) * (x(2) * x(9) - x(8) * x(3))&
        &            + x(7) * (x(2) * x(6) - x(5) * x(3)))
    else
      block
        integer(IK) :: i, j, k, ipiv(D)
#ifdef USE_REAL32
        call SGETRF(D, D, x, D, ipiv, j)
#else
        call DGETRF(D, D, x, D, ipiv, j)
#endif
        ipiv(1) = COUNT([(ipiv(i) == i, i=1, D)])
        j = 1
        k = D + 1
        do i = 1, D
          if (x(j) < ZERO) ipiv(1) = ipiv(1) + 1
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


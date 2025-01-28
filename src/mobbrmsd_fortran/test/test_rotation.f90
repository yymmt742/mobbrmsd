program main
  use mod_dimspec_functions, only: D, setup_dimension
  use mod_params, only: R8, RK, IK, ONE => RONE, ZERO => RZERO
#ifdef USE_REAL32
  use mod_mobbrmsd_lapack, only: SGESVD, SGETRF
#else
  use mod_mobbrmsd_lapack, only: DGESVD, DGETRF
#endif
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
  call setup_dimension(4) ! for xd
  call z%init('test rotation')
!
  call test1(1, 100)
  call test1(2, 100)
  call test1(3, 100)
  call test1(10, 100)
  call test1(20, 100)
  call test1(100, 100)
!
  call z%finish_and_terminate(passing_score=0.99_R8)
!
contains
  subroutine test1(n, n_test)
    integer, intent(in)   :: n, n_test
    real(RK)              :: Y(D, n), X(D, n)
    real(RK)              :: cov(D, D), g
    real(RK)              :: rot(D, D), krot(D, D), nrm, sd, kd, rd, sm, inxy
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
      inxy = inner(D * n, X, MATMUL(krot, Y))
      call z%assert_almost_equal(inxy, 1.0_RK, 'X = YR '//cn//ct, place=place)
      !call z%assert_almost_equal([X], [MATMUL(krot, Y)], 'X = YR '//cn//ct, place=place)
!
      if (D <= n) call z%assert_is_eye(MATMUL(rot, krot), 'S@RT = I '//cn//ct, place=place)
      call z%assert_is_eye(MATMUL(krot, TRANSPOSE(krot)), 'R@RT = I '//cn//ct, place=place)
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
      rd = nrm * SUM((X - MATMUL(krot, Y))**2)
      kd = nrm * (g - SUM(cov * krot) * 2.0_RK)
      call z%assert_greater_equal(SUM(cov * krot), SUM(cov * SO()), 'CR >= CQ ')
      call z%assert_almost_equal(rd, kd, 'rd vs kd ', place=place)
      call z%assert_almost_equal(sd, rd, 'sd vs rd ', place=place)
      call z%assert_almost_equal(sm, rd, 'sdmin-sd ', place=place)
      call z%assert_is_eye(MATMUL(krot, TRANSPOSE(krot)), 'R@RT = I ', place=place)
    end do
!
  end subroutine test1
!
  pure function inner(n, X, Y) result(res)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: X(n), Y(n)
    real(RK)                :: res
    res = DOT_PRODUCT(X, Y) / SQRT(DOT_PRODUCT(X, X) * DOT_PRODUCT(Y, Y))
  end function inner
!
end program main


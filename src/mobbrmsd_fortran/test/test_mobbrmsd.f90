program main
  use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
  use mod_dimspec_functions, only: D, setup_dimension
  use mod_params, only: R8, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mobbrmsd
  use mod_mobbrmsd_state
  use mod_mobbrmsd_mst
  use mod_mobbrmsd_batch_run
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
#ifdef USE_REAL32
  integer, parameter :: place = 1
#else
  integer, parameter :: place = 6
#endif
!
#ifdef USE_DIM2
  integer, parameter :: ntest_def = 20
#elif USE_DIM3
  integer, parameter :: ntest_def = 10
#else
  integer, parameter :: ntest_def = 1
#endif
!
  call setup_dimension(4)
  call u%init('test mobbrmsd for (n,M,S)=(1,1,1)')
  call test1(1, 1, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(1,2,1)')
  call test1(1, 2, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(1,3,1)')
  call test1(1, 3, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(4,1,1)')
  call test1(4, 1, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(4,3,1)')
  call test1(4, 2, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(4,3,1)')
  call test1(4, 3, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(4,1,2)')
  call test1(4, 1, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd for (n,M,S)=(4,2,2)')
  call test1(4, 2, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd for (n,M,S)=(40,6,1)')
  call test1(40, 6, 1, [0])
!
  call u%init('test mobbrmsd for {(n,M,S)}={(5,1,1), (5,4,1)}')
  call test2(5, 1, 1, [0], 5, 4, 1, [0])
  call u%init('test mobbrmsd for {(n,M,S)}={(4,2,1), (5,2,1)}')
  call test2(4, 2, 1, [0], 5, 2, 1, [0])
  call u%init('test mobbrmsd for {(n,M,S)}={(8,2,1), (4,2,2)}')
  call test2(8, 2, 1, [0], 4, 2, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd for {(n,M,S)}={(24,3,1), (24,4,1)}')
  call test2(24, 3, 1, [0], 24, 4, 1, [0])
!
  call u%init('test mobbrmsd repeat for {(n,M,S)}={(4,4,2)}')
  call test3(4, 4, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd repeat for {(n,M,S)}={(4,8,1)}')
  call test3(4, 8, 1, [0])
!
  call u%init('test mobbrmsd difflim for {(n,M,S)}={(4,4,2)}')
  call test4(4, 4, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd difflim for {(n,M,S)}={(4,8,1)}')
  call test4(4, 8, 1, [0])
#ifdef USE_DIMX
#else
  call u%init('test mobbrmsd difflim for {(n,M,S)}={(8,10,1)}')
  call test4(8, 10, 1, [0])
#endif
!
  call u%init('test mobbrmsd absolute difflim for {(n,M,S)}={(4,4,2)}')
  call test5(4, 4, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd absolute difflim for {(n,M,S)}={(4,8,1)}')
  call test5(4, 8, 1, [0])
#ifdef USE_DIMX
#else
  call u%init('test mobbrmsd absolute difflim for {(n,M,S)}={(8,10,1)}')
  call test5(8, 10, 1, [0])
#endif
!
  call u%init('test mobbrmsd cutoff for {(n,M,S)}={(4,4,2)}')
  call test6(4, 4, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd cutoff for {(n,M,S)}={(4,8,1)}')
  call test6(4, 8, 1, [0])
  call u%init('test mobbrmsd cutoff for {(n,M,S)}={(8,10,1)}')
  call test6(8, 10, 1, [0])
!
! call u%init('test mobbrmsd min_span_tree for {(n,M,S)}={(4,10,1)}, n_target=4')
! call test7(4, 10, 1, [0], 4, ntest_def)
! call u%init('test mobbrmsd min_span_tree for {(n,M,S)}={(4,10,1)}, n_target=10')
! call test7(4, 8, 1, [0], 10, ntest_def * 2)
  call u%init('test mobbrmsd min_span_tree for {(n,M,S)}={(4,4,1)}, n_target=100')
  call test7(4, 6, 1, [0], 100, 1)
  !call test7(4, 4, 1, [0], 100, ntest_def * 2)
!
  call u%finish_and_terminate(passing_score=0.95_R8)
!
contains
!
  subroutine test1(n, m, s, sym)
    integer, intent(in)    :: n, m, s, sym(n * (s - 1))
    type(mobbrmsd)         :: mobb
    type(mobbrmsd_state)   :: stat
    type(mol_block_input), allocatable :: inp(:)
    real(RK)               :: X(D, n, m), Y(D, n, m), Z(D, n, m)
    real(RK)               :: R(D, D)
    real(RK)               :: sd
    real(RK), allocatable  :: W(:)
    integer(IK)            :: i
!
    call mol_block_input_add_molecule(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X = sample(n, m)
    call centering(n, m, X)
    Y = X
!
    allocate (W(mobbrmsd_memsize(mobb)))
!
    do i = 1, 20
      call mobbrmsd_run(mobb, stat, X, Y, W, get_rotation=.true.)
      call u%assert_almost_equal( &
     &       mobbrmsd_state_squared_deviation(stat), &
     &       brute_sd(n, m, s, sym, X, Y), &
     &       'minrmsd value', &
     &       place=place &
     &      )
      Z = Y
      call mobbrmsd_state_rotation_matrix(stat, R)
      call mobbrmsd_swap_and_rotation(mobb, stat, Z)
      sd = SUM((X - Z)**2)
      call u%assert_is_eye(MATMUL(R, TRANSPOSE(R)), "RRT = I", place=place)
      call u%assert_almost_equal(mobbrmsd_state_squared_deviation(stat), &
     &                           sd, 'rotation', &
     &                           place=place)
      Y = 0.9 * Y + 0.1 * sample(n, m)
      call centering(n, m, Y)
    end do
!
    deallocate (inp)
!
  end subroutine test1
!
  subroutine test2(n1, m1, s1, sym1, n2, m2, s2, sym2)
    integer, intent(in)   :: n1, m1, s1, sym1(n1 * (s1 - 1))
    integer, intent(in)   :: n2, m2, s2, sym2(n2 * (s2 - 1))
    type(mol_block_input), allocatable :: inp(:)
    type(mobbrmsd)        :: mobb
    type(mobbrmsd_state)  :: stat
    real(RK)              :: brute
    real(RK)              :: X1(D, n1, m1), X2(D, n2, m2)
    real(RK)              :: Y1(D, n1, m1), Y2(D, n2, m2)
    real(RK)              :: X(D, m1 * n1 + m2 * n2)
    real(RK)              :: Y(D, m1 * n1 + m2 * n2)
    real(RK)              :: Z(D, m1 * n1 + m2 * n2)
    real(RK), allocatable :: W(:)
    real(RK)              :: R(D, D), sd
    integer(IK)           :: i
!
    call mol_block_input_add_molecule(inp, n1, m1, sym=RESHAPE(sym1, [n1, s1 - 1]))
    call mol_block_input_add_molecule(inp, n2, m2, sym=RESHAPE(sym2, [n2, s2 - 1]))
    mobb = mobbrmsd(inp)
!
    X1 = sample(n1, m1)
    X2 = sample(n2, m2)
    X = RESHAPE([X1, X2], SHAPE(X))
    call centering(SIZE(X, 2), X)
    Y1 = X1
    Y2 = X2
!
    allocate (W(mobbrmsd_memsize(mobb)))
!
    do i = 1, 10
      Y = RESHAPE([Y1, Y2], SHAPE(Y))
      call centering(SIZE(Y, 2), Y)
      call mobbrmsd_run(mobb, stat, X, Y, W, get_rotation=.true.)
      brute = brute_sd_double(n1, m1, s1, sym1, n2, m2, s2, sym2, X1, Y1, X2, Y2)
      call u%assert_almost_equal(mobbrmsd_state_squared_deviation(stat), brute, 'minrmsd value', place=place)
      Z = Y
      call mobbrmsd_state_rotation_matrix(stat, R)
      call mobbrmsd_swap_and_rotation(mobb, stat, Z)
      sd = SUM((X - Z)**2)
      call u%assert_is_eye(MATMUL(R, TRANSPOSE(R)), "RRT = I", place=place)
      call u%assert_almost_equal(sd, &
     &                           brute, 'rotation', &
     &                           place=place)
      Y1 = 0.5 * Y1 + 0.5 * sample(n1, m1)
      Y2 = 0.5 * Y2 + 0.5 * sample(n2, m2)
    end do
!
  end subroutine test2
!
  subroutine test3(n, m, s, sym)
    integer, intent(in)    :: n, m, s, sym(n * (s - 1))
    type(mobbrmsd)         :: mobb
    type(mobbrmsd_state)   :: stat
    type(mobbrmsd_input)   :: inp
    real(RK)               :: X(D, n, m), Y(D, n, m)
    real(RK), allocatable  :: W(:)
    real(RK)               :: sd, brute, sd2
    integer(IK)            :: i
!
    call mobbrmsd_input_add_molecule(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X = sample(n, m)
    Y = sample(n, m)
!
    allocate (W(mobbrmsd_memsize(mobb)))
!
    call mobbrmsd_run(mobb, stat, X, Y, W, maxeval=0)
!
    do i = 1, 10
      call mobbrmsd_restart(mobb, stat, W, maxeval=0)
      call print_stat(stat, mobbrmsd_state_autovariance(stat))
    end do
!
    call mobbrmsd_restart(mobb, stat, W)
    call print_stat(stat, mobbrmsd_state_autovariance(stat))
    sd = mobbrmsd_state_squared_deviation(stat)
    brute = brute_sd(n, m, s, sym, X, Y)
!
    call mobbrmsd_run(mobb, stat, X, Y)
    sd2 = mobbrmsd_state_squared_deviation(stat)
    call print_stat(stat, mobbrmsd_state_autovariance(stat))
    call u%assert_almost_equal(sd, brute, 'minrmsd value', place=place)
    call u%assert_almost_equal(sd, sd2, 'vs at once   ', place=place)
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
!
  end subroutine test3
!
  subroutine test4(n, m, s, sym)
    integer, intent(in)   :: n, m, s, sym(n * (s - 1))
    type(mobbrmsd)        :: mobb
    type(mobbrmsd_state)  :: stat
    type(mobbrmsd_input)  :: inp
    real(RK)              :: X(D, n, m), Y(D, n, m)
    real(RK), allocatable :: W(:)
!
    call mobbrmsd_input_add_molecule(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X = sample(n, m)
    Y = sample(n, m)
!
    allocate (W(mobbrmsd_memsize(mobb)))
!
    call mobbrmsd_run(mobb, stat, X, Y, W, difflim=1.0_RK)
    call print_stat(stat, 1.0_RK * mobbrmsd_state_autovariance(stat))
    call mobbrmsd_restart(mobb, stat, W, difflim=0.2_RK)
    call print_stat(stat, 0.2_RK * mobbrmsd_state_autovariance(stat))
    call mobbrmsd_restart(mobb, stat, W, difflim=0.1_RK)
    call print_stat(stat, 0.1_RK * mobbrmsd_state_autovariance(stat))
    call mobbrmsd_restart(mobb, stat, W, difflim=0.01_RK)
    call print_stat(stat, 0.01_RK * mobbrmsd_state_autovariance(stat))
    call mobbrmsd_restart(mobb, stat, W, difflim=0.001_RK)
    call print_stat(stat, 0.001_RK * mobbrmsd_state_autovariance(stat))
    call mobbrmsd_restart(mobb, stat, W)
    call print_stat(stat, 0.0_RK)
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
!
  end subroutine test4
!
  subroutine test5(n, m, s, sym)
    integer, intent(in)   :: n, m, s, sym(n * (s - 1))
    type(mobbrmsd)        :: mobb
    type(mobbrmsd_state)  :: stat
    type(mobbrmsd_input)  :: inp
    real(RK)              :: X(D, n, m), Y(D, n, m)
    real(RK), allocatable :: W(:)
!
    call mobbrmsd_input_add_molecule(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X = sample(n, m)
    Y = sample(n, m)
!
    allocate (W(mobbrmsd_memsize(mobb)))
!
    call mobbrmsd_run(mobb, stat, X, Y, W, difflim=1.0_RK, difflim_absolute=.true.)
    call print_stat(stat, 1.0_RK)
    call mobbrmsd_restart(mobb, stat, W, difflim=0.2_RK, difflim_absolute=.true.)
    call print_stat(stat, 0.2_RK)
    call mobbrmsd_restart(mobb, stat, W, difflim=0.1_RK, difflim_absolute=.true.)
    call print_stat(stat, 0.1_RK)
    call mobbrmsd_restart(mobb, stat, W, difflim=0.01_RK, difflim_absolute=.true.)
    call print_stat(stat, 0.01_RK)
    call mobbrmsd_restart(mobb, stat, W, difflim=0.001_RK, difflim_absolute=.true.)
    call print_stat(stat, 0.001_RK)
    call mobbrmsd_restart(mobb, stat, W)
    call print_stat(stat, 0.0_RK)
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
!
  end subroutine test5
!
  subroutine test6(n, m, s, sym)
    integer, intent(in)   :: n, m, s, sym(n * (s - 1))
    type(mobbrmsd)        :: mobb
    type(mobbrmsd_state)  :: stat
    type(mobbrmsd_input)  :: inp
    real(RK)              :: X(D, n, m), Y(D, n, m)
    real(RK), allocatable :: W(:)
!
    call mobbrmsd_input_add_molecule(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X = sample(n, m)
    Y = sample(n, m)
!
    allocate (W(mobbrmsd_memsize(mobb)))
!
    call mobbrmsd_run(mobb, stat, X, Y, W, maxeval=0)
    call mobbrmsd_restart(mobb, stat, W, cutoff=0.0_RK)
    call print_stat(stat, 0.0_RK)
    call mobbrmsd_restart(mobb, stat, W, cutoff=0.1_RK)
    call print_stat(stat, 0.1_RK)
    call mobbrmsd_restart(mobb, stat, W, cutoff=0.2_RK)
    call print_stat(stat, 0.2_RK)
    call mobbrmsd_restart(mobb, stat, W, cutoff=0.3_RK)
    call print_stat(stat, 0.3_RK)
    call mobbrmsd_run(mobb, stat, X, Y, W, cutoff=0.4_RK)
    call print_stat(stat, 0.4_RK)
    call mobbrmsd_restart(mobb, stat, W)
    call print_stat(stat, 99.0_RK)
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
!
  end subroutine test6
!
  subroutine test7(n, m, s, sym, n_target, n_test)
!$  use omp_lib
    integer, intent(in)    :: n, m, s, sym(n * (s - 1)), n_target, n_test
    type(mobbrmsd)         :: mobb
    type(mobbrmsd_state)   :: state(n_target * (n_target - 1) / 2)
    type(mobbrmsd_input)   :: inp
    real(RK)               :: X(D, n, m, n_target)
    real(RK), allocatable  :: W(:)
    integer(IK)            :: edges(2, n_target - 1)
    integer(IK)            :: redges(2, n_target - 1)
    real(RK)               :: weights(n_target - 1)
    real(RK)               :: refer(n_target, n_target)
    integer(IK)            :: i, j, k, itest, nerr
!
    call mobbrmsd_input_add_molecule(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    !$omp parallel
    if (omp_get_thread_num() == 0) k = MAX(n_target * (n_target - 1) / 2, omp_get_num_threads())
    !$omp end parallel
    allocate (W(mobbrmsd_memsize(mobb) * k))
!
    do itest = 1, n_test
      X(:, :, :, 1) = sample(n, m)
      do i = 2, n_target / 4
        X(:, :, :, i) = 0.9 * X(:, :, :, i - 1) + 0.1 * sample(n, m)
      end do
!
      X(:, :, :, n_target / 4 + 1) = 0.5 * X(:, :, :, n_target / 4) + 0.5 * sample(n, m)
      do i = n_target / 4 + 2, 3 * n_target / 4
        X(:, :, :, i) = 0.9 * X(:, :, :, i - 1) + 0.1 * sample(n, m)
      end do
!
      X(:, :, :, 3 * n_target / 4 + 1) = 0.5 * X(:, :, :, 3 * n_target / 4) + 0.5 * sample(n, m)
      do i = 3 * n_target / 4 + 2, n_target
        X(:, :, :, i) = 0.9 * X(:, :, :, i - 1) + 0.1 * sample(n, m)
      end do
!
      call mobbrmsd_min_span_tree(n_target, mobb, X, W, &
     &                            edges=edges, weights=weights)
!
      call mobbrmsd_batch_tri_run(n_target, mobb, state, X, W)
!
      k = 0
      do j = 1, n_target
        do i = 1, j - 1
          k = k + 1
          refer(i, j) = mobbrmsd_state_rmsd(state(k))
          refer(j, i) = refer(i, j)
        end do
        refer(j, j) = 0.0_RK
      end do
!
      call min_span_tree(n_target, refer, redges)
      call u%assert(is_same_graph(n_target, edges, redges), 'is_same_graph')
      if (.not. is_same_graph(n_target, edges, redges)) then
        nerr = 0
        do i = 1, n_target - 1
          do j = 1, n_target - 1
            if (ALL(edges(:, i) == redges(:, j))) exit
            if (j == n_target - 1) then
              nerr = nerr + 1
              edges(:, nerr) = edges(:, i)
            end if
          end do
        end do
        nerr = 0
        do i = 1, n_target - 1
          do j = 1, n_target - 1
            if (ALL(edges(:, j) == redges(:, i))) exit
            if (j == n_target - 1) then
              nerr = nerr + 1
              redges(:, nerr) = redges(:, i)
            end if
          end do
        end do
        do i = 1, nerr
          print *, edges(:, i), redges(:, i)
        end do
      end if
      FLUSH (OUTPUT_UNIT)
      FLUSH (ERROR_UNIT)
    end do
!
  end subroutine test7
!
  subroutine print_stat(stat, cutoff)
    type(mobbrmsd_state), intent(in) :: stat
    real(RK), intent(in)             :: cutoff
    print'(A,F9.3,I8,*(F9.3))', "# ", &
   &  cutoff, &
   &  mobbrmsd_state_n_eval(stat), &
   &  EXP(mobbrmsd_state_log_eval_ratio(stat)), &
   &  mobbrmsd_state_upperbound(stat), &
   &  mobbrmsd_state_lowerbound(stat), &
   &  mobbrmsd_state_bbgap(stat), &
   &  mobbrmsd_state_rmsd(stat)
  end subroutine print_stat
!
  pure subroutine min_span_tree(n, r, e)
    integer(IK), intent(in)    :: n
    real(RK), intent(in)       :: r(n, n)
    integer(IK), intent(inout) :: e(2, n - 1)
    real(RK)                   :: c
    integer(IK)                :: f(n)
    integer(IK)                :: i, j, k, pi, pj
    do concurrent(i=1:n)
      f(i) = i
    end do
    pi = 0
    pj = 0
    do k = 1, n - 1
      c = 999._RK
      do j = 1, k
        do i = k + 1, n
          if (c > r(f(i), f(j))) then
            c = r(f(i), f(j))
            pi = i
            pj = j
          end if
        end do
      end do
      if (f(pi) < f(pj)) then
        e(:, k) = [f(pi), f(pj)]
      else
        e(:, k) = [f(pj), f(pi)]
      end if
      j = f(k + 1)
      f(k + 1) = f(pi)
      f(pi) = j
    end do
  end subroutine min_span_tree
!
  pure function is_same_graph(n, e1, e2) result(res)
    integer(IK), intent(in) :: n, e1(2, n - 1), e2(2, n - 1)
    logical                 :: res, t
    integer(IK)             :: i, j
    res = .false.
    do i = 1, n - 2
      do j = i + 1, n - 1
        if (ALL(e1(:, i) == e1(:, j))) return
      end do
    end do
    do i = 1, n - 2
      do j = i + 1, n - 1
        if (ALL(e2(:, i) == e2(:, j))) return
      end do
    end do
    do i = 1, n - 1
      do j = 1, n - 1
        t = ALL(e1(:, i) == e2(:, j))
        if (t) exit
      end do
      if (t) cycle
      return
    end do
    res = .true.
  end function is_same_graph
!
end program main


program main
  use mod_dimspec_functions, only: D
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mobbrmsd
  use mod_mobbrmsd_state
  use mod_mobbrmsd_mst
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
#ifdef USE_REAL32
  integer, parameter :: place = 3
#else
  integer, parameter :: place = 7
#endif
!
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
  call u%init('test mobbrmsd cutoff for {(n,M,S)}={(4,4,2)}')
  call test4(4, 4, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd cutoff for {(n,M,S)}={(4,8,1)}')
  call test4(4, 8, 1, [0])
!
  call u%init('test mobbrmsd min_span_tree for {(n,M,S)}={(4,10,1)}, n_target=10')
  call test5(4, 8, 1, [0], 10)
  call u%init('test mobbrmsd min_span_tree for {(n,M,S)}={(4,4,1)}, n_target=100')
  call test5(4, 4, 1, [0], 100)
!
  call u%finish_and_terminate()
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
      call mobbrmsd_run(mobb, stat, X, Y, W)
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
      call mobbrmsd_run(mobb, stat, X, Y, W)
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
      print'(I8, *(f16.9))', mobbrmsd_state_n_eval(stat), &
     &                       EXP(mobbrmsd_state_log_eval_ratio(stat)), &
     &                       mobbrmsd_state_upperbound(stat), &
     &                       mobbrmsd_state_lowerbound(stat)
    end do
!
    call mobbrmsd_restart(mobb, stat, W)
    print'(I8, *(f16.9))', mobbrmsd_state_n_eval(stat), &
   &                       EXP(mobbrmsd_state_log_eval_ratio(stat)), &
   &                       mobbrmsd_state_upperbound(stat), &
   &                       mobbrmsd_state_lowerbound(stat)
    sd = mobbrmsd_state_squared_deviation(stat)
    brute = brute_sd(n, m, s, sym, X, Y)
!
    call mobbrmsd_run(mobb, stat, X, Y)
    sd2 = mobbrmsd_state_squared_deviation(stat)
    print'(I8, *(f16.9))', mobbrmsd_state_n_eval(stat), &
   &                       EXP(mobbrmsd_state_log_eval_ratio(stat)), &
   &                       mobbrmsd_state_upperbound(stat), &
   &                       mobbrmsd_state_lowerbound(stat)
    call u%assert_almost_equal(sd, brute, 'minrmsd value', place=place)
    call u%assert_almost_equal(sd, sd2, 'vs at once   ', place=place)
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
    call mobbrmsd_run(mobb, stat, X, Y, W, maxeval=0)
    call mobbrmsd_restart(mobb, stat, W, cutoff=0.0_RK)
    print'(I8, *(f16.9))', mobbrmsd_state_n_eval(stat), &
   &                       EXP(mobbrmsd_state_log_eval_ratio(stat)), &
   &                       mobbrmsd_state_upperbound(stat), &
   &                       mobbrmsd_state_lowerbound(stat), &
   &                       SQRT((mobbrmsd_state_autovariance(stat) + 2 * mobbrmsd_state_lowerbound(stat)) / (n * m)), &
   &                       mobbrmsd_state_rmsd(stat)
    call mobbrmsd_restart(mobb, stat, W, cutoff=0.1_RK)
    print'(I8, *(f16.9))', mobbrmsd_state_n_eval(stat), &
   &                       EXP(mobbrmsd_state_log_eval_ratio(stat)), &
   &                       mobbrmsd_state_upperbound(stat), &
   &                       mobbrmsd_state_lowerbound(stat), &
   &                       SQRT((mobbrmsd_state_autovariance(stat) + 2 * mobbrmsd_state_lowerbound(stat)) / (n * m)), &
   &                       mobbrmsd_state_rmsd(stat)
    call mobbrmsd_restart(mobb, stat, W, cutoff=0.2_RK)
    print'(I8, *(f16.9))', mobbrmsd_state_n_eval(stat), &
   &                       EXP(mobbrmsd_state_log_eval_ratio(stat)), &
   &                       mobbrmsd_state_upperbound(stat), &
   &                       mobbrmsd_state_lowerbound(stat), &
   &                       SQRT((mobbrmsd_state_autovariance(stat) + 2 * mobbrmsd_state_lowerbound(stat)) / (n * m)), &
   &                       mobbrmsd_state_rmsd(stat)
    call mobbrmsd_restart(mobb, stat, W, cutoff=0.3_RK)
    print'(I8, *(f16.9))', mobbrmsd_state_n_eval(stat), &
   &                       EXP(mobbrmsd_state_log_eval_ratio(stat)), &
   &                       mobbrmsd_state_upperbound(stat), &
   &                       mobbrmsd_state_lowerbound(stat), &
   &                       SQRT((mobbrmsd_state_autovariance(stat) + 2 * mobbrmsd_state_lowerbound(stat)) / (n * m)), &
   &                       mobbrmsd_state_rmsd(stat)
    call mobbrmsd_run(mobb, stat, X, Y, W, cutoff=0.4_RK)
    print'(I8, *(f16.9))', mobbrmsd_state_n_eval(stat), &
   &                       EXP(mobbrmsd_state_log_eval_ratio(stat)), &
   &                       mobbrmsd_state_upperbound(stat), &
   &                       mobbrmsd_state_lowerbound(stat), &
   &                       SQRT((mobbrmsd_state_autovariance(stat) + 2 * mobbrmsd_state_lowerbound(stat)) / (n * m)), &
   &                       mobbrmsd_state_rmsd(stat)
    call mobbrmsd_restart(mobb, stat, W)
    print'(I8, *(f16.9))', mobbrmsd_state_n_eval(stat), &
   &                       EXP(mobbrmsd_state_log_eval_ratio(stat)), &
   &                       mobbrmsd_state_upperbound(stat), &
   &                       mobbrmsd_state_lowerbound(stat), &
   &                       SQRT((mobbrmsd_state_autovariance(stat) + 2 * mobbrmsd_state_lowerbound(stat)) / (n * m)), &
   &                       mobbrmsd_state_rmsd(stat)
!
  end subroutine test4
!
  subroutine test5(n, m, s, sym, n_target)
!$  use omp_lib
    integer, intent(in)    :: n, m, s, sym(n * (s - 1)), n_target
    type(mobbrmsd)         :: mobb
    type(mobbrmsd_state)   :: state(n_target, n_target)
    type(mobbrmsd_input)   :: inp
    real(RK)               :: X(D, n, m, n_target)
    real(RK), allocatable  :: W(:)
    integer(IK)            :: edges(2, n_target - 1)
    real(RK)               :: weights(n_target - 1)
    integer(IK)            :: i, n_job
!
    call mobbrmsd_input_add_molecule(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X(:, :, :, 1) = sample(n, m)
    do i = 2, n_target / 3
      X(:, :, :, i) = 0.9 * X(:, :, :, i - 1) + 0.1 * sample(n, m)
    end do
!
    X(:, :, :, n_target / 3 + 1) = sample(n, m)
    do i = n_target / 3 + 2, 2 * n_target / 3
      X(:, :, :, i) = 0.9 * X(:, :, :, i - 1) + 0.1 * sample(n, m)
    end do
!
    X(:, :, :, 2 * n_target / 3 + 1) = sample(n, m)
    do i = 2 * n_target / 3 + 2, n_target
      X(:, :, :, i) = 0.9 * X(:, :, :, i - 1) + 0.1 * sample(n, m)
    end do
!
    !$omp parallel
    n_job = omp_get_num_threads()
    !$omp end parallel
!
    allocate (W(mobbrmsd_memsize(mobb) * n_job))
!
    call mobbrmsd_min_span_tree(n_target, mobb, state, X, W, &
   &                            edges=edges, weights=weights)
!
    do i = 1, n_target - 1
      print'(2i4, *(f9.3))', edges(:, i), weights(i), &
     &                       mobbrmsd_state_upperbound(state(edges(1, i), edges(2, i))), &
     &                       mobbrmsd_state_lowerbound(state(edges(1, i), edges(2, i))), &
     &                       mobbrmsd_state_upperbound(state(edges(1, i), edges(2, i)))  &
     &                     - mobbrmsd_state_lowerbound(state(edges(1, i), edges(2, i)))
    end do
!
  end subroutine test5
!
end program main


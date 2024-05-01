program main
  use mod_dimspec_functions, only: D
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mobbrmsd
  use mod_mobbrmsd_mst
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
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
  call u%init('test mol_sample')
  call test6(3, 4, 1, [0], 1000)
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(n, m, s, sym)
    integer, intent(in)    :: n, m, s, sym(n * (s - 1))
    type(mobbrmsd)         :: mobb
    type(mol_block_input), allocatable :: inp(:)
    real(RK)               :: X(D, n, m), Y(D, n, m)
    real(RK), allocatable  :: W(:)
    integer(IK)            :: i
!
    call mol_block_input_add(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X = sample(n, m)
    Y = X
!
    allocate (W(mobb%h%memsize()))
!
    do i = 1, 20
      call mobbrmsd_run(mobb%h, mobb%s, X, Y, W)
      call u%assert_almost_equal(mobb%s%squared_deviation(), &
     &                           brute_sd(n, m, s, sym, X, Y), 'minrmsd value')
      Y = 0.5 * Y + 0.5 * sample(n, m)
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
    real(RK)              :: brute
    real(RK)              :: X1(D, n1, m1), X2(D, n2, m2)
    real(RK)              :: Y1(D, n1, m1), Y2(D, n2, m2)
    real(RK)              :: Y(D, m1 * n1 + m2 * n2)
    real(RK), allocatable :: W(:)
    integer(IK)           :: i
!
    call mol_block_input_add(inp, n1, m1, sym=RESHAPE(sym1, [n1, s1 - 1]))
    call mol_block_input_add(inp, n2, m2, sym=RESHAPE(sym2, [n2, s2 - 1]))
    mobb = mobbrmsd(inp)
!
    X1 = sample(n1, m1)
    X2 = sample(n2, m2)
    Y1 = X1
    Y2 = X2
!
    allocate (W(mobb%h%memsize()))
!
    do i = 1, 10
      Y = RESHAPE([Y1, Y2], SHAPE(Y))
      call mobbrmsd_run(mobb%h, mobb%s, [X1, X2], Y, W)
      brute = brute_sd_double(n1, m1, s1, sym1, n2, m2, s2, sym2, X1, Y1, X2, Y2)
      call u%assert_almost_equal(mobb%s%squared_deviation(), brute, 'minrmsd value')
      Y1 = 0.5 * Y1 + 0.5 * sample(n1, m1)
      Y2 = 0.5 * Y2 + 0.5 * sample(n2, m2)
    end do
!
  end subroutine test2
!
  subroutine test3(n, m, s, sym)
    integer, intent(in)    :: n, m, s, sym(n * (s - 1))
    type(mobbrmsd)         :: mobb
    type(mol_block_input), allocatable :: inp(:)
    real(RK)               :: X(D, n, m), Y(D, n, m)
    real(RK), allocatable  :: W(:)
    real(RK)               :: sd, brute, sd2
    integer(IK)            :: i
!
    call mol_block_input_add(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X = sample(n, m)
    Y = sample(n, m)
!
    allocate (W(mobb%h%memsize()))
!
    call mobbrmsd_run(mobb%h, mobb%s, X, Y, W, maxeval=0)
!
    do i = 1, 10
      call mobbrmsd_restart(mobb%h, mobb%s, W, maxeval=0)
      print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), &
     &                       mobb%s%upperbound(), mobb%s%lowerbound()
    end do
!
    call mobbrmsd_restart(mobb%h, mobb%s, W)
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), &
   &                       mobb%s%upperbound(), mobb%s%lowerbound()
    sd = mobb%s%squared_deviation()
    brute = brute_sd(n, m, s, sym, X, Y)
!
    call mobbrmsd_run(mobb%h, mobb%s, X, Y)
    sd2 = mobb%s%squared_deviation()
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), &
   &                       mobb%s%upperbound(), mobb%s%lowerbound()
    call u%assert_almost_equal(sd, brute, 'minrmsd value')
    call u%assert_almost_equal(sd, sd2, 'vs at once   ')
!
    deallocate (inp)
!
  end subroutine test3
!
  subroutine test4(n, m, s, sym)
    integer, intent(in)    :: n, m, s, sym(n * (s - 1))
    type(mobbrmsd)         :: mobb
    type(mol_block_input), allocatable :: inp(:)
    real(RK)               :: X(D, n, m), Y(D, n, m)
    real(RK), allocatable  :: W(:)
!
    call mol_block_input_add(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X = sample(n, m)
    Y = sample(n, m)
!
    allocate (W(mobb%h%memsize()))
!
    call mobbrmsd_run(mobb%h, mobb%s, X, Y, W, maxeval=0)
    call mobbrmsd_restart(mobb%h, mobb%s, W, cutoff=0.0_RK)
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), &
   &                       mobb%s%upperbound(), mobb%s%lowerbound(), &
   &                       SQRT((mobb%s%autovariance() + 2 * mobb%s%lowerbound()) / (n * m)), &
   &                       mobb%s%rmsd()
    call mobbrmsd_restart(mobb%h, mobb%s, W, cutoff=0.1_RK)
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), &
   &                       mobb%s%upperbound(), mobb%s%lowerbound(), &
   &                       SQRT((mobb%s%autovariance() + 2 * mobb%s%lowerbound()) / (n * m)), &
   &                       mobb%s%rmsd()
    call mobbrmsd_restart(mobb%h, mobb%s, W, cutoff=0.2_RK)
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), &
   &                       mobb%s%upperbound(), mobb%s%lowerbound(), &
   &                       SQRT((mobb%s%autovariance() + 2 * mobb%s%lowerbound()) / (n * m)), &
   &                       mobb%s%rmsd()
    call mobbrmsd_restart(mobb%h, mobb%s, W, cutoff=0.3_RK)
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), &
   &                       mobb%s%upperbound(), mobb%s%lowerbound(), &
   &                       SQRT((mobb%s%autovariance() + 2 * mobb%s%lowerbound()) / (n * m)), &
   &                       mobb%s%rmsd()
    call mobbrmsd_run(mobb%h, mobb%s, X, Y, W, cutoff=0.4_RK)
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), &
   &                       mobb%s%upperbound(), mobb%s%lowerbound(), &
   &                       SQRT((mobb%s%autovariance() + 2 * mobb%s%lowerbound()) / (n * m)), &
   &                       mobb%s%rmsd()
    call mobbrmsd_restart(mobb%h, mobb%s, W)
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), &
   &                       mobb%s%upperbound(), mobb%s%lowerbound(), &
   &                       SQRT((mobb%s%autovariance() + 2 * mobb%s%lowerbound()) / (n * m)), &
   &                       mobb%s%rmsd()
!
    deallocate (inp)
!
  end subroutine test4
!
  subroutine test5(n, m, s, sym, n_target)
!$  use omp_lib
    integer, intent(in)    :: n, m, s, sym(n * (s - 1)), n_target
    type(mobbrmsd)         :: mobb
    type(mobbrmsd_state)   :: state(n_target, n_target)
    type(mol_block_input), allocatable :: inp(:)
    real(RK)               :: X(D, n, m, n_target)
    real(RK), allocatable  :: W(:)
    integer(IK)            :: edges(2, n_target - 1)
    real(RK)               :: weights(n_target - 1)
    integer(IK)            :: i, n_job
!
    call mol_block_input_add(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
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
    allocate (W(mobb%h%memsize() * n_job))
!
    call mobbrmsd_min_span_tree(n_target, mobb%h, state, X, W, &
   &                            edges=edges, weights=weights, &
   &                            verbose=.true.)
!
    do i = 1, n_target - 1
      print'(2i4, *(f9.3))', edges(:, i), weights(i), &
     &                        state(edges(1, i), edges(2, i))%upperbound(), &
     &                        state(edges(1, i), edges(2, i))%lowerbound(), &
     &                        state(edges(1, i), edges(2, i))%upperbound()  &
     &                      - state(edges(1, i), edges(2, i))%lowerbound()
    end do
!
    deallocate (inp)
!
  end subroutine test5
!
  subroutine test6(n_apm, n_mol, n_sym, sym, n_target)
    integer, intent(in)    :: n_apm, n_mol, n_sym, sym(n_apm * (n_sym - 1)), n_target
    type(mobbrmsd)         :: mobb
    type(mol_block_input), allocatable :: inp(:)
    real(RK)               :: XY(D, n_apm, n_mol, 2)
    real(RK)               :: a, b, mean, var
    real(RK), allocatable  :: W(:)
    integer(IK)            :: i, j, k
!
    call mol_block_input_add(inp, n_apm, n_mol, sym=RESHAPE(sym, [n_apm, n_sym - 1]))
    mobb = mobbrmsd(inp)
    allocate (W(mobb%h%memsize()))
    b = 0.0_RK
    do k = 1, 3
      a = 0.0_RK
      do j = 1, 10
        mean = ZERO
        var = ZERO
        do i = 1, n_target
          XY = RESHAPE(mol_sample(n_apm, n_mol * 2, a, b), SHAPE(XY))
          call mobbrmsd_run(mobb%h, mobb%s, XY(1, 1, 1, 1), XY(1, 1, 1, 2), W)
          mean = mean + mobb%s%n_eval()
          var = var + mobb%s%n_eval()**2
        end do
        a = a + 0.1_RK
        print '(2f6.1, 2f16.9)', a, b, mean / n_target, SQRT(var / n_target - mean**2 / (n_target**2))
      end do
      print *
      b = b + 0.5_RK
    end do
!
    deallocate (inp)
!
  end subroutine test6
!
end program main


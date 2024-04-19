program main
  use blas_lapack_interface, only: D
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mobbrmsd
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test mobbrmsd for (n,M,S)=(1,1,1)')
  call test1(1, 1, 1, [0])
! call u%init('test mobbrmsd for (n,M,S)=(1,2,1)')
! call test1(1, 2, 1, [0])
! call u%init('test mobbrmsd for (n,M,S)=(1,3,1)')
! call test1(1, 3, 1, [0])
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
  call u%init('test mobbrmsd repeat for {(n,M,S)}={(4,3,1)}')
  call test3(4, 8, 1, [0])
!
  call u%init('test mobbrmsd min_span_tree for {(n,M,S)}={(4,4,1)}, n_target=10')
  call test4(4, 4, 1, [0], 10)
  call u%init('test mobbrmsd min_span_tree for {(n,M,S)}={(4,6,1)}, n_target=10')
  call test4(4, 10, 1, [0], 10)
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
      call u%assert_almost_equal(w(1), brute_sd(n, m, s, sym, X, Y), 'minrmsd value')
      Y = 0.8 * Y + 0.2 * sample(n, m)
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
      call u%assert_almost_equal(w(1), brute, 'minrmsd value')
      Y1 = 0.8 * Y1 + 0.2 * sample(n1, m1)
      Y2 = 0.8 * Y2 + 0.2 * sample(n2, m2)
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
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), mobb%s%upperbound(), mobb%s%lowerbound()
!
    do i = 1, 10
      call mobbrmsd_restart(mobb%h, mobb%s, W, maxeval=i * 500)
      print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), mobb%s%upperbound(), mobb%s%lowerbound()
    end do
!
    call mobbrmsd_restart(mobb%h, mobb%s, W)
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), mobb%s%upperbound(), mobb%s%lowerbound()
    call u%assert_almost_equal(mobb%s%squared_deviation(), brute_sd(n, m, s, sym, X, Y), 'minrmsd value')
!
    call mobbrmsd_run(mobb%h, mobb%s, X, Y, W)
    print'(I8, *(f16.9))', mobb%s%n_eval(), EXP(mobb%s%log_eval_ratio()), mobb%s%upperbound(), mobb%s%lowerbound()
!
    deallocate (inp)
!
  end subroutine test3
!
  subroutine test4(n, m, s, sym, n_target)
    integer, intent(in)    :: n, m, s, sym(n * (s - 1)), n_target
    type(mobbrmsd)         :: mobb
    type(mobbrmsd_state)   :: state(n_target, n_target)
    type(mol_block_input), allocatable :: inp(:)
    real(RK)               :: X(D, n, m, n_target)
    real(RK), allocatable  :: W(:)
    integer(IK)            :: edges(2, n_target - 1)
    real(RK)               :: weights(n_target - 1)
    integer(IK)            :: i
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
    allocate (W(mobb%h%memsize() * mobbrmsd_num_threads()))
!
    call mobbrmsd_min_span_tree(n_target, mobb%h, state, X, W, &
   &                            edges=edges, weights=weights)
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
  end subroutine test4
!
end program main


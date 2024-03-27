program main
  use mod_params, only: D, setup_dimension, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_rotation
  use mod_bb_block
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call setup_dimension(2)
  call u%init('test bb_block d=2')
  call test0()
  call test1(8, 4, 2, [3,4,2,1], [1,2,1,2], [2, 3, 5, 4, 8, 6, 7, 1])
  call test1(8, 4, 1, [3,4,2,1], [1,1,1,1], [0])
!
  call setup_dimension(3)
  call u%init('test bb_block d=3')
  call test0()
  call test1(8, 4, 2, [3,4,2,1], [1,2,1,2], [2, 3, 5, 4, 8, 6, 7, 1])
  call test1(8, 4, 1, [3,4,2,1], [1,1,1,1], [0])
!
  call setup_dimension(4)
  call u%init('test bb_block d=4')
  call test0()
  call test1(8, 4, 2, [3,4,2,1], [1,2,1,2], [2, 3, 5, 4, 8, 6, 7, 1])
  call test1(8, 4, 1, [3,4,2,1], [1,1,1,1], [0])
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    type(bb_block)         :: bm
    integer(IK), parameter :: m = 8
    integer(IK), parameter :: n = 3
    integer(IK), parameter :: s = 2
    integer(IK), parameter :: sym(m) = [2, 3, 4, 5, 6, 7, 8, 1]
    real(RK)               :: X(D, m, n), Y(D, m, n), Z(D, m, n)
    real(RK), allocatable  :: W(:)
    integer(IK)            :: sb(n + 1)
    integer(IK)            :: i
!
    bm = bb_block(8, 3, sym=RESHAPE(sym, [m, 1]))
!
    allocate (W(bb_block_memsize(bm%q) + bb_block_worksize(bm%q)))
!
    X = RESHAPE(sample(D, m * n), SHAPE(X))
    Y(:, :, :2) = X(:, :, 2:)
    Y(:, :, 3) = X(:, :, 1)
!
    do i = 1, 20
      call run(bm%q, bm%s, X, Y, W, sb)
      call u%assert_almost_equal(W(1), brute_sd(m, n, s, sym, X, Y), 'minrmsd value')
      Z = Y
      call bb_block_swap_y(bm%q, sb, Z)
      call u%assert_almost_equal(W(1), sd(m, n, X, Z),               'swap sd value')
      Y = 0.8 * Y + RESHAPE(sample(D, m * n), SHAPE(Y)) * 0.2
    end do
!
  end subroutine test0
!
  subroutine test1(m, n, s, per, map, sym)
    integer(IK), intent(in) :: m, n, s, per(n), map(n), sym(m * (s - 1))
    type(bb_block)          :: bm
    real(RK)                :: X(D, m, n), Y(D, m, n), Z(D, m, n)
    real(RK), allocatable   :: W(:)
    integer(IK)             :: sb(n + 1)
    integer(IK)             :: i
!
    bm = bb_block(m, n, RESHAPE(sym, [m, s - 1]))
    allocate (W(bb_block_memsize(bm%q) + bb_block_worksize(bm%q)))
!
    X = RESHAPE(sample(D, m * n), SHAPE(X))
    Y = swp(m, n, s, per, map, sym, X)
    do i = 1, 20
      call run(bm%q, bm%s, X, Y, W, sb)
      call u%assert_almost_equal(W(1), brute_sd(m, n, s, sym, X, Y), 'minrmsd value')
      Z = Y
      call bb_block_swap_y(bm%q, sb, Z)
      call u%assert_almost_equal(W(1), sd(m, n, X, Z),               'swap sd value')
      Y = 0.8 * Y + RESHAPE(sample(D, m * n), SHAPE(Y)) * 0.2
    end do
!
  end subroutine test1
!
  subroutine run(q, s, X, Y, W, sb)
    integer(IK), intent(in)    :: q(*)
    integer(IK), intent(inout) :: s(*)
    real(RK), intent(in)       :: X(*), Y(*)
    real(RK), intent(inout)    :: W(:)
    integer(IK), intent(inout) :: sb(*)
    integer(IK)                :: nmol
    real(RK)                   :: ub
!
    ub = 999.9_RK
    nmol = bb_block_nmol(q)
    call bb_block_setup(q, X, Y, s, W, zfill=.TRUE.)
!
    do
      call bb_block_expand(ub, q, s, W)
!
      if (.not. bb_block_queue_is_empty(q, s) &
        & .and. bb_block_queue_is_bottom(q, s)) then
        if (bb_block_current_value(q, s, w) < ub) sb(:1 + nmol) = s(:1 + nmol)
        ub = MIN(bb_block_current_value(q, s, w), ub)
      end if
!
      call bb_block_leave(ub, q, s, W)
      if(bb_block_queue_is_empty(q, s)) exit
    end do
!
    W(1) = ub
!
  end subroutine run
!
end program main

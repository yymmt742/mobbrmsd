program main
  use mod_dimspec_functions, only: D
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_rotation
  use mod_bb_block
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
  call u%init('test bb_block [1,1,1]')
  call test1(1, 1, 1, [0])
  call u%init('test bb_block [1,2,1]')
  call test1(1, 2, 1, [0])
  call u%init('test bb_block [8,4,1]')
  call test1(8, 4, 1, [0])
  call u%init('test bb_block [8,4,2]')
  call test1(8, 4, 2, [2, 3, 5, 4, 8, 6, 7, 1])
!
  call u%finish_and_terminate()
!
contains
  subroutine test1(n, m, s, sym)
    integer(IK), intent(in) :: n, m, s, sym(n * (s - 1))
    type(bb_block)          :: bm
    real(RK)                :: X(D, n, m), Y(D, n, m), Z(D, n, m)
    real(RK)                :: CX(D), CY(D)
    real(RK), allocatable   :: W(:)
    integer(IK)             :: sb(m)
    integer(IK)             :: i
!
    bm = bb_block(n, m, RESHAPE(sym, [n, s - 1]))
    allocate (W(bb_block_memsize(bm%q) + bb_block_worksize(bm%q)))
!
    X = sample(n, m)
    Y = X
    CX = SUM(RESHAPE(X, [D, m * n]), 2) / (m * n)
    do i = 1, 20
      CY = SUM(RESHAPE(Y, [D, m * n]), 2) / (m * n)
      call run(bm%q, bm%s, X, Y, CX, CY, W, sb)
      call u%assert_almost_equal(W(1), brute_sd(n, m, s, sym, X, Y), 'minrmsd value', place=place)
      call u%assert_almost_equal(W(2), autovar(n * m, X, Y), 'auto variance', place=place)
      Z = Y
      call bb_block_swap_y(bm%q, bm%s, sb, Z)
      call u%assert_almost_equal(W(1), sd(m * n, X, Z), 'swap sd value', place=place)
      Y = 0.5 * Y + sample(n, m) * 0.5
    end do
  end subroutine test1
!
  subroutine run(q, s, X, Y, CX, CY, W, sb)
    integer(IK), intent(in)    :: q(*)
    integer(IK), intent(inout) :: s(*)
    real(RK), intent(in)       :: X(*), Y(*), CX(*), CY(*)
    real(RK), intent(inout)    :: W(*)
    integer(IK), intent(inout) :: sb(*)
    real(RK)                   :: ub, g
    ub = 999.9_RK
    call bb_block_setup(q, X, Y, CX, CY, s, W, zfill=.true.)
    g = bb_block_autocorr(q, W)
    do
      call bb_block_expand(ub, q, s, W)
      if (bb_block_is_bottom(q, s)) then
        if (bb_block_current_value(q, s, w) < ub) then
          ub = bb_block_current_value(q, s, w)
          call bb_block_save_state(q, s, sb)
        end if
      end if
      call bb_block_closure(ub, q, s, W)
      if (bb_block_tree_is_empty(q, s)) exit
    end do
    W(1) = g + ub + ub
    W(2) = g
  end subroutine run
end program main


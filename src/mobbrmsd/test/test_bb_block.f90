program main
  use blas_lapack_interface, only: D
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_rotation
  use mod_bb_block
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test bb_block')
  call test0()
  call u%init('test bb_block manual')
  call test1(1, 1, 1, [1], [1], [0])
  call test1(1, 2, 1, [1, 2], [1, 1], [0])
  call test1(8, 4, 1, [3, 4, 2, 1], [1, 1, 1, 1], [0])
  call test1(8, 4, 2, [3, 4, 2, 1], [1, 2, 1, 2], [2, 3, 5, 4, 8, 6, 7, 1])
! call test1(8, 8, 1, [3, 4, 2, 1, 5, 7, 6, 8], [1, 1, 1, 1, 1, 1, 1, 1], [0])
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    type(bb_block)         :: bm
    integer(IK), parameter :: m = 3
    integer(IK), parameter :: n = 8
    integer(IK), parameter :: s = 2
    integer(IK), parameter :: sym(n) = [2, 3, 4, 5, 6, 7, 8, 1]
    real(RK)               :: X(D, n, m), Y(D, n, m), Z(D, n, m)
    real(RK)               :: CX(D), CY(D)
    real(RK)               :: G, C(D, D), R(D, D), V(rotation_worksize()), sxz, rxz
    real(RK), allocatable  :: W(:)
    integer(IK)            :: sb(m)
    integer(IK)            :: i
!
    bm = bb_block(n, m, sym=RESHAPE(sym, [n, 1]))
!
    allocate (W(bb_block_memsize(bm%q) + bb_block_worksize(bm%q)))
!
    X = sample(n, m)
    Y(:, :, :2) = X(:, :, 2:)
    Y(:, :, 3) = X(:, :, 1)
    CX = SUM(RESHAPE(X, [D, n * m]), 2) / (n * m)
!
    do i = 1, 20
      CY = SUM(RESHAPE(Y, [D, n * m]), 2) / (n * m)
      call run(bm%q, bm%s, X, Y, CX, CY, W, sb)
      call u%assert_almost_equal(W(1), brute_sd(n, m, s, sym, X, Y), 'minrmsd value')
      call u%assert_almost_equal(W(2), SUM(X * X + Y * Y), 'autocorr     ')
      Z = Y
      call bb_block_swap_y(bm%q, sb, Z)
      sxz = sd(m * n, X, Z)
      call u%assert_almost_equal(W(1), sxz, 'swap sd value')
!
      G = ZERO
      C = ZERO
      call bb_block_covmat_add(bm%q, sb, W, G, C)
      call estimate_rotation(G, C, R, V)
      rxz = SUM((X - RESHAPE(MATMUL(TRANSPOSE(R), RESHAPE(Z, [D, m * n])), [D, n, m]))**2)
!
      call u%assert_almost_equal(sxz, rxz, 'rotmat value ')
!
      Y = 0.5 * Y + sample(n, m) * 0.5
    end do
!
  end subroutine test0
!
  subroutine test1(n, m, s, per, map, sym)
    integer(IK), intent(in) :: n, m, s, per(m), map(m), sym(n * (s - 1))
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
    Y = swp(n, m, s, per, map, sym, X)
    CX = SUM(RESHAPE(X, [D, m * n]), 2) / (m * n)
    do i = 1, 20
      CY = SUM(RESHAPE(Y, [D, m * n]), 2) / (m * n)
      call run(bm%q, bm%s, X, Y, CX, CY, W, sb)
      call u%assert_almost_equal(W(1), brute_sd(n, m, s, sym, X, Y), 'minrmsd value')
      call u%assert_almost_equal(W(2), SUM(X * X + Y * Y), 'autocorr     ')
      Z = Y
      call bb_block_swap_y(bm%q, sb, Z)
      call u%assert_almost_equal(W(1), sd(m * n, X, Z), 'swap sd value')
      Y = 0.5 * Y + sample(n, m) * 0.5
    end do
!
  end subroutine test1
!
  subroutine run(q, s, X, Y, CX, CY, W, sb)
    integer(IK), intent(in)    :: q(*)
    integer(IK), intent(inout) :: s(*)
    real(RK), intent(in)       :: X(*), Y(*), CX(*), CY(*)
    real(RK), intent(inout)    :: W(*)
    integer(IK), intent(inout) :: sb(*)
    real(RK)                   :: ub, g
!
    ub = 999.9_RK
    call bb_block_setup(q, X, Y, CX, CY, s, W, zfill=.true.)
    g = bb_block_autocorr(q, W)
!
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
!
    W(1) = g + ub + ub
    W(2) = g
!
  end subroutine run
!
end program main


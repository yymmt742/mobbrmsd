program main
  use mod_params, only: D, setup_dimension, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_rotation_matrix
  use mod_bb_block
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test bb_block')
!
  call setup_dimension(3)
  call test0()
! call test1()
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
    real(RK)               :: X(D, m, n), Y(D, m, n)
    real(RK), allocatable  :: W(:)
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
      call run(bm%q, bm%s, X, Y, W)
      print'(2f9.3)', W(1), brute_sd(m, n, s, sym, X, Y)
      print*
      Y = 0.8 * Y + RESHAPE(sample(D, m * n), SHAPE(Y)) * 0.2
    end do
!
  end subroutine test0
!
  subroutine test1()
    type(bb_block)         :: bm
    integer(IK), parameter :: m = 8
    integer(IK), parameter :: n = 3
    real(RK)               :: X(D, m * n), Y(D, m * n)
    real(RK), allocatable  :: W(:)
    integer(IK)            :: i
!
    bm = bb_block(8, 3)
    print'(4I4)',bm%q
    print*
    print'(4I4)',bm%s
!
    allocate (W(bb_block_memsize(bm%q) + bb_block_worksize(bm%q)))
!
    X = sample(D, m * n)
    Y(:, m + 1:m * n) = X(:, :m * (n - 1))
    Y(:, :m) = X(:, m * (n - 1) + 1:m * n)
    do i = 1, 100
      call run(bm%q, bm%s, X, Y, W)
      Y = 0.8 * Y + sample(D, m * n) * 0.2
    end do
!
  end subroutine test1
!
  subroutine run(q, s, X, Y, W)
    integer(IK), intent(in)    :: q(*)
    integer(IK), intent(inout) :: s(*)
    real(RK), intent(in)       :: X(*), Y(*)
    real(RK), intent(inout)    :: W(:)
    real(RK)                   :: ub
!
    ub = 999.9_RK
    call bb_block_setup(q, X, Y, s, W, zfill=.TRUE.)
!
    do
      call bb_block_expand(ub, q, s, W)
      if(bb_block_queue_is_empty(q, s)) exit
      print '(2L4,2f9.3,*(I3))', bb_block_queue_is_empty(q, s), &
        &                       bb_block_queue_is_bottom(q, s), &
        &                       bb_block_current_value(q, s, w), &
        &                       bb_block_lowest_value(q, s, w), &
        &                       s(:4)
      ub = bb_block_current_value(q, s, w)
      call bb_block_leave(ub, q, s, W)
    end do
!
    W(1) = ub
!
  end subroutine run
!
! subroutine test3()
!   integer, parameter    :: s = 1
!   integer, parameter    :: m = 5, n = 5, g = 3
!   integer, parameter    :: mn = m * n
!   type(mol_block)       :: b = mol_block(0, 2, m, n, g)
!   type(s_matrix_list)   :: dm
!   type(mol_block_list)  :: blk
!   type(mol_symmetry)    :: ms(s)
!   real(RK)              :: X(d, mn), Y(d, mn)
!   real(RK)              :: R1(6, 2, 2, 2), R2(6, 2, 2, 2)
!   real(RK), allocatable :: W(:)
!   integer               :: i, j, k
!
!   ms(1) = mol_symmetry(RESHAPE([2, 1, 3, 4, 5], [m, 1]))
!   blk = mol_block_list(s, [b])
!   dm = s_matrix_list(blk, 1)
!   allocate (w(dm%memsize()))
!   W(:)=999
!
!   X = sample(d, mn)
!   Y = sample(d, mn)
!
!   call dm%eval(ms, X, Y, W)
!
!   do k = 1, 2
!   do j = 1, 2
!   do i = 1, 2
!     R1(1, i, j, k) = sd(d, X, swp(d, m, n, [1, 2, 3], [i, j, k] - 1, ms(1), Y))
!     R1(2, i, j, k) = sd(d, X, swp(d, m, n, [1, 3, 2], [i, j, k] - 1, ms(1), Y))
!     R1(3, i, j, k) = sd(d, X, swp(d, m, n, [2, 1, 3], [i, j, k] - 1, ms(1), Y))
!     R1(4, i, j, k) = sd(d, X, swp(d, m, n, [2, 3, 1], [i, j, k] - 1, ms(1), Y))
!     R1(5, i, j, k) = sd(d, X, swp(d, m, n, [3, 1, 2], [i, j, k] - 1, ms(1), Y))
!     R1(6, i, j, k) = sd(d, X, swp(d, m, n, [3, 2, 1], [i, j, k] - 1, ms(1), Y))
!     R2(1, i, j, k) = pe(dm, d, 3, [1, 2, 3, 1, 2, 3, 1, 2, 3], [1, 1, 1], [i, j, k] - 1, W)
!     R2(2, i, j, k) = pe(dm, d, 3, [1, 2, 3, 1, 2, 3, 1, 3, 2], [1, 2, 1], [i, j, k] - 1, W)
!     R2(3, i, j, k) = pe(dm, d, 3, [1, 2, 3, 2, 1, 3, 2, 1, 3], [2, 1, 1], [i, j, k] - 1, W)
!     R2(4, i, j, k) = pe(dm, d, 3, [1, 2, 3, 2, 1, 3, 2, 3, 1], [2, 2, 1], [i, j, k] - 1, W)
!     R2(5, i, j, k) = pe(dm, d, 3, [1, 2, 3, 3, 1, 2, 3, 1, 2], [3, 1, 1], [i, j, k] - 1, W)
!     R2(6, i, j, k) = pe(dm, d, 3, [1, 2, 3, 3, 1, 2, 3, 2, 1], [3, 2, 1], [i, j, k] - 1, W)
!   end do
!   end do
!   end do
!   call u%assert_almost_equal([R1 - R2], ZERO, 'R1 - R2')
!
! end subroutine test3
!
! function pe(dm, d, l, iper, jper, isym, W) result(res)
!   type(s_matrix_list), intent(in) :: dm
!   integer, intent(in)             :: d, l, iper(l, l), jper(l), isym(l)
!   real(RK), intent(in)            :: W(*)
!   real(RK)                        :: C(d,d), H, T, LF, LB, res
!   integer                         :: i
!   H = W(dm%h)
!   call copy(d * d, W(dm%c), C)
!   do i = 1, l
!     call dm%partial_eval(i, iper(1, i), jper(i), isym(i), W, T, H, C, LF=LF, LB=LB)
!   end do
!   res = LF
!   !res = T
! end function pe
!
end program main

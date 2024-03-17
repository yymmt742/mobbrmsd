program main
  use mod_params, only: D, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_bb_block
  use mod_branch_and_bound
  use mod_unittest
  use mod_testutil
  implicit none
  type(unittest) :: u
! integer, parameter :: NTEST=25
! integer            :: itest
!
  call u%init('test branch_and_bound')
!
  call test0()
!
! do itest = 1, NTEST
!   call test1()
! end do
!
! call test2()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    type(bb_block)         :: blk(2)
    type(branch_and_bound) :: b
    integer(IK), parameter :: m = 8
    integer(IK), parameter :: n = 3
    real(RK)               :: X(D, m * n), Y(D, m * n)
!
    blk(1) = bb_block(8, 3, sym=RESHAPE([2, 3, 4, 5, 6, 7, 8, 1], [8, 1]))
    blk(2) = bb_block(4, 5, sym=RESHAPE([1, 3, 2, 4], [4, 1]))
    b = branch_and_bound(blk)

    X = sample(D, m * n)
    Y(:, m + 1:m * n) = X(:, :m * (n - 1))
    Y(:, :m) = X(:, m * (n - 1) + 1:m * n)

    print'(4i4)',blk(2)%q
    print*
    print'(4i4)',b%q
    print*
    print'(4i4)',b%s
    print*
    print*, branch_and_bound_memsize(b%q)
    print*, branch_and_bound_worksize(b%q)
!   call bb_block_setup(bm%q, X, Y, bm%x)
!
!   call bb_block_expand(bm%q, bm%x, bm%w)
!   call bb_block_select_top_node(bm%q, bm%x, 999.0_RK)
!   print *, bb_block_queue_is_empty(bm%q), &
!     &      bb_block_queue_is_bottom(bm%q), &
!     &      bb_block_current_value(bm%q, bm%x)
!   call bb_block_expand(bm%q, bm%x, bm%w)
!   call bb_block_select_top_node(bm%q, bm%x, 999.0_RK)
!   print *, bb_block_queue_is_empty(bm%q), &
!     &      bb_block_queue_is_bottom(bm%q), &
!     &      bb_block_current_value(bm%q, bm%x)
!   call bb_block_expand(bm%q, bm%x, bm%w)
!   call bb_block_select_top_node(bm%q, bm%x, 999.0_RK)
!   ub = bb_block_current_value(bm%q, bm%x)
!   print *, bb_block_queue_is_empty(bm%q), &
!     &      bb_block_queue_is_bottom(bm%q), &
!     &      bb_block_current_value(bm%q, bm%x)
!
  end subroutine test0
!
! subroutine test1()
!   integer, parameter     :: l = 1
!   integer, parameter     :: s = 2
!   integer, parameter     :: m = 5, n = 8, g = 3
!   integer, parameter     :: mn = m * n
!   type(mol_block)        :: b = mol_block(0, s, m, n, g)
!   type(branch_and_bound) :: bra
!   type(mol_block_list)   :: blk
!   type(mol_symmetry)     :: ms(l)
!   real(RK)               :: X(d, mn), Y(d, mn), isd, msd
!   real(RK), allocatable  :: W(:)
!   integer                :: i, j, k
!
!   ms(1) = mol_symmetry(RESHAPE([2, 3, 1, 4, 5], [m, 1]))
!   blk = mol_block_list(l, [b])
!
!   X = sample(d, mn)
!   Y = sample(d, mn)
!
!   bra = branch_and_bound(blk, ms)
!
!   allocate (W(bra%memsize))
!   call bra%setup(X, Y, W)
!   call bra%run(W, .true.)
!   print'(*(f16.3))', w(bra%lncmb), w(bra%nsrch), EXP(w(bra%ratio))
!
!   msd = 999D0
!   do k=0,s-1
!   do j=0,s-1
!   do i=0,s-1
!     isd = sd(X, swp(m, n, [1, 2, 3], [i, j, k], ms(1), Y)) ; msd = MIN(msd, isd)
!     isd = sd(X, swp(m, n, [1, 3, 2], [i, j, k], ms(1), Y)) ; msd = MIN(msd, isd)
!     isd = sd(X, swp(m, n, [2, 1, 3], [i, j, k], ms(1), Y)) ; msd = MIN(msd, isd)
!     isd = sd(X, swp(m, n, [2, 3, 1], [i, j, k], ms(1), Y)) ; msd = MIN(msd, isd)
!     isd = sd(X, swp(m, n, [3, 1, 2], [i, j, k], ms(1), Y)) ; msd = MIN(msd, isd)
!     isd = sd(X, swp(m, n, [3, 2, 1], [i, j, k], ms(1), Y)) ; msd = MIN(msd, isd)
!   enddo
!   enddo
!   enddo
!
!   Y = RESHAPE(W(bra%yp:bra%yp + d * mn), [d, mn])
!
!   call u%assert_almost_equal(msd, W(bra%upperbound),             'branchcut vs brute')
!   call u%assert_almost_equal(SUM((X - Y)**2), W(bra%upperbound), 'swap a            ')
!   call u%assert_almost_equal(sd(X, Y), W(bra%upperbound),     'swap b            ')
!
! end subroutine test1
!
! subroutine test2()
!   integer, parameter     :: s = 3
!   integer, parameter     :: m1 = 5, n1 = 3, g1 = 2
!   integer, parameter     :: m2 = 3, n2 = 4, g2 = 4
!   integer, parameter     :: m3 = 7, n3 = 6, g3 = 5
!   integer, parameter     :: mn = m1 * n1 + m2 * n2 + m3 * n3
!   type(mol_block)        :: b(3) = [mol_block(0, 3, m1, n1, g1), &
!                                  &  mol_block(0, 1, m2, n2, g2), &
!                                  &  mol_block(0, 2, m3, n3, g3)]
!   type(branch_and_bound) :: bra
!   type(mol_block_list)   :: blk
!   type(mol_symmetry)     :: ms(s)
!   real(RK)               :: X(d, mn), Y(d, mn)
!   real(RK), allocatable  :: W(:)
!   integer                :: i
!
!   ms(1) = mol_symmetry(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m1, 2]))
!   ms(2) = mol_symmetry(RESHAPE([(i, i=1,0)], [0, 1]))
!   ms(3) = mol_symmetry(RESHAPE([7, 6, 5, 4, 3, 2, 1], [m3, 1]))
!   blk = mol_block_list(s, b)
!
!   X = sample(D, mn)
!   Y = sample(D, mn)
!
!   bra = branch_and_bound(blk, ms)
!   allocate (W(bra%memsize))
!   call bra%setup(X, Y, W)
!   call bra%run(W, .true.)
!   print'(*(f16.3))', w(bra%lncmb), w(bra%nsrch), EXP(w(bra%ratio))
!   Y = RESHAPE(W(bra%yp:bra%yp + d * mn), [d, mn])
!   call u%assert_almost_equal(sd(X, Y), W(bra%upperbound), 'multiple swap')
!
! end subroutine test2
!
! pure function swp(m, n, per, sym, ms, X) result(res)
!   integer(IK), intent(in)        :: m, n, per(:), sym(:)
!   type(mol_symmetry), intent(in) :: ms
!   real(RK), intent(in)           :: X(d, m, n)
!   real(RK)                       :: tmp(d, m, n), res(d, m * n)
!   integer(IK)                    :: i
!   tmp = X
!   do i = 1, SIZE(per)
!     tmp(:, :, per(i)) = X(:, :, i)
!     call ms%swap(D, tmp(:, :, per(i)), sym(i))
!   end do
!   res = RESHAPE(tmp, [D, m * n])
! end function swp
!
! pure function sd(X, Y) result(res)
!   real(RK), intent(in)    :: X(:, :), Y(:, :)
!   real(RK)                :: C(D, D), R(D, D), W(100), res
!   C = MATMUL(Y, TRANSPOSE(X))
!   call estimate_rotation_matrix(SUM(X * X) + SUM(Y * Y), C, R, W)
!   res = SUM(X**2) + SUM(Y**2) - 2 * SUM(C * R)
! end function sd
!
! function sample(d, n) result(res)
!   integer, intent(in)  :: d, n
!   real(RK)             :: cnt(d)
!   real(RK)             :: res(d, n)
!   integer              :: i
!   call RANDOM_NUMBER(res)
!   cnt = SUM(res, 2) / n
!   do concurrent(i=1:n)
!     res(:, i) = res(:, i) - cnt
!   enddo
! end function sample
!
end program main

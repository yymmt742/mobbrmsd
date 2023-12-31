program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_estimate_rotation_matrix
  use mod_molecular_rotation
  use mod_branch_and_prune
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer, parameter :: NTEST=25
  integer            :: itest
!
  call u%init('test branch_and_prune')
  do itest = 1, NTEST
    call test1()
  end do
  do itest = 1, NTEST
    call test2()
  end do
!
    call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter          :: d = 3
    integer, parameter          :: l = 1
    integer, parameter          :: s = 2
    integer, parameter          :: m = 5, n = 8, f = 5, g = 3
    integer, parameter          :: mn = m * n
    type(mol_block)             :: b = mol_block(0, s, m, n, f, g)
    type(branch_and_prune)      :: bra
    type(mol_block_list)        :: blk
    type(molecular_rotation)    :: rot(l)
    real(RK)                    :: X(d, mn), Y(d, mn), Z(d, mn), isd, msd
    real(RK), allocatable       :: W(:)
    integer                     :: i, j, k
!
    rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5], [m, 1]))
    blk = mol_block_list(d, l, [b])
!
    X = sample(d, mn)
    Y = sample(d, mn)
!
    bra = branch_and_prune(blk, 1, rot)
    allocate (W(bra%memsize()))
    call bra%setup(X, Y, W)
    call bra%run(W)
!
    msd = 999D0
    do k=0,s-1
    do j=0,s-1
    do i=0,s-1
      isd = sd(d, X, swp(d, m, n, [1, 2, 3], [i, j, k], rot(1), Y)) ; msd = MIN(msd, isd)
      isd = sd(d, X, swp(d, m, n, [1, 3, 2], [i, j, k], rot(1), Y)) ; msd = MIN(msd, isd)
      isd = sd(d, X, swp(d, m, n, [2, 1, 3], [i, j, k], rot(1), Y)) ; msd = MIN(msd, isd)
      isd = sd(d, X, swp(d, m, n, [2, 3, 1], [i, j, k], rot(1), Y)) ; msd = MIN(msd, isd)
      isd = sd(d, X, swp(d, m, n, [3, 1, 2], [i, j, k], rot(1), Y)) ; msd = MIN(msd, isd)
      isd = sd(d, X, swp(d, m, n, [3, 2, 1], [i, j, k], rot(1), Y)) ; msd = MIN(msd, isd)
    enddo
    enddo
    enddo
!
    Z = Y
    call bra%swap(Z)
    call u%assert_almost_equal(msd,         bra%upperbound(W), 'branchcut vs brute')
    call u%assert_almost_equal(sd(d, X, Z), bra%upperbound(W), 'swap')
!
  end subroutine test1
!
  subroutine test2()
    integer, parameter          :: d = 3
    integer, parameter          :: s = 3
    integer, parameter          :: m1 = 5, n1 = 3, f1 = 3, g1 = 2
    integer, parameter          :: m2 = 3, n2 = 4, f2 = 2, g2 = 2
    integer, parameter          :: m3 = 7, n3 = 5, f3 = 7, g3 = 3
    integer, parameter          :: mn = m1 * n1 + m2 * n2 + m3 * n3
    type(mol_block)             :: b(3) = [mol_block(0, 3, m1, n1, f1, g1), &
                                        &  mol_block(0, 1, m2, n2, f2, g2), &
                                        &  mol_block(0, 2, m3, n3, f3, g3)]
    type(branch_and_prune)      :: bra
    type(mol_block_list)        :: blk
    type(molecular_rotation)    :: rot(s)
    real(RK)                    :: X(d, mn), Y(d, mn), Z(d, mn)
    real(RK), allocatable       :: W(:)
    integer                     :: i
!
    rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m1, 2]))
    rot(2) = molecular_rotation(RESHAPE([(i, i=1,0)], [0, 1]))
    rot(3) = molecular_rotation(RESHAPE([7, 6, 5, 4, 3, 2, 1], [m3, 1]))
    blk = mol_block_list(d, s, b)
!
    X = sample(d, mn)
    Y = sample(d, mn)
!
    bra = branch_and_prune(blk, 1, rot)
    allocate (W(bra%memsize()))
    call bra%setup(X, Y, W)
    call bra%run(W)
!
    Z=Y
    call bra%swap(Z)
    call u%assert_almost_equal(sd(d, X, Z), bra%upperbound(W), 'multiple swap')
!
  end subroutine test2
!
  pure function swp(d, m, n, per, sym, rot, X) result(res)
    integer(IK), intent(in) :: d, m, n, per(:), sym(:)
    type(molecular_rotation), intent(in) :: rot
    real(RK), intent(in)    :: X(d, m, n)
    real(RK)                :: tmp(d, m, n), res(d, m * n)
    integer(IK)             :: i
    tmp = X
    do i = 1, SIZE(per)
      tmp(:, :, per(i)) = X(:, :, i)
      call rot%swap(d, tmp(:, :, per(i)), sym(i))
    end do
    res = RESHAPE(tmp, [d, m * n])
  end function swp
!
  pure function sd(d, X, Y) result(res)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: X(:, :), Y(:, :)
    real(RK)                :: C(d, d), R(d, d), W(100), res
    C = MATMUL(Y, TRANSPOSE(X))
    call estimate_rotation_matrix(d, SUM(X * X) + SUM(Y * Y), C, R, W)
    res = SUM(X**2) + SUM(Y**2) - 2 * SUM(C * R)
  end function sd
!
  function sample(d, n) result(res)
    integer, intent(in)  :: d, n
    real(RK)             :: cnt(d)
    real(RK)             :: res(d, n)
    integer              :: i
    call RANDOM_NUMBER(res)
    cnt = SUM(res, 2) / n
    do concurrent(i=1:n)
      res(:, i) = res(:, i) - cnt
    enddo
  end function sample
!
end program main

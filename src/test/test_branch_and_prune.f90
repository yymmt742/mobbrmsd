program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_Kabsch
  use mod_molecular_rotation
  use mod_branch_and_prune
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test node')
  call test1()
! call test2()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter          :: d = 3
    integer, parameter          :: s = 1
    integer, parameter          :: m = 5, n = 8, f = 5, g = 2
    integer, parameter          :: mn = m * n
    type(mol_block)             :: b(1) = [mol_block(0, 2, m, n, f, g)]
    type(branch_and_prune)      :: bra
    type(mol_block_list)        :: blk
    type(molecular_rotation)    :: rot(s)
    real(RK)                    :: X(d, mn), Y(d, mn), Z(d, mn)
    real(RK), allocatable       :: W(:)
!
!   rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m, 0]))
    rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5], [m, 1]))
!   rot(1) = molecular_rotation(RESHAPE([(i, i=1, 0)], [m, 0]))
    blk = mol_block_list(d, s, b)
!
    X = sample(d,mn)
    Y = 0.0D0 * X + 1.0D0 * sample(d, mn)
    !Y = sample(d, mn)
!
    bra = branch_and_prune(blk, 1, rot)
    allocate (W(bra%memsize()))
    call bra%setup(X, Y, W)
    call bra%run(W)
!
    Z=Y
    call bra%swap(Z)
!print'(3f9.3)', Z(:,:15)-Y(:,:15)
!print*
    print'(F9.3)', bra%upperbound(W)
    print'(f9.3)', sd(d, X, Z)
!
    Z = swp(d, m, n, [1, 2], [0, 0], rot(1), Y)
    print'(4I4, f9.3)', 1, 2, 0, 0, sd(d, X, Z)
!print'(3f9.3)', Z(:,:15)-Y(:,:15)
!print*
    Z = swp(d, m, n, [1, 2], [1, 0], rot(1), Y)
    print'(4I4, f9.3)', 1, 2, 1, 0, sd(d, X, Z)
!print'(3f9.3)', Z(:,:15)-Y(:,:15)
!print*
    Z = swp(d, m, n, [1, 2], [0, 1], rot(1), Y)
    print'(4I4, f9.3)', 1, 2, 0, 1, sd(d, X, Z)
!print'(3f9.3)', Z(:,:15)-Y(:,:15)
!print*
    Z = swp(d, m, n, [1, 2], [1, 1], rot(1), Y)
    print'(4I4, f9.3)', 1, 2, 1, 1, sd(d, X, Z)
!print'(3f9.3)', Z(:,:15)-Y(:,:15)
!print*
    Z = swp(d, m, n, [2, 1], [0, 0], rot(1), Y)
    print'(4I4, f9.3)', 2, 1, 0, 0, sd(d, X, Z)
!print'(3f9.3)', Z(:,:15)-Y(:,:15)
!print*
    Z = swp(d, m, n, [2, 1], [1, 0], rot(1), Y)
    print'(4I4, f9.3)', 2, 1, 1, 0, sd(d, X, Z)
!print'(3f9.3)', Z(:,:15)-Y(:,:15)
!print*
    Z = swp(d, m, n, [2, 1], [0, 1], rot(1), Y)
    print'(4I4, f9.3)', 2, 1, 0, 1, sd(d, X, Z)
!print'(3f9.3)', Z(:,:15)-Y(:,:15)
!print*
    Z = swp(d, m, n, [2, 1], [1, 1], rot(1), Y)
    print'(4I4, f9.3)', 2, 1, 1, 1, sd(d, X, Z)
!print'(3f9.3)', Z(:,:15)-Y(:,:15)
!print*
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
    real(RK)                    :: C(d, d), R(d, d)
    real(RK), allocatable       :: W(:)
    integer                     :: i
!
    rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m1, 2]))
    rot(2) = molecular_rotation(RESHAPE([(i, i=1,0)], [0, 1]))
    rot(3) = molecular_rotation(RESHAPE([7, 6, 5, 4, 3, 2, 1], [m3, 1]))
    blk = mol_block_list(d, s, b)
!
    X = sample(d,mn)
    Y = 0.3D0 * X + 0.7D0 * sample(d, mn)
    !Y = sample(d, mn)
!
    bra = branch_and_prune(blk, 1, rot)
    allocate (W(bra%memsize()))
    call bra%setup(X, Y, W)
    call bra%run(W)
!
    Z=Y
    call bra%swap(Z)
print'(3f9.3)', Z(:,:15)-Y(:,:15)
print*
print'(3f9.3)', Z(:,16:27)-Y(:,16:27)
print*
print'(3f9.3)', Z(:,28:62)-Y(:,28:62)
print*
!
    print'(F9.3)', bra%upperbound(W)
!
    C = MATMUL(Z, TRANSPOSE(X))
    call Kabsch(d, C, R, W)
    print'(3f9.3)', R
    print'(F9.3)', SUM(X**2) + SUM(Y**2) - 2 * SUM(C * R)
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
      tmp(:, :, i) = X(:, :, per(i))
      call rot%swap(d, tmp(:, :, i), sym(per(i)))
    end do
    res = RESHAPE(tmp, [d, m * n])
  end function swp
!
  pure function sd(d, X, Y) result(res)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: X(:, :), Y(:, :)
    real(RK)                :: C(d, d), R(d, d), W(100), res
    C = MATMUL(Y, TRANSPOSE(X))
    call Kabsch(d, C, R, W)
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

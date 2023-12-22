program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_molecular_rotation
  use mod_branch_and_prune
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test node')
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
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
    real(RK)                    :: X(d, mn), Y(d, mn)
    real(RK), allocatable       :: W(:)
    integer                     :: i
!
    rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m1, 2]))
    rot(2) = molecular_rotation(RESHAPE([(i, i=1,0)], [0, 1]))
    rot(3) = molecular_rotation(RESHAPE([7, 6, 5, 4, 3, 2, 1], [m3, 1]))
    blk = mol_block_list(d, s, b)
!
    X = sample(d,mn)
    !Y = sample(d, mn)
    Y = 1.0D0 * X + 0.0D0 * sample(d, mn)
!
    bra = branch_and_prune(blk, rot, 1)
    allocate (W(bra%memsize()))
    call bra%setup(X, Y, W)
    call bra%run(W)
!
  end subroutine test1
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

program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_molecular_rotation
  use mod_d_matrix
  use mod_tree
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test node')
  call test1()
! call test2()
! call test3()
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
    type(mol_block)             :: b(3) = [mol_block(0, m1, n1, f1, g1), &
                                        &  mol_block(0, m2, n2, f2, g2), &
                                        &  mol_block(0, m3, n3, f3, g3)]
    type(mol_block_list)        :: blk
    type(molecular_rotation)    :: rot(s)
    real(RK)                    :: X(d, mn), Y(d, mn)
    real(RK), allocatable       :: W(:)
    type(d_matrix_list)         :: dmat
    type(node)                  :: a
    type(breadth)               :: z1, z2, z3, z4, z5, z6, z7, z8
    integer                     :: i, ml
!
    rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m1, 2]))
    rot(2) = molecular_rotation(RESHAPE([(i, i=1,0)], [0, 1]))
    rot(3) = molecular_rotation(RESHAPE([7, 6, 5, 4, 3, 2, 1], [m3, 1]))
    blk = mol_block_list(d, s, b, rot)
    dmat = d_matrix_list(s, blk%b, rot, d, 1)
    allocate(W(dmat%memsize()))
!
    X = sample(d,mn)
    Y = 0.9D0 * X + 0.1D0 * sample(d, mn)
!
    call dmat%eval(X, Y, W)
!
!   a = node(blk, [x], [y])
!   print'(*(f16.3))',a%lowerbound
!   z1 = a%generate_breadth(rot, [x], [y])
!   print'(*(f16.3))',z1%nodes%lowerbound
!   ml = MINLOC(z1%nodes%lowerbound, 1)
!   z2 = z1%nodes(ml)%generate_breadth(rot, [x], [y])
!   print'(*(f16.3))',z2%nodes%lowerbound
!   ml = MINLOC(z2%nodes%lowerbound, 1)
!   z3 = z2%nodes(ml)%generate_breadth(rot, [x], [y])
!   print'(*(f16.3))',z3%nodes%lowerbound
!   ml = MINLOC(z3%nodes%lowerbound, 1)
!   z4 = z3%nodes(ml)%generate_breadth(rot, [x], [y])
!   print'(*(f16.3))',z4%nodes%lowerbound
!
!   ml = MINLOC(z4%nodes%lowerbound, 1)
!   z5 = z4%nodes(ml)%generate_breadth(rot, [x], [y])
!   print'(*(f16.3))',z5%nodes%lowerbound
!
!   ml = MINLOC(z5%nodes%lowerbound, 1)
!   z6 = z5%nodes(ml)%generate_breadth(rot, [x], [y])
!   print'(*(f16.3))',z6%nodes%lowerbound
!
!   ml = MINLOC(z6%nodes%lowerbound, 1)
!   z7 = z6%nodes(ml)%generate_breadth(rot, [x], [y])
!   print'(*(f16.3))',z7%nodes%lowerbound
!
!   ml = MINLOC(z7%nodes%lowerbound, 1)
!   z8 = z7%nodes(ml)%generate_breadth(rot, [x], [y])
!   print'(*(f16.3))',z8%nodes%lowerbound
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

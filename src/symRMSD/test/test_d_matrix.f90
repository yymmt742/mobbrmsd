program main
  use mod_params, only: D, DD, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_mol_symmetry
  use mod_estimate_rotation_matrix
  use mod_d_matrix
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test d_matrix')
!
  call setup_dimension(3)
!
  call test0()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    type(mol_block)       :: b(3)
    type(mol_symmetry)    :: ms(3)
    type(c_matrix)        :: a(3)
    real(RK)              :: X(3, 5 * 3 + 3 * 4 + 8 * 3)
    real(RK)              :: Y(3, 5 * 3 + 3 * 2 + 8 * 5)
    real(RK), allocatable :: W(:)
!
    X = sample(3, SIZE(X, 2))
    Y = sample(3, SIZE(Y, 2))
!
    b(1) = mol_block(1, 5, 3, 3)
    b(2) = mol_block(3, 3, 4, 2)
    b(3) = mol_block(2, 8, 3, 5)
    call mol_block_list_init(b)
!
    print*,b(1)
    ms(2) = mol_symmetry(RESHAPE([2, 3, 1, 3, 1, 2], [3, 2]))
    ms(3) = mol_symmetry(RESHAPE([1, 2, 3, 4, 5, 6, 7, 8], [8, 1]))
!
    a = c_matrix(b)
    print*,memsize_c_matrix(a)
    print*,worksize_c_matrix(a)
    a(1)%pc = 1
    a(2)%pc = memsize_c_matrix(a(1)) + 1
    a(3)%pc = memsize_c_matrix(a(1)) + memsize_c_matrix(a(2)) + 1
    allocate (W(SUM(memsize_c_matrix(a)) + SUM(worksize_c_matrix(a))))
    W(:) = 999
    call c_matrix_eval(a(1), b(1), ms(1), X, Y, W)
    call c_matrix_eval(a(2), b(2), ms(2), X, Y, W)
    call c_matrix_eval(a(3), b(3), ms(3), X, Y, W)
    print'(10f9.3)',W
!
  end subroutine test0
end program main

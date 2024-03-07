program main
  use mod_params, only: setup_dimension, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_rotation_matrix
  use mod_c_matrix
  use mod_testutil
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
    type(mol_block_tuple) :: b(3)
    type(c_matrix)        :: c(3)
    real(RK)              :: X(3, 5 * 3 + 3 * 4 + 8 * 3)
    real(RK)              :: Y(3, 5 * 3 + 3 * 2 + 8 * 5)
    real(RK), allocatable :: W(:)
    integer(IK)           :: x1, x2, x3
    integer(IK)           :: p1, p2, p3, w1, w2, w3
!
    X = sample(3, SIZE(X, 2))
    Y = sample(3, SIZE(Y, 2))
!
    b(1) = mol_block_tuple(5, 3)
    b(2) = mol_block_tuple(3, 2, sym=RESHAPE([2, 3, 1, 3, 1, 2], [3, 2]))
    b(3) = mol_block_tuple(8, 3, sym=RESHAPE([1, 2, 3, 4, 5, 6, 7, 8], [8, 1]))
!
    x1 = 1
    x2 = x1 + mol_block_natm(b(1)%b)
    x3 = x2 + mol_block_natm(b(2)%b)
!
    print*,mol_block_nsym(b%b)
    c = c_matrix(b%b)
    print*,c_matrix_memsize(c)
    print*,c_matrix_worksize(c)
    p1 = 1
    w1 = p1 + c_matrix_memsize(c(1))
    p2 = w1
    w2 = p2 + c_matrix_memsize(c(2))
    p3 = w2
    w3 = p3 + c_matrix_memsize(c(3))
    allocate (W(SUM(c_matrix_memsize(c)) + c_matrix_worksize(c(3))))
    W(:) = 999
    call c_matrix_eval(c(1), b(1)%b, b(1)%w, X(1, x1), Y(1, x1), W(p1), W(w1))
    call c_matrix_eval(c(2), b(2)%b, b(2)%w, X(1, x2), Y(1, x2), W(p2), W(w2))
    call c_matrix_eval(c(3), b(3)%b, b(3)%w, X(1, x3), Y(1, x3), W(p3), W(w3))
    print'(10f5.1)',W(p1:p1+89)
    print*
    print'(14f5.1)',W(p2:p2+111)
    print*
    print'(19f5.1)',W(p3:p3+170)
    print*
!
  end subroutine test0
!
end program main

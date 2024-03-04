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
! call test1()
! call test2()
! call test3()
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
    integer(IK)           :: p
!
    X = sample(3, SIZE(X, 2))
    Y = sample(3, SIZE(Y, 2))
!
    b(1) = mol_block_tuple(5, 3)
    b(2) = mol_block_tuple(3, 2, sym=RESHAPE([2, 3, 1, 3, 1, 2], [3, 2]))
    b(3) = mol_block_tuple(8, 3, sym=RESHAPE([1, 2, 3, 4, 5, 6, 7, 8], [8, 1]))
!
    p = 1 + mol_block_total_size(b(1)%b)
    call mol_block_set_pointer(b(2)%b, p, 1)
    p = 1 + mol_block_total_size(b(2)%b)
    call mol_block_set_pointer(b(3)%b, p, 1)
!
print*,mol_block_nsym(b%b)
    c = c_matrix(b%b)
    print*,c_matrix_memsize(c)
    print*,c_matrix_worksize(c)
    c(1)%p = 1
    c(1)%w = c(1)%p + c_matrix_memsize(c(1))
    c(2)%p = c(1)%w
    c(2)%w = c(2)%p + c_matrix_memsize(c(2))
    c(3)%p = c(2)%w
    c(3)%w = c(3)%p + c_matrix_memsize(c(3))
    allocate (W(SUM(c_matrix_memsize(c)) + c_matrix_worksize(c(3))))
    W(:) = 999
    call c_matrix_eval(c(1), b(1)%b, b(1)%w, X, Y, W, W)
    call c_matrix_eval(c(2), b(2)%b, b(2)%w, X, Y, W, W)
    call c_matrix_eval(c(3), b(3)%b, b(3)%w, X, Y, W, W)
    print'(10f5.1)',W(c(1)%p:c(1)%p+89)
    print*
    print'(14f5.1)',W(c(2)%p:c(2)%p+111)
    print*
    print'(19f5.1)',W(c(3)%p:c(3)%p+170)
    print*
!
  end subroutine test0
!
end program main

program main
  use mod_params, only: setup_dimension, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_c_matrix
  use mod_f_matrix
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test f_matrix')
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
    type(f_matrix)        :: f(3)
    real(RK)              :: X(3, 5 * 3 + 3 * 4 + 8 * 3)
    real(RK)              :: Y(3, 5 * 3 + 3 * 2 + 8 * 5)
    real(RK), allocatable :: W(:)
!
    X = sample(3, SIZE(X, 2))
    Y = sample(3, SIZE(Y, 2))
!
    b(1) = mol_block_tuple(5, 3)
    b(2) = mol_block_tuple(3, 2, sym=RESHAPE([2, 3, 1, 3, 1, 2], [3, 2]))
    b(3) = mol_block_tuple(8, 3, sym=RESHAPE([1, 2, 3, 4, 5, 6, 7, 8], [8, 1]))
!
    c = c_matrix(b%b)
    f = f_matrix(b%b)
!
    print*,c_matrix_memsize(c)
    print*,c_matrix_worksize(c)
    print*,f_matrix_memsize(f)
    print*,f_matrix_worksize(f)
!
    c(1)%p = 1
    c(1)%w = c(1)%p + c_matrix_memsize(c(1))
    f(1)%p = c(1)%w
    f(1)%w = f(1)%p + f_matrix_memsize(f(1))
    c(2)%p = f(1)%w
    c(2)%w = c(2)%p + c_matrix_memsize(c(2))
    f(2)%p = c(2)%w
    f(2)%w = f(2)%p + f_matrix_memsize(f(2))
    c(3)%p = f(2)%w
    c(3)%w = c(3)%p + c_matrix_memsize(c(3))
    f(3)%p = c(3)%w
    f(3)%w = f(3)%p + f_matrix_memsize(f(3))
!
    allocate (W(SUM(c_matrix_memsize(c) + c_matrix_worksize(c) &
   &              + f_matrix_memsize(f) + f_matrix_worksize(f))))
!
    W(:) = 999
    call c_matrix_eval(c(1), b(1)%b, b(1)%w, X, Y, W, W)
    call f_matrix_eval(f(1), b(1)%b, c(1), W, W)
    call c_matrix_eval(c(2), b(2)%b, b(2)%w, X, Y, W, W)
    call f_matrix_eval(f(2), b(2)%b, c(2), W, W)
    call c_matrix_eval(c(3), b(3)%b, b(3)%w, X, Y, W, W)
    call f_matrix_eval(f(3), b(3)%b, c(3), W, W)
!
    print'(3f9.3)', W(f(1)%p:f(1)%p + f_matrix_memsize(f(1)) - 1)
    print'(4f9.3)', W(f(2)%p:f(2)%p + f_matrix_memsize(f(2)) - 1)
    print'(5f9.3)', W(f(3)%p:f(3)%p + f_matrix_memsize(f(3)) - 1)
!
  end subroutine test0
!
end program main


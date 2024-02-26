program main
  use mod_params, only: setup_dimension, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_mol_symmetry
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
    type(mol_block)       :: b(3)
    type(mol_symmetry)    :: ms(3)
    type(c_matrix)        :: c(3)
    type(f_matrix)        :: s(3)
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
    ms(2) = mol_symmetry(RESHAPE([2, 3, 1, 3, 1, 2], [3, 2]))
    ms(3) = mol_symmetry(RESHAPE([1, 2, 3, 4, 5, 6, 7, 8], [8, 1]))
!
    c = c_matrix(b)
    s = f_matrix(b)
    print*,memsize_c_matrix(c)
    print*,worksize_c_matrix(c)
    print*,memsize_f_matrix(s)
    print*,worksize_f_matrix(s)
!
    c(1)%p = 1
    c(1)%w = c(1)%p + memsize_c_matrix(c(1))
    s(1)%p = c(1)%w
    s(1)%w = s(1)%p + memsize_f_matrix(s(1))
    c(2)%p = s(1)%w
    c(2)%w = c(2)%p + memsize_c_matrix(c(2))
    s(2)%p = c(2)%w
    s(2)%w = s(2)%p + memsize_f_matrix(s(2))
    c(3)%p = s(2)%w
    c(3)%w = c(3)%p + memsize_c_matrix(c(3))
    s(3)%p = c(3)%w
    s(3)%w = s(3)%p + memsize_f_matrix(s(3))
!
    allocate (W(SUM(memsize_c_matrix(c) + worksize_c_matrix(c) &
   &              + memsize_f_matrix(s) + worksize_f_matrix(s))))
!
    W(:) = 999
    call c_matrix_eval(c(1), b(1), ms(1), X, Y, W, W)
    call f_matrix_eval(s(1), b(1), W(c(1)%p), W, W)
    call c_matrix_eval(c(2), b(2), ms(2), X, Y, W, W)
    call f_matrix_eval(s(2), b(2), W(c(2)%p), W, W)
    call c_matrix_eval(c(3), b(3), ms(3), X, Y, W, W)
    call f_matrix_eval(s(3), b(3), W(c(3)%p), W, W)
!
    print'(3f9.3)', W(s(1)%p:s(1)%p + memsize_f_matrix(s(1)) - 1)
    print'(4f9.3)', W(s(2)%p:s(2)%p + memsize_f_matrix(s(2)) - 1)
    print'(5f9.3)', W(s(3)%p:s(3)%p + memsize_f_matrix(s(3)) - 1)
!
  end subroutine test0
!
end program main


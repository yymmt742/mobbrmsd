program main
  use mod_params, only: D, RK, IK, ONE => RONE, ZERO => RZERO
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
  call test0(9, 1, 1, [1])
  call test0(4, 1, 2, [1,3,4,2])
  call test0(4, 3, 2, [1,3,4,2])
  call test0(3, 3, 3, [1,3,2,2,3,1])
  call test0(1, 4, 1, [1])
  call test0(3, 4, 1, [1])
!
  call test1(0.0_RK)
  call test1(1.0_RK)
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0(m, n, s, sym)
    integer(IK), intent(in) :: m, n, s, sym(m * (s - 1))
    type(mol_block)         :: b
    type(c_matrix)          :: c
    type(f_matrix)          :: f
    real(RK)                :: X(D, m, n), Y(D, m, n)
    real(RK), allocatable   :: CX(:), CW(:), FX(:), FW(:)
!
    X = sample(m, n)
    Y = X
    b = mol_block(m, n)
    c = c_matrix(b%q)
    f = f_matrix(b%q)
!
    allocate (CX(c_matrix_memsize(c%q)))
    allocate (CW(c_matrix_worksize(c%q)))
    allocate (FX(f_matrix_memsize(f%q)))
    allocate (FW(f_matrix_worksize(f%q)))
    CX(:) = 999
    CW(:) = 999
    FX(:) = 999
    FW(:) = 999
!
    call c_matrix_eval(c%q, b%q, X, Y, CX, CW)
    call f_matrix_eval(f%q, c%q, CX, FX, FW)
!
    print'(4f9.3)', FX
    print*
!   print'(4f9.3)', FW
!
  end subroutine test0
!
  subroutine test1(s)
    real(RK), intent(in)  :: s
    type(mol_block)       :: b(3)
    type(c_matrix)        :: c(3)
    type(f_matrix)        :: f(3)
    real(RK)              :: X(D, 5 * 3 + 3 * 2 + 8 * 3)
    real(RK)              :: Y(D, 5 * 3 + 3 * 2 + 8 * 3)
    real(RK), allocatable :: W(:)
    integer(IK)           :: nx
    integer(IK)           :: x1, x2, x3
    integer(IK)           :: c1, c2, c3
    integer(IK)           :: f1, f2, f3
    integer(IK)           :: v1, v2, v3
    integer(IK)           :: w1, w2, w3
!
    X = sample(SIZE(X, 2))
    Y = (ONE - s) * X + s * sample(SIZE(Y, 2))
!
    b(1) = mol_block(5, 3)
    b(2) = mol_block(3, 2, sym=RESHAPE([2, 3, 1, 3, 1, 2], [3, 2]))
    b(3) = mol_block(8, 3, sym=RESHAPE([1, 2, 3, 4, 5, 6, 7, 8], [8, 1]))
!
    c(1) = c_matrix(b(1)%q)
    c(2) = c_matrix(b(2)%q)
    c(3) = c_matrix(b(3)%q)
    f(1) = f_matrix(b(1)%q)
    f(2) = f_matrix(b(2)%q)
    f(3) = f_matrix(b(3)%q)
!
    print *, c_matrix_memsize(c(1)%q), c_matrix_memsize(c(2)%q), c_matrix_memsize(c(3)%q)
    print *, c_matrix_worksize(c(1)%q), c_matrix_worksize(c(2)%q), c_matrix_worksize(c(3)%q)
    print *, f_matrix_memsize(f(1)%q), f_matrix_memsize(f(2)%q), f_matrix_memsize(f(3)%q)
    print *, f_matrix_worksize(f(1)%q), f_matrix_worksize(f(2)%q), f_matrix_worksize(f(3)%q)
!
    x1 = 1
    x2 = x1 + mol_block_natm(b(1)%q)
    x3 = x2 + mol_block_natm(b(2)%q)
!
    c1 = 1
    v1 = c1 + c_matrix_memsize(c(1)%q)
    f1 = v1
    w1 = f1 + f_matrix_memsize(f(1)%q)
    c2 = w1
    v2 = c2 + c_matrix_memsize(c(2)%q)
    f2 = v2
    w2 = f2 + f_matrix_memsize(f(2)%q)
    c3 = w2
    v3 = c3 + c_matrix_memsize(c(3)%q)
    f3 = v3
    w3 = f3 + f_matrix_memsize(f(3)%q)
!
    nx = c_matrix_memsize(c(1)%q) &
   &   + f_matrix_memsize(f(1)%q) &
   &   + c_matrix_memsize(c(2)%q) &
   &   + f_matrix_memsize(f(2)%q) &
   &   + c_matrix_memsize(c(3)%q) &
   &   + f_matrix_memsize(f(3)%q) &
   &   + MAX(c_matrix_worksize(c(3)%q), f_matrix_worksize(f(3)%q))
    allocate (w(nx))
    W(:) = 999
!
    call c_matrix_eval(c(1)%q, b(1)%q, X(1, x1), Y(1, x1), W(c1), W(v1))
    call f_matrix_eval(f(1)%q, c(1)%q, W(c1), W(f1), W(w1))
    call c_matrix_eval(c(2)%q, b(2)%q, X(1, x2), Y(1, x2), W(c2), W(v2))
    call f_matrix_eval(f(2)%q, c(2)%q, W(c2), W(f2), W(w2))
    call c_matrix_eval(c(3)%q, b(3)%q, X(1, x3), Y(1, x3), W(c3), W(v3))
    call f_matrix_eval(f(3)%q, c(3)%q, W(c3), W(f3), W(w3))
!
    print'(3f9.3)', W(f1:f1 + f_matrix_memsize(f(1)%q) - 1)
    print'(2f9.3)', W(f2:f2 + f_matrix_memsize(f(2)%q) - 1)
    print'(3f9.3)', W(f3:f3 + f_matrix_memsize(f(3)%q) - 1)
    print*
!
  end subroutine test1
!
end program main


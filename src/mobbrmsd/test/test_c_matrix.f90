program main
  use blas_lapack_interface, only : D, DD
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_rotation
  use mod_c_matrix
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test c_matrix')
!
  call test0(0, 0, 1, [0])
  call test0(5, 1, 1, [0])
  call test0(5, 1, 2, [4, 5, 1, 2, 3])
  call test0(5, 3, 1, [0])
  call test0(2, 1, 2, [2, 1])
  call test0(3, 5, 3, [1, 3, 2, 3, 2, 1])
!
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0(m, n, s, sym)
    integer(IK), intent(in) :: m, n, s, sym(m * (s - 1))
    type(mol_block)         :: b
    type(c_matrix)          :: c
    real(RK)                :: X(D, m, n)
    real(RK)                :: Y(D, m, n)
    real(RK), allocatable   :: Z(:), W(:)
    integer(IK)             :: i, j, p
!
    X = sample(m, n)
    Y = sample(m, n)
!
    b = mol_block(m, n, sym=RESHAPE(sym, [m, (s - 1)]))
    c = c_matrix(b%q)
    allocate (Z(c_matrix_memsize(c%q)))
    allocate (W(c_matrix_worksize(c%q)))
    Z(:) = 999
    W(:) = 999
!
    call c_matrix_eval(c%q, b%q, X, Y, Z, W)
!
    p = 1
    do j = 1, n
      do i = 1, n
        call check_gcov(m, s, sym, Z(p), X(1, 1, i), Y(1, 1, j))
        p = p + 1 + DD * s
      end do
    end do
!
  end subroutine test0
!
  subroutine test1()
    type(mol_block)       :: b(3)
    type(c_matrix)        :: c(3)
    real(RK)              :: X(D, 5 * 3 + 3 * 4 + 8 * 3)
    real(RK)              :: Y(D, 5 * 3 + 3 * 2 + 8 * 5)
    real(RK), allocatable :: W(:)
    integer(IK)           :: nw, x1, x2, x3
    integer(IK)           :: p1, p2, p3, w1, w2, w3
!
    X = sample(SIZE(X, 2))
    Y = sample(SIZE(Y, 2))
!
    b(1) = mol_block(5, 3)
    b(2) = mol_block(3, 2, sym=RESHAPE([2, 3, 1, 3, 1, 2], [3, 2]))
    b(3) = mol_block(8, 3, sym=RESHAPE([1, 2, 6, 3, 7, 4, 5, 8], [8, 1]))
!
    x1 = 1
    x2 = x1 + mol_block_natm(b(1)%q)
    x3 = x2 + mol_block_natm(b(2)%q)
!
    print *, mol_block_nsym(b(1)%q), mol_block_nsym(b(2)%q), mol_block_nsym(b(3)%q)
    c(1) = c_matrix(b(1)%q)
    c(2) = c_matrix(b(2)%q)
    c(3) = c_matrix(b(3)%q)
!
    print *, c_matrix_memsize(c(1)%q), c_matrix_memsize(c(2)%q), c_matrix_memsize(c(3)%q)
    print *, c_matrix_worksize(c(1)%q), c_matrix_worksize(c(2)%q), c_matrix_worksize(c(3)%q)
    p1 = 1
    w1 = p1 + c_matrix_memsize(c(1)%q)
    p2 = w1
    w2 = p2 + c_matrix_memsize(c(2)%q)
    p3 = w2
    w3 = p3 + c_matrix_memsize(c(3)%q)
!
    nw = c_matrix_memsize(c(1)%q) &
   &   + c_matrix_memsize(c(2)%q) &
   &   + c_matrix_memsize(c(3)%q) &
   &   + c_matrix_worksize(c(3)%q)
    allocate (W(nw))
    W(:) = 999
!
    call c_matrix_eval(c(1)%q, b(1)%q, X(1, x1), Y(1, x1), W(p1), W(w1))
    call c_matrix_eval(c(2)%q, b(2)%q, X(1, x2), Y(1, x2), W(p2), W(w2))
    call c_matrix_eval(c(3)%q, b(3)%q, X(1, x3), Y(1, x3), W(p3), W(w3))
!
    print'(10f5.1)', W(p1:p1 + c_matrix_memsize(c(1)%q) - 1)
    print*
    print'(10f5.1)', W(p2:p2 + c_matrix_memsize(c(2)%q) - 1)
    print*
    print'(10f5.1)', W(p3:p3 + c_matrix_memsize(c(3)%q) - 1)
    print*
!
  end subroutine test1
!
  subroutine check_gcov(m, s, sym, C, X, Y)
    integer(IK), intent(in) :: m, s, sym(m, s - 1)
    real(RK), intent(in)    :: C(*), X(D, m), Y(D, m)
    integer(IK)             :: i
!
    call u%assert_almost_equal(C(1), SUM(X * X) + SUM(Y * Y),          'auto variance')
    call u%assert_almost_equal(C(2:1 + DD), [MATMUL(Y, TRANSPOSE(X))], 'covariance 1 ')
!
    do i = 1, s - 1
      call u%assert_almost_equal(C(2 + DD * i:1 + DD * (i + 1)), [MATMUL(Y(:, sym(:, i)), TRANSPOSE(X))], 'covariance i ')
    end do
!
  end subroutine check_gcov
!
end program main

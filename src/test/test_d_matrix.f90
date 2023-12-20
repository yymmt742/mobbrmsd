program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_molecular_rotation
  use mod_d_matrix
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test d_matrix')
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter         :: d = 3
    integer, parameter         :: s = 3
    integer, parameter         :: m = 5
    integer, parameter         :: n = 7
    integer, parameter         :: f = 3
    integer, parameter         :: g = 5
    integer, parameter         :: mn = m * n
    integer, parameter         :: swp(m, s - 1) = RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m, s - 1])
    type(mol_block), parameter :: b = mol_block(1, m, n, f, g)
    type(molecular_rotation)   :: rot
    type(d_matrix)             :: a
    real(RK)                   :: X(d, mn), Y(d, mn)
    real(RK)                   :: LF, LB, H, C(d, d), R(d, d)
    real(RK), allocatable      :: w(:)
!
    rot = molecular_rotation(swp)
    a = d_matrix(1, 1, d, s, b)
    X = sample(d, mn)
    Y = 0.9D0 * MATMUL(SO3(), X) + 0.1D0 * sample(d, mn)
    print *, d_matrix_memsize(a)
    allocate (w(d_matrix_memsize(a)))
    call d_matrix_eval(a, rot, X, Y, W)
    call d_matrix_partial_eval(a, 1, 1, 1, [2,3,4,5], W, LF, LB, H, C, R)
print'(3f9.3)',LF, LB, LF+LB
    call d_matrix_partial_eval(a, 2, 2, 1, [3,4,5], W, LF, LB, H, C, R)
print'(3f9.3)',LF, LB, LF+LB
    call d_matrix_partial_eval(a, 3, 3, 1, [4,5], W, LF, LB, H, C, R)
print'(3f9.3)',LF, LB, LF+LB
    call d_matrix_partial_eval(a, 4, 4, 1, [5], W, LF, LB, H, C, R)
print'(3f9.3)',LF, LB, LF+LB
    call d_matrix_partial_eval(a, 5, 5, 1, [5], W, LF, LB, H, C, R)
print'(3f9.3)',LF, LB, LF+LB
print*
    call d_matrix_partial_eval(a, 1, 1, 2, [2,3,4,5], W, LF, LB, H, C, R)
print'(3f9.3)',LF, LB, LF+LB
    call d_matrix_partial_eval(a, 1, 2, 1, [1,3,4,5], W, LF, LB, H, C, R)
print'(3f9.3)',LF, LB, LF+LB
    call d_matrix_partial_eval(a, 1, 2, 2, [1,3,4,5], W, LF, LB, H, C, R)
print'(3f9.3)',LF, LB, LF+LB
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
  function SO3() result(res)
    real(RK) :: a(3), res(3, 3)
    call RANDOM_NUMBER(a)
    a = a / SQRT(DOT_PRODUCT(a, a))
    res(:, 1) = [a(1) * a(1), a(1) * a(2) - a(3), a(1) * a(3) + a(2)]
    res(:, 2) = [a(1) * a(2) + a(3), a(2) * a(2), a(2) * a(3) - a(1)]
    res(:, 3) = [a(1) * a(3) - a(2), a(2) * a(3) + a(1), a(3) * a(3)]
  end function SO3
!
end program main

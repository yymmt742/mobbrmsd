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
    integer, parameter         :: s = 2
    integer, parameter         :: m = 5
    integer, parameter         :: n = 7
    integer, parameter         :: f = 3
    integer, parameter         :: g = 5
    integer, parameter         :: mn = m * n
    integer, parameter         :: swp(m,s) = RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m, s])
    type(mol_block), parameter :: b = mol_block(1, m, n, f, g)
    type(molecular_rotation)   :: rot
    type(d_matrix)             :: a
    real(RK)                   :: X(d, mn), Y(d, mn)
    real(RK), allocatable      :: w(:)
!
    rot = molecular_rotation(swp)
    a = d_matrix(1, d, s, b)
    X = sample(d, mn)
    Y = 0.8D0 * X + 0.2D0 * sample(d, mn)
    print *, d_matrix_memsize(a)
    allocate (w(d_matrix_memsize(a)))
    call d_matrix_eval(a, rot, X, Y, W)
    print*, d_matrix_partial_lower_bound(a, 0, [1], W)
    print*, d_matrix_partial_lower_bound(a, 1, [5], W)
    print*, d_matrix_partial_lower_bound(a, 2, [1, 2], W)
    print*, d_matrix_partial_lower_bound(a, 2, [4, 5], W)
    print*, d_matrix_partial_lower_bound(a, 3, [1, 3, 5], W)
    print*, d_matrix_partial_lower_bound(a, 3, [3, 4, 5], W)
    print*, d_matrix_partial_lower_bound(a, 4, [1, 2, 3, 5], W)
    print*, d_matrix_partial_lower_bound(a, 4, [2, 3, 4, 5], W)
    print*, d_matrix_partial_lower_bound(a, 5, [1, 2, 3, 4, 5], W)
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

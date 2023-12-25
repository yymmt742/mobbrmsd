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
  call test2()
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
    type(mol_block), parameter :: b = mol_block(1, s, m, n, f, g)
    type(molecular_rotation)   :: rot
    type(d_matrix)             :: a
    real(RK)                   :: X(d, mn), Y(d, mn)
    real(RK)                   :: LT, LF, LB, H, C(d, d), R(d, d)
    real(RK), allocatable      :: w(:)
!
    rot = molecular_rotation(swp)
    C = 0D0
    H = 0D0
    a = d_matrix(1, d, b)
    X = sample(d, mn)
    Y = 0.9D0 * MATMUL(SO3(), X) + 0.1D0 * sample(d, mn)
    print *, d_matrix_memsize(a)
    allocate (w(d_matrix_memsize(a)))
    call d_matrix_eval(a, rot, X, Y, W)
    call d_matrix_partial_eval(a, 1, 1, 1, [2, 3, 4, 5], W, LT, H, C, LF, LB, R)
    print'(4f9.3)', LF, LB, LF + LB, LT
    call d_matrix_partial_eval(a, 2, 2, 1, [3, 4, 5], W, LT, H, C, LF, LB, R)
    print'(4f9.3)', LF, LB, LF + LB, LT
    call d_matrix_partial_eval(a, 3, 3, 1, [4, 5], W, LT, H, C, LF, LB, R)
    print'(4f9.3)', LF, LB, LF + LB, LT
    call d_matrix_partial_eval(a, 4, 4, 1, [5], W, LT, H, C, LF, LB, R)
    print'(4f9.3)', LF, LB, LF + LB, LT
    call d_matrix_partial_eval(a, 5, 5, 1, [5], W, LT, H, C, LF, LB, R)
    print'(4f9.3)', LF, LB, LF + LB, LT
    print *
    call d_matrix_partial_eval(a, 1, 1, 2, [2, 3, 4, 5], W, LT, H, C, LF, LB, R)
    print'(4f9.3)', LF, LB, LF + LB, LT
    call d_matrix_partial_eval(a, 1, 2, 1, [1, 3, 4, 5], W, LT, H, C, LF, LB, R)
    print'(4f9.3)', LF, LB, LF + LB, LT
    call d_matrix_partial_eval(a, 1, 2, 2, [1, 3, 4, 5], W, LT, H, C, LF, LB, R)
    print'(4f9.3)', LF, LB, LF + LB, LT
!
  end subroutine test1
!
  subroutine test2()
    integer, parameter       :: l = 3
    integer, parameter       :: d = 3
    integer, parameter       :: s = 3
    integer, parameter       :: m = 5
    integer, parameter       :: n = 7
    integer, parameter       :: f = 3
    integer, parameter       :: g = 5
    integer, parameter       :: mnl = (3 * m + 6) * n
    integer, parameter       :: swp(m, s - 1) = RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m, s - 1])
    integer                  :: perm(3 * g - 6)
    type(mol_block_list)     :: blk
    type(mol_block)          :: b(l)
    type(molecular_rotation) :: rot(l)
    type(d_matrix_list)      :: a
    real(RK)                 :: X(d, mnl), Y(d, mnl)
    real(RK), allocatable    :: w(:)
    real(RK)                 :: C(d * d), R(d * d), H, LT, LF, LB
    integer                  :: i, j, k
!
    do concurrent(i=1:l)
      rot(i) = molecular_rotation(swp(:, :i - 1))
    end do
!
    do i = 1, l
      b(i) = mol_block(0, rot(i)%n_sym() + 1, m + i, n, f, g - i)
    end do
!
    k = 0
    do j = 1, l
      do i = 1, b(j)%g
        k = k + 1
        perm(k) = i
      end do
    end do
    print*,perm
!
    blk = mol_block_list(d, l, b)
    a = d_matrix_list(blk, 1)
!
    X = sample(d, mnl)
    Y = 0.99D0 * MATMUL(SO3(), X) + 0.01D0 * sample(d, mnl)
    print *, a%memsize()
    print *, d_matrix_memsize(a%m)
!
    allocate (w(a%memsize()))
!
    call a%eval(rot, X, Y, W)
    print *, W(1)
    print'(3f9.3)', W(2:10)
    print'(*(f9.3))', W(a%o:a%o+l-1)
    print*
!
    H = W(a%h)
    C = W(a%c:a%c + a%dd - 1)
    call a%partial_eval(1, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(2, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(3, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(4, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(5, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(6, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(7, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(8, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(9, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    print*
!
    H = W(a%h)
    C = W(a%c:a%c + a%dd - 1)
    call a%partial_eval(1, perm, 4, 1, W, LT, H, C, LF, LB, R)
    perm(:4) = [4,1,2,3]
    print'(3f9.3)',H, LT
    call a%partial_eval(2, perm, 3, 1, W, LT, H, C, LF, LB, R)
    perm(:4) = [4,3,1,2]
    print'(3f9.3)',H, LT
    call a%partial_eval(3, perm, 2, 1, W, LT, H, C, LF, LB, R)
    perm(:4) = [4,3,2,1]
    print'(3f9.3)',H, LT
    call a%partial_eval(4, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(5, perm, 3, 1, W, LT, H, C, LF, LB, R)
    perm(5:7) = [3,1,2]
    print'(3f9.3)',H, LT
    call a%partial_eval(6, perm, 2, 1, W, LT, H, C, LF, LB, R)
    perm(5:7) = [3,2,1]
    print'(3f9.3)',H, LT
    call a%partial_eval(7, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
    call a%partial_eval(8, perm, 2, 1, W, LT, H, C, LF, LB, R)
    perm(8:9) = [2,1]
    print'(3f9.3)',H, LT
    call a%partial_eval(9, perm, 1, 1, W, LT, H, C, LF, LB, R)
    print'(3f9.3)',H, LT
!
  end subroutine test2
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

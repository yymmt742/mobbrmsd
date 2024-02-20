program main
  use mod_params, only: setup_dimension, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_mol_symmetry
  use mod_estimate_rotation_matrix
  use mod_c_matrix
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
    print*,c_matrix_memsize(a)
    print*,c_matrix_worksize(a)
    call c_matrix_init(a(1))
    call c_matrix_init(a(2), c_matrix_memsize(a(1)) + 1)
    call c_matrix_init(a(3), c_matrix_memsize(a(1)) + c_matrix_memsize(a(2)) + 1)
    allocate (W(SUM(c_matrix_memsize(a)) + SUM(c_matrix_worksize(a))))
    W(:) = 999
    call c_matrix_eval(a(1), ms(1), X, Y, W)
    call c_matrix_eval(a(2), ms(2), X, Y, W)
    call c_matrix_eval(a(3), ms(3), X, Y, W)
    print'(10f9.3)',W
!
  end subroutine test0
!
! subroutine test1()
!   integer, parameter         :: s = 2
!   integer, parameter         :: m = 5
!   integer, parameter         :: n = 7
!   integer, parameter         :: g = 5
!   integer, parameter         :: mn = m * n
!   integer, parameter         :: swp(m, s - 1) = RESHAPE([1, 2, 3, 4, 5], [m, s - 1])
!   type(mol_block), parameter :: b = mol_block(1, s, m, n, g)
!   type(mol_symmetry)         :: ms
!   type(d_matrix)             :: a
!   real(RK)                   :: X(d, mn), Y(d, mn)
!   real(RK)                   :: LT, LF, LB, H, C(d, d)
!   real(RK), allocatable      :: w(:)
!
!   ms = mol_symmetry(swp)
!   C = 0D0
!   H = 0D0
!   a = d_matrix(1, 28, b)
!   X = sample(d, mn)
!   Y = 0.9D0 * MATMUL(SO3(), X) + 0.1D0 * sample(d, mn)
!   print *, d_matrix_memsize(a)
!   allocate (w(d_matrix_memsize(a)))
!   call d_matrix_eval(a, ms, X, Y, W)
!   call d_matrix_partial_eval(a, 1, 1, 1, [2, 3, 4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call d_matrix_partial_eval(a, 2, 2, 1, [3, 4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call d_matrix_partial_eval(a, 3, 3, 1, [4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call d_matrix_partial_eval(a, 4, 4, 1, [5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call d_matrix_partial_eval(a, 5, 5, 1, [5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   print *
!   call d_matrix_partial_eval(a, 1, 1, 2, [2, 3, 4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call d_matrix_partial_eval(a, 1, 2, 1, [1, 3, 4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call d_matrix_partial_eval(a, 1, 2, 2, [1, 3, 4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!
! end subroutine test1
!
! subroutine test2()
!   integer, parameter    :: l = 3
!   integer, parameter    :: s = 3
!   integer, parameter    :: m = 5
!   integer, parameter    :: n = 7
!   integer, parameter    :: f = 3
!   integer, parameter    :: g = 5
!   integer, parameter    :: mnl = (3 * m + 6) * n
!   integer, parameter    :: swp(m, s - 1) = RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m, s - 1])
!   integer               :: perm(3 * g - 6)
!   type(mol_block_list)  :: blk
!   type(mol_block)       :: b(l)
!   type(mol_symmetry)    :: ms(l)
!   type(d_matrix_list)   :: a
!   real(RK)              :: X(d, mnl), Y(d, mnl)
!   real(RK), allocatable :: w(:)
!   real(RK)              :: C(d * d), H, LT, LF, LB
!   integer               :: i, j, k
!
!   do concurrent(i=1:l)
!     ms(i) = mol_symmetry(swp(:, :i - 1))
!   end do
!
!   do i = 1, l
!     b(i) = mol_block(0, ms(i)%n_sym() + 1, m + i, n, g - i)
!   end do
!
!   k = 0
!   do j = 1, l
!     do i = 1, b(j)%g
!       k = k + 1
!       perm(k) = i
!     end do
!   end do
!   print'(10I4)',perm
!
!   blk = mol_block_list(l, b)
!   a = d_matrix_list(blk, 1)
!
!   X = sample(d, mnl)
!   Y = 0.99D0 * MATMUL(SO3(), X) + 0.01D0 * sample(d, mnl)
!   print *, a%memsize()
!   print *, d_matrix_memsize(a%m)
!
!   allocate (w(a%memsize()))
!
!   call a%eval(ms, X, Y, W)
!   print *, W(1)
!   print'(3f9.3)', W(2:10)
!   print'(*(f9.3))', W(a%o:a%o+l-1)
!   print*
!
!   H = W(a%h)
!   C = W(a%c:a%c + DD - 1)
!   call a%partial_eval(1, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(2, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(3, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(4, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(5, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(6, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(7, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(8, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(9, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   print*
!
!   H = W(a%h)
!   C = W(a%c:a%c + DD - 1)
!   call a%partial_eval(1, perm, 4, 1, W, LT, H, C, LF, LB)
!   perm(:4) = [4,1,2,3]
!   print'(3f9.3)',H, LT
!   call a%partial_eval(2, perm, 3, 1, W, LT, H, C, LF, LB)
!   perm(:4) = [4,3,1,2]
!   print'(3f9.3)',H, LT
!   call a%partial_eval(3, perm, 2, 1, W, LT, H, C, LF, LB)
!   perm(:4) = [4,3,2,1]
!   print'(3f9.3)',H, LT
!   call a%partial_eval(4, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(5, perm, 3, 1, W, LT, H, C, LF, LB)
!   perm(5:7) = [3,1,2]
!   print'(3f9.3)',H, LT
!   call a%partial_eval(6, perm, 2, 1, W, LT, H, C, LF, LB)
!   perm(5:7) = [3,2,1]
!   print'(3f9.3)',H, LT
!   call a%partial_eval(7, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!   call a%partial_eval(8, perm, 2, 1, W, LT, H, C, LF, LB)
!   perm(8:9) = [2,1]
!   print'(3f9.3)',H, LT
!   call a%partial_eval(9, perm, 1, 1, W, LT, H, C, LF, LB)
!   print'(3f9.3)',H, LT
!
! end subroutine test2
!
! subroutine test3()
!   integer, parameter    :: s = 1
!   integer, parameter    :: m = 5, n = 5, g = 3
!   integer, parameter    :: mn = m * n
!   type(mol_block)       :: b = mol_block(0, 2, m, n, g)
!   type(d_matrix_list)   :: dm
!   type(mol_block_list)  :: blk
!   type(mol_symmetry)    :: ms(s)
!   real(RK)              :: X(d, mn), Y(d, mn)
!   real(RK)              :: R1(6, 2, 2, 2), R2(6, 2, 2, 2)
!   real(RK), allocatable :: W(:)
!   integer               :: i, j, k
!
!   ms(1) = mol_symmetry(RESHAPE([2, 1, 3, 4, 5], [m, 1]))
!   blk = mol_block_list(s, [b])
!   dm = d_matrix_list(blk, 1)
!   allocate (w(dm%memsize()))
!   W(:)=999
!
!   X = sample(d, mn)
!   Y = sample(d, mn)
!
!   call dm%eval(ms, X, Y, W)
!
!   do k = 1, 2
!   do j = 1, 2
!   do i = 1, 2
!     R1(1, i, j, k) = sd(d, X, swp(d, m, n, [1, 2, 3], [i, j, k] - 1, ms(1), Y))
!     R1(2, i, j, k) = sd(d, X, swp(d, m, n, [1, 3, 2], [i, j, k] - 1, ms(1), Y))
!     R1(3, i, j, k) = sd(d, X, swp(d, m, n, [2, 1, 3], [i, j, k] - 1, ms(1), Y))
!     R1(4, i, j, k) = sd(d, X, swp(d, m, n, [2, 3, 1], [i, j, k] - 1, ms(1), Y))
!     R1(5, i, j, k) = sd(d, X, swp(d, m, n, [3, 1, 2], [i, j, k] - 1, ms(1), Y))
!     R1(6, i, j, k) = sd(d, X, swp(d, m, n, [3, 2, 1], [i, j, k] - 1, ms(1), Y))
!     R2(1, i, j, k) = pe(dm, d, 3, [1, 2, 3, 1, 2, 3, 1, 2, 3], [1, 1, 1], [i, j, k] - 1, W)
!     R2(2, i, j, k) = pe(dm, d, 3, [1, 2, 3, 1, 2, 3, 1, 3, 2], [1, 2, 1], [i, j, k] - 1, W)
!     R2(3, i, j, k) = pe(dm, d, 3, [1, 2, 3, 2, 1, 3, 2, 1, 3], [2, 1, 1], [i, j, k] - 1, W)
!     R2(4, i, j, k) = pe(dm, d, 3, [1, 2, 3, 2, 1, 3, 2, 3, 1], [2, 2, 1], [i, j, k] - 1, W)
!     R2(5, i, j, k) = pe(dm, d, 3, [1, 2, 3, 3, 1, 2, 3, 1, 2], [3, 1, 1], [i, j, k] - 1, W)
!     R2(6, i, j, k) = pe(dm, d, 3, [1, 2, 3, 3, 1, 2, 3, 2, 1], [3, 2, 1], [i, j, k] - 1, W)
!   end do
!   end do
!   end do
!   call u%assert_almost_equal([R1 - R2], ZERO, 'R1 - R2')
!
! end subroutine test3
!
! function pe(dm, d, l, iper, jper, isym, W) result(res)
!   type(d_matrix_list), intent(in) :: dm
!   integer, intent(in)             :: d, l, iper(l, l), jper(l), isym(l)
!   real(RK), intent(in)            :: W(*)
!   real(RK)                        :: C(d,d), H, T, LF, LB, res
!   integer                         :: i
!   H = W(dm%h)
!   call copy(d * d, W(dm%c), C)
!   do i = 1, l
!     call dm%partial_eval(i, iper(1, i), jper(i), isym(i), W, T, H, C, LF=LF, LB=LB)
!   end do
!   res = LF
!   !res = T
! end function pe
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
! function SO3() result(res)
!   real(RK) :: a(3), res(3, 3)
!   call RANDOM_NUMBER(a)
!   a = a / SQRT(DOT_PRODUCT(a, a))
!   res(:, 1) = [a(1) * a(1), a(1) * a(2) - a(3), a(1) * a(3) + a(2)]
!   res(:, 2) = [a(1) * a(2) + a(3), a(2) * a(2), a(2) * a(3) - a(1)]
!   res(:, 3) = [a(1) * a(3) - a(2), a(2) * a(3) + a(1), a(3) * a(3)]
! end function SO3
!
! pure function swp(d, m, n, per, sym, ms, X) result(res)
!   integer(IK), intent(in) :: d, m, n, per(:), sym(:)
!   type(mol_symmetry), intent(in) :: ms
!   real(RK), intent(in)    :: X(d, m, n)
!   real(RK)                :: tmp(d, m, n), res(d, m * n)
!   integer(IK)             :: i
!   tmp = X
!   do i = 1, SIZE(per)
!     tmp(:, :, per(i)) = X(:, :, i)
!     call ms%swap(d, tmp(:, :, per(i)), sym(i))
!   end do
!   res = RESHAPE(tmp, [d, m * n])
! end function swp
!
! pure function sd(d, X, Y) result(res)
!   integer(IK), intent(in) :: d
!   real(RK), intent(in)    :: X(:, :), Y(:, :)
!   real(RK)                :: C(d, d), R(d, d), W(100), res
!   C = MATMUL(Y, TRANSPOSE(X))
!   call estimate_rotation_matrix(SUM(X * X) + SUM(Y * Y), C, R, W)
!   res = SUM(X**2) + SUM(Y**2) - 2 * SUM(C * R)
! end function sd
!
! pure subroutine copy(d, source, dest)
!   integer(IK), intent(in) :: d
!   real(RK), intent(in)    :: source(*)
!   real(RK), intent(inout) :: dest(*)
!   integer(IK)             :: i
!   do concurrent(i=1:d)
!     dest(i) = source(i)
!   end do
! end subroutine copy
!
end program main

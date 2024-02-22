program main
  use mod_params, only: D, setup_dimension, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_mol_symmetry
  use mod_rotation_matrix
  use mod_lowerbound
  use mod_bb_manager
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test lowerbound')
!
  call setup_dimension(3)
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
    real(RK)         :: GC(10, 15), S(15), W(100), G, C(9)
    type(mol_block)  :: b
    type(bb_manager) :: bm
    integer(IK)      :: i
!
    b = mol_block(1, 8, 3, 5)
!
    do i = 1, 15
      GC(:, i) = gcov(D, 8)
      call estimate_sdmin(GC(1, i), GC(2, i), W)
      S(i) = W(1)
    end do
!
    G = ZERO
    C = ZERO
!
    call lowerbound(0, b, G, C, S, W)
    print*,0, W(1)
    G = G + GC(1, 1)
    C = C + GC(2:, 1)
    call lowerbound(1, b, G, C, S([7,8,9,10,12,13,14,15]), W)
    print*,1, W(1)
    G = G + GC(1, 7)
    C = C + GC(2:, 7)
    call lowerbound(2, b, G, C, S([13,14,15]), W)
    print*,2, W(1)
    G = G + GC(1, 13)
    C = C + GC(2:, 13)
    call lowerbound(3, b, G, C, S, W)
    print*,3, W(1)
!
    G = ZERO
    C = ZERO
    G = G + GC(1, 3)
    C = C + GC(2:, 3)
    call lowerbound(1, b, G, C, S([6,7,9,10,11,12,14,15]), W)
    print*,1, W(1)
    G = G + GC(1, 9)
    C = C + GC(2:, 9)
    call lowerbound(2, b, G, C, S([11,12,15]), W)
    print*,2, W(1)
    G = G + GC(1, 15)
    C = C + GC(2:, 15)
    call lowerbound(3, b, G, C, S, W)
    print*,3, W(1)
!
  end subroutine test0
!
! subroutine test1()
!   type(mol_block)       :: b
!   type(mol_symmetry)    :: ms
!   type(c_matrix)        :: cm
!   real(RK)              :: G, C(DD), R(5, 3)
!   real(RK)              :: X(D, 8 * 3), Y(D, 8 * 5)
!   real(RK), allocatable :: w(:)
!
!   b = mol_block(2, 8, 3, 5)
!   ms = mol_symmetry(RESHAPE([2, 3, 1, 4, 5, 6, 7, 8], [8, 1]))
!   cm = c_matrix(b)
!   allocate (W(memsize_c_matrix(cm) + worksize_c_matrix(cm)))
!   W(:) = 999
!   call c_matrix_eval(cm, b, ms, X, Y, W)
!   G = 0D0
!   C = 0D0
!
!   a = s_matrix(1, 28, b)
!   X = sample(d, mn)
!   Y = 0.9D0 * MATMUL(SO3(), X) + 0.1D0 * sample(d, mn)
!   print *, s_matrix_memsize(a)
!   allocate (w(s_matrix_memsize(a)))
!   call s_matrix_eval(a, ms, X, Y, W)
!   call s_matrix_partial_eval(a, 1, 1, 1, [2, 3, 4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call s_matrix_partial_eval(a, 2, 2, 1, [3, 4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call s_matrix_partial_eval(a, 3, 3, 1, [4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call s_matrix_partial_eval(a, 4, 4, 1, [5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call s_matrix_partial_eval(a, 5, 5, 1, [5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   print *
!   call s_matrix_partial_eval(a, 1, 1, 2, [2, 3, 4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call s_matrix_partial_eval(a, 1, 2, 1, [1, 3, 4, 5], W, LT, H, C, LF, LB)
!   print'(4f9.3)', LF, LB, LF + LB, LT
!   call s_matrix_partial_eval(a, 1, 2, 2, [1, 3, 4, 5], W, LT, H, C, LF, LB)
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
!   type(s_matrix_list)   :: a
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
!   a = s_matrix_list(blk, 1)
!
!   X = sample(d, mnl)
!   Y = 0.99D0 * MATMUL(SO3(), X) + 0.01D0 * sample(d, mnl)
!   print *, a%memsize()
!   print *, s_matrix_memsize(a%m)
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
!   type(s_matrix_list)   :: dm
!   type(mol_block_list)  :: blk
!   type(mol_symmetry)    :: ms(s)
!   real(RK)              :: X(d, mn), Y(d, mn)
!   real(RK)              :: R1(6, 2, 2, 2), R2(6, 2, 2, 2)
!   real(RK), allocatable :: W(:)
!   integer               :: i, j, k
!
!   ms(1) = mol_symmetry(RESHAPE([2, 1, 3, 4, 5], [m, 1]))
!   blk = mol_block_list(s, [b])
!   dm = s_matrix_list(blk, 1)
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
!   type(s_matrix_list), intent(in) :: dm
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
end program main

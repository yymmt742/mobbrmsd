program main
  use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_unittest
  use mod_tree
  implicit none
  type(unittest) :: u
#ifdef USE_REAL32
  integer(IK), parameter :: place = 3
#else
  integer(IK), parameter :: place = 7
#endif
!
  call u%init('test tree')
  call test1()
  call test2()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    type(tree)            :: t
    real(RK), allocatable :: w(:)
    integer               :: s(4)
    integer(IK)           :: i
!
!   0    1  o
!           |_____________
!           | | | | | | | |
!   1    8  o o o o o o o o
!           |_________
!           | | | | | |
!   2   48  o o o o o o
!           |_____
!           | | | |
!   3  192  o o o o
!           |_
!           | |
!   4  384  o o
!     ----
!      633
!
!   mmap      1         2         3         4
!    1234567890123456789012345678901234567890
!   [1_1_1_1_1_1_1_1_2_2_2_2_2_2_3_3_3_3_4_4_]
!
    t = tree(4, 2)
!
    call u%assert_equal(tree_n_depth(t%q), 4, 'n_depth    [4,3,2,1]')
    call u%assert_equal(tree_nnodes(t%q), 20, 'nnodes     [4,3,2,1]')
    call u%assert_equal(NINT(EXP(tree_log_ncomb(t%q))), 632, 'log_ncomb  [4,3,2,1]')
    call u%assert_almost_equal(tree_ncomb_frac(t%q), 6.32_RK, 'ncomb_frac [4,3,2,1]', place=place)
    call u%assert_equal(tree_ncomb_exp(t%q), 2, 'ncomb_exp  [4,3,2,1]')
!
    allocate (W(2 * tree_nnodes(t%q)))
    do i = 1, SIZE(w)
      w(i) = -i
    end do
!
    call tree_reset(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
!
    call u%assert_equal(tree_queue_pointer(t%q, t%s), 1, 'tree_queue_pointer  ')
    call u%assert_equal(tree_current_pointer(t%q, t%s), 8, 'tree_current_pointer')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0), 1, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 1), 2, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 0), 3, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 1), 4, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 2, 0), 5, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 2, 1), 6, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 3, 0), 7, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 3, 1), 8, 'tree_node_pointer   ')
!
    call tree_descend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
!
    call u%assert_equal(tree_queue_pointer(t%q, t%s), 9, 'tree_queue_pointer  ')
    call u%assert_equal(tree_current_pointer(t%q, t%s), 14, 'tree_current_pointer')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0), 9, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 1), 10, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 0), 11, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 1), 12, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 2, 0), 13, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 2, 1), 14, 'tree_node_pointer   ')
!
    call tree_descend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
!
    call u%assert_equal(tree_queue_pointer(t%q, t%s), 15, 'tree_queue_pointer  ')
    call u%assert_equal(tree_current_pointer(t%q, t%s), 18, 'tree_current_pointer')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0), 15, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 1), 16, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 0), 17, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 1), 18, 'tree_node_pointer   ')
!
    call tree_descend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
!
    call u%assert_equal(tree_queue_pointer(t%q, t%s), 19, 'tree_queue_pointer  ')
    call u%assert_equal(tree_current_pointer(t%q, t%s), 20, 'tree_current_pointer')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0), 19, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 1), 20, 'tree_node_pointer   ')
!
    call tree_current_sequence(t%q, t%s, s)
    call u%assert_equal(s, [7, 5, 3, 1], 'current_sequence    ')
    call tree_current_permutation(t%q, t%s, s)
    call u%assert_equal(s, [4, 3, 2, 1], 'current_permutation ')
    call tree_current_mapping(t%q, t%s, s)
    call u%assert_equal(s, [1, 1, 1, 1], 'current_mapping     ')
!
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
    call tree_current_sequence(t%q, t%s, s)
    call u%assert_equal(s, [7, 5, 3, 0], 'current_sequence    ')
    call tree_current_permutation(t%q, t%s, s)
    call u%assert_equal(s, [4, 3, 2, 1], 'current_permutation ')
    call tree_current_mapping(t%q, t%s, s)
    call u%assert_equal(s, [1, 1, 1, 0], 'current_mapping     ')
!
    call tree_ascend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
    call tree_descend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
    call tree_current_sequence(t%q, t%s, s)
    call u%assert_equal(s, [7, 5, 2, 1], 'current_sequence    ')
    call tree_current_permutation(t%q, t%s, s)
    call u%assert_equal(s, [4, 3, 2, 1], 'current_permutation ')
    call tree_current_mapping(t%q, t%s, s)
    call u%assert_equal(s, [1, 1, 0, 1], 'current_mapping     ')
!
    call tree_ascend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
    call tree_descend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
    call tree_current_sequence(t%q, t%s, s)
    call u%assert_equal(s, [7, 5, 1, 1], 'current_sequence    ')
    call tree_current_permutation(t%q, t%s, s)
    call u%assert_equal(s, [4, 3, 1, 2], 'current_permutation ')
    call tree_current_mapping(t%q, t%s, s)
    call u%assert_equal(s, [1, 1, 1, 1], 'current_mapping     ')
!
    call tree_ascend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
    call tree_descend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
    call tree_current_sequence(t%q, t%s, s)
    call u%assert_equal(s, [7, 5, 0, 1], 'current_sequence    ')
    call tree_current_permutation(t%q, t%s, s)
    call u%assert_equal(s, [4, 3, 1, 2], 'current_permutation ')
    call tree_current_mapping(t%q, t%s, s)
    call u%assert_equal(s, [1, 1, 0, 1], 'current_mapping     ')
!
    call tree_ascend(t%q, t%s)
    call tree_ascend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
    call tree_descend(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 2, 999._RK, w)
    call tree_current_sequence(t%q, t%s, s)
    call u%assert_equal(s, [7, 4, 3, 1], 'current_sequence    ')
    call tree_current_permutation(t%q, t%s, s)
    call u%assert_equal(s, [4, 3, 2, 1], 'current_permutation ')
    call tree_current_mapping(t%q, t%s, s)
    call u%assert_equal(s, [1, 0, 1, 1], 'current_mapping     ')
!
  end subroutine test1
!
  subroutine test2()
    type(tree)            :: t
    real(RK)              :: ub
    real(RK), allocatable :: w(:)
    integer(IK)           :: prm(4)
    integer(IK)           :: p, q, r
!
    p = 1
    t = tree(4, 1)
    ub = 999.0_RK
    allocate (W(tree_nnodes(t%q)))
    call RANDOM_NUMBER(w)
!
    do
      do
        call tree_select_top_node(t%q, t%s, 1, ub, w)
        if (tree_queue_is_bottom(t%q, t%s) .or. tree_queue_is_empty(t%q, t%s)) exit
        r = tree_current_pointer(t%q, t%s)
        call tree_descend(t%q, t%s)
        p = p + 1
        q = tree_queue_pointer(t%q, t%s)
        call RANDOM_NUMBER(w(q:q + 4 - p))
        w(q:q + 4 - p) = w(q:q + 4 - p) * 1.0 + w(r)
      end do
      if (tree_queue_is_bottom(t%q, t%s) .and. tree_queue_is_selected(t%q, t%s)) then
        call tree_current_permutation(t%q, t%s, prm)
        print '(2f9.3,3L4,*(I4))', ub, w(tree_current_pointer(t%q, t%s)), tree_queue_is_empty(t%q, t%s),&
       &         tree_queue_is_explored(t%q, t%s), tree_queue_is_unexplored(t%q, t%s),&
       &         prm, tree_current_pointer(t%q, t%s)
        ub = MIN(ub, w(tree_current_pointer(t%q, t%s)))
      end if
      do
        if (tree_queue_is_left(t%q, t%s, 1, ub, W) .or. tree_queue_is_root(t%q, t%s)) exit
        call tree_ascend(t%q, t%s)
        p = p - 1
      end do
      if (tree_queue_is_empty(t%q, t%s) .and. tree_queue_is_root(t%q, t%s)) exit
    end do
!
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
!
  end subroutine test2
!
! subroutine dump(n, w)
!   integer(IK), intent(in) :: n
!   real(RK), intent(in)    :: w(*)
!   integer(IK)             :: i, p
!   p = 1
!   do i = 1, n
!     print'(*(F9.3))', w(p:p + n - i)
!     p = p + n - i + 1
!   end do
! end subroutine dump
!
end program main


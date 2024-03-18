program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_unittest
  use mod_tree
  implicit none
  type(unittest) :: u
!
  call u%init('test tree')
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    type(mol_block)       :: b
    type(tree)            :: t
    real(RK), allocatable :: w(:)
    integer(IK)           :: i
!
!   0    1  o____________________________
!           |                | | | | | | |
!   1    8  o_______________ o o o o o o o
!           |       | | | | |
!   2   48  o______ o o o o o
!           |  | | |
!   3  192  o_ o o o
!           | |
!   4  384  o o
!     ----
!      633
!
!   mmap      1         2         3         4         5         6
!    123456789012345678901234567890123456789012345678901234567890
!   [1_1_1_1_1_1_1_1_2__2__2__2__2__2__3___3___3___3___4____4____]
!
    b = mol_block(10, 4, sym=RESHAPE([1, 2, 4, 3], [4, 1]))
    t = tree(b%q, memsize1)
!
    call u%assert_equal(tree_n_depth(t%q),                 4, 'n_depth    [4,3,2,1]')
    call u%assert_equal(tree_memsize(t%q),                60, 'memsize    [4,3,2,1]')
    call u%assert_equal(NINT(EXP(tree_log_ncomb(t%q))),  632, 'log_ncomb  [4,3,2,1]')
    call u%assert_almost_equal(tree_ncomb_frac(t%q), 6.32_RK, 'ncomb_frac [4,3,2,1]')
    call u%assert_equal(tree_ncomb_exp(t%q),               2, 'ncomb_exp  [4,3,2,1]')
!
    call u%assert_equal(tree_queue_pointer(t%q, t%s),      1, 'tree_queue_pointer  ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0), 1, 'tree_node_pointer   ')
!
    allocate(W(tree_memsize(t%q)))
    do i = 1, SIZE(w)
      w(i) = -i
    end do
!
    call tree_reset(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
!
    call u%assert_equal(tree_queue_pointer(t%q, t%s),      1, 'tree_queue_pointer  ')
    call u%assert_equal(tree_current_pointer(t%q, t%s),   15, 'tree_current_pointer')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0), 1, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 1), 3, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 0), 5, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 1), 7, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 2, 0), 9, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 2, 1),11, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 3, 0),13, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 3, 1),15, 'tree_node_pointer   ')
!
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
!
    call u%assert_equal(tree_queue_pointer(t%q, t%s),     17, 'tree_queue_pointer  ')
    call u%assert_equal(tree_current_pointer(t%q, t%s),   32, 'tree_current_pointer')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0),17, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 1),20, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 0),23, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 1),26, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 2, 0),29, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 2, 1),32, 'tree_node_pointer   ')
!
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
!
    call u%assert_equal(tree_queue_pointer(t%q, t%s),     35, 'tree_queue_pointer  ')
    call u%assert_equal(tree_current_pointer(t%q, t%s),   47, 'tree_current_pointer')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0),35, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 1),39, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 0),43, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 1, 1),47, 'tree_node_pointer   ')
!
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
!
    call u%assert_equal(tree_queue_pointer(t%q, t%s),     51, 'tree_queue_pointer  ')
    call u%assert_equal(tree_current_pointer(t%q, t%s),   56, 'tree_current_pointer')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0),51, 'tree_node_pointer   ')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 1),56, 'tree_node_pointer   ')
!
    call u%assert_equal(tree_current_sequence(t%q, t%s),    [8,6,4,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q, t%s), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q, t%s),     [1,1,1,1], 'current_mapping     ')
!
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call u%assert_equal(tree_current_sequence(t%q, t%s),    [8,6,4,1], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q, t%s), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q, t%s),     [1,1,1,0], 'current_mapping     ')
!
    call tree_leave(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call u%assert_equal(tree_current_sequence(t%q, t%s),    [8,6,3,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q, t%s), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q, t%s),     [1,1,0,1], 'current_mapping     ')
!
    call tree_leave(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call u%assert_equal(tree_current_sequence(t%q, t%s),    [8,6,2,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q, t%s), [4,3,1,2], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q, t%s),     [1,1,1,1], 'current_mapping     ')
!
    call tree_leave(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call u%assert_equal(tree_current_sequence(t%q, t%s),    [8,6,1,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q, t%s), [4,3,1,2], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q, t%s),     [1,1,0,1], 'current_mapping     ')
!
    call tree_leave(t%q, t%s)
    call tree_leave(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call u%assert_equal(tree_current_sequence(t%q, t%s),    [8,5,4,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q, t%s), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q, t%s),     [1,0,1,1], 'current_mapping     ')
!
  end subroutine test1
!
  pure function memsize1(b, p) result(res)
    integer(IK), intent(in)     :: b(*)
    integer(IK), intent(in)     :: p
    integer(IK)                 :: res
    res = 1 + p
  end function memsize1
!
end program main


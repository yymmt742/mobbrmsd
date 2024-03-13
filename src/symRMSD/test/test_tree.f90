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
!   -1234567890123456789012345678901234567890123456789012345678901
!   [01-1-1-1-1-1-1-1-2--2--2--2--2--2--3---3---3---3---4----4----]
!
    b = mol_block(10, 4, sym=RESHAPE([1, 2, 4, 3], [4, 1]))
    print*,mol_block_nsym(b%q)
    t = tree(b%q, memsize1)
    print'(3I3)', t%q(:2)
    print'(3I3)', t%q(3:)
    print'(2I3)', t%s(1)
    print'(2I3)', t%s(2:)
!
    call u%assert_equal(tree_n_depth(t%q),                 4, 'n_depth    [4,3,2,1]')
    call u%assert_equal(tree_memsize(t%q),                61, 'memsize    [4,3,2,1]')
    call u%assert_equal(NINT(EXP(tree_log_ncomb(t%q))),  633, 'log_ncomb  [4,3,2,1]')
    call u%assert_almost_equal(tree_ncomb_frac(t%q), 6.33_RK, 'ncomb_frac [4,3,2,1]')
    call u%assert_equal(tree_ncomb_exp(t%q),               2, 'ncomb_exp  [4,3,2,1]')
!
    call u%assert_equal(tree_current_level(t%s),           0, 'tree_current_level  ')
    call u%assert_equal(tree_root_pointer(t%q),            1, 'tree_root_pointer   ')
    call u%assert_equal(tree_queue_pointer(t%q, t%s),      1, 'tree_queue_pointer  ')
    call u%assert_equal(tree_current_pointer(t%q, t%s),    1, 'tree_current_pointer')
    call u%assert_equal(tree_node_pointer(t%q, t%s, 0, 0), 1, 'tree_node_pointer   ')
    print *, tree_current_sequence(t%q, t%s)
    print *, tree_current_permutation(t%q, t%s)
    print *, tree_current_mapping(t%q, t%s)
    print*, tree_queue_is_empty(t%q, t%s), tree_queue_is_bottom(t%q, t%s)
!
    allocate(W(tree_memsize(t%q)))
    do i = 1, SIZE(w)
      w(i) = -i
    end do
!
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    call tree_expand(t%q, t%s)
    call tree_select_top_node(t%q, t%s, 999._RK, w)
    print'(2I3)', t%s(1)
    print'(2I3)', t%s(2:)
    print *, tree_current_level(t%s)
    print *, tree_root_pointer(t%q)
    print *, tree_queue_pointer(t%q, t%s)
    print *, tree_current_pointer(t%q, t%s)
    print *, tree_node_pointer(t%q, t%s, 0, 0)
    print *, tree_current_sequence(t%q, t%s)
    print *, tree_current_permutation(t%q, t%s)
    print *, tree_current_mapping(t%q, t%s)
    print*, tree_queue_is_empty(t%q, t%s), tree_queue_is_bottom(t%q, t%s)
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


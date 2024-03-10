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
    type(mol_block) :: b
    type(tree)      :: t
    integer(IK)     :: i
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
    print'(4I3)', t%q(:2)
    print'(4I3)', t%q(3:)
!
    call u%assert_equal(tree_memsize(t%q),                61, 'memsize             ')
    call u%assert_equal(NINT(EXP(tree_log_ncomb(t%q))),  633, 'log_ncomb  [4,3,2,1]')
    call u%assert_almost_equal(tree_ncomb_frac(t%q), 6.33_RK, 'ncomb_frac [4,3,2,1]')
    call u%assert_equal(tree_ncomb_exp(t%q),               2, 'ncomb_exp  [4,3,2,1]')
!
  print *, tree_n_depth(t%q)
  print *, tree_current_level(t%q)
  print *, tree_root_pointer(t%q)
  print *, tree_queue_pointer(t%q)
  print *, tree_current_pointer(t%q)
  print *, tree_node_pointer(t%q, 0, 0)
  print *, tree_current_sequence(t%q)
  print*, tree_queue_is_empty(t%q), tree_queue_is_bottom(t%q)
!
!   allocate(W(tree_memsize(t%t, t%q)))
    do i = 1, SIZE(t%x)
      t%x(i) = -i
    end do
!
    call tree_expand(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
print'(4I3)', t%q(:2)
print'(4I3)', t%q(3:)
  print *, tree_n_depth(t%q)
  print *, tree_current_level(t%q)
  print *, tree_root_pointer(t%q)
  print *, tree_queue_pointer(t%q)
  print *, tree_current_pointer(t%q)
  print *, tree_node_pointer(t%q, 0, 0)
  print *, tree_current_sequence(t%q)
  print*, tree_queue_is_empty(t%q), tree_queue_is_bottom(t%q)
!
    call tree_expand(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call tree_expand(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call tree_expand(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
print'(4I3)', t%q(:2)
print'(4I3)', t%q(3:)
  print *, tree_n_depth(t%q)
  print *, tree_current_level(t%q)
  print *, tree_root_pointer(t%q)
  print *, tree_queue_pointer(t%q)
  print *, tree_current_pointer(t%q)
  print *, tree_node_pointer(t%q, 0, 0)
  print *, tree_current_sequence(t%q)
  print *, tree_current_permutation(t%q)
  print *, tree_current_mapping(t%q)
  print*, tree_queue_is_empty(t%q), tree_queue_is_bottom(t%q)
!
    call u%assert_equal(tree_current_sequence(t%q),    [8,6,4,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q),     [1,1,1,1], 'current_mapping     ')
!
    call tree_select_top_node(t%q, 999._RK, t%x)
    call u%assert_equal(tree_current_sequence(t%q),    [8,6,4,1], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q),     [1,1,1,0], 'current_mapping     ')
!
    call tree_leave(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call tree_expand(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call u%assert_equal(tree_current_sequence(t%q),    [8,6,3,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q),     [1,1,0,1], 'current_mapping     ')
!
    call tree_leave(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call tree_expand(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call u%assert_equal(tree_current_sequence(t%q),    [8,6,2,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q), [4,3,1,2], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q),     [1,1,1,1], 'current_mapping     ')
!
    call tree_leave(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call tree_expand(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call u%assert_equal(tree_current_sequence(t%q),    [8,6,1,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q), [4,3,1,2], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q),     [1,1,0,1], 'current_mapping     ')
!
    call tree_leave(t%q)
    call tree_leave(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call tree_expand(t%q)
    call tree_select_top_node(t%q, 999._RK, t%x)
    call u%assert_equal(tree_current_sequence(t%q),    [8,5,4,2], 'current_sequence    ')
    call u%assert_equal(tree_current_permutation(t%q), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(tree_current_mapping(t%q),     [1,0,1,1], 'current_mapping     ')
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


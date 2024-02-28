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
    real(RK), allocatable :: W(:)
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
    b = mol_block(2, 10, 4, 6)
    t = tree(b, memsize1)
    t%p = 1
!
    call u%assert_equal(t%memsize(),                    61, 'memsize             ')
    call u%assert_equal(NINT(EXP(t%log_ncomb())),      633, 'log_ncomb  [4,3,2,1]')
    call u%assert_almost_equal(t%ncomb_frac(),     6.33_RK, 'ncomb_frac [4,3,2,1]')
    call u%assert_equal(t%ncomb_exp(),                   2, 'ncomb_exp  [4,3,2,1]')
!
    allocate(W(t%memsize()))
    do i = 1, SIZE(W)
      W(i) = -i
    end do
!
    call t%expand()
    call t%select_top_node(999._RK, W)
    call t%expand()
    call t%select_top_node(999._RK, W)
    call t%expand()
    call t%select_top_node(999._RK, W)
    call t%expand()
    call t%select_top_node(999._RK, W)
!
    call u%assert_equal(t%current_sequence(),    [8,6,4,2], 'current_sequence    ')
    call u%assert_equal(t%current_permutation(), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(t%current_mapping(),     [1,1,1,1], 'current_mapping     ')
!
    call t%select_top_node(999._RK, W)
    call u%assert_equal(t%current_sequence(),    [8,6,4,1], 'current_sequence    ')
    call u%assert_equal(t%current_permutation(), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(t%current_mapping(),     [1,1,1,0], 'current_mapping     ')
!
    call t%leave()
    call t%select_top_node(999._RK, W)
    call t%expand()
    call t%select_top_node(999._RK, W)
    call u%assert_equal(t%current_sequence(),    [8,6,3,2], 'current_sequence    ')
    call u%assert_equal(t%current_permutation(), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(t%current_mapping(),     [1,1,0,1], 'current_mapping     ')
!
    call t%leave()
    call t%select_top_node(999._RK, W)
    call t%expand()
    call t%select_top_node(999._RK, W)
    call u%assert_equal(t%current_sequence(),    [8,6,2,2], 'current_sequence    ')
    call u%assert_equal(t%current_permutation(), [4,3,1,2], 'current_permutation ')
    call u%assert_equal(t%current_mapping(),     [1,1,1,1], 'current_mapping     ')
!
    call t%leave()
    call t%select_top_node(999._RK, W)
    call t%expand()
    call t%select_top_node(999._RK, W)
    call u%assert_equal(t%current_sequence(),    [8,6,1,2], 'current_sequence    ')
    call u%assert_equal(t%current_permutation(), [4,3,1,2], 'current_permutation ')
    call u%assert_equal(t%current_mapping(),     [1,1,0,1], 'current_mapping     ')
!
    call t%leave()
    call t%leave()
    call t%select_top_node(999._RK, W)
    call t%expand()
    call t%select_top_node(999._RK, W)
    call t%expand()
    call t%select_top_node(999._RK, W)
    call u%assert_equal(t%current_sequence(),    [8,5,4,2], 'current_sequence    ')
    call u%assert_equal(t%current_permutation(), [4,3,2,1], 'current_permutation ')
    call u%assert_equal(t%current_mapping(),     [1,0,1,1], 'current_mapping     ')
!
!   call u%assert_equal(node_pointer(t, q, 0, 0), [1, 9,15,23], 'node_pointer [0,0]     ')
!   call u%assert_equal(node_pointer(t, q, 0, 1), [2,10,17,27], 'node_pointer [0,1]     ')
!   call u%assert_equal(node_pointer(t, q, 1, 0), [3,11,19,31], 'node_pointer [1,0]     ')
!   call u%assert_equal(node_pointer(t, q, 1, 1), [4,12,21,35], 'node_pointer [1,1]     ')
!   call u%assert_equal(node_pointer(t, q, 2, 0), [5,13,23,39], 'node_pointer [2,0]     ')
!   call u%assert_equal(node_pointer(t, q, 2, 1), [6,14,25,43], 'node_pointer [2,1]     ')
!   call u%assert_equal(node_pointer(t, q, 3, 0), [7,15,27,47], 'node_pointer [3,0]     ')
!   call u%assert_equal(node_pointer(t, q, 3, 1), [8,16,29,51], 'node_pointer [3,1]     ')
!
  end subroutine test1
!
  pure function memsize1(b, p) result(res)
    type(mol_block), intent(in) :: b
    integer(IK), intent(in)     :: p
    integer(IK)                 :: res
    res = 1 + p
  end function memsize1
!
  pure function memsize2(b, p) result(res)
    type(mol_block), intent(in) :: b
    integer(IK), intent(in)     :: p
    integer(IK)                 :: res
    res = 1 + p**2
  end function memsize2
!
end program main


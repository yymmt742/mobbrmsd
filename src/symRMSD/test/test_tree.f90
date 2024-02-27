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
    real(RK), allocatable :: W(:)
    integer(IK)           :: i
!
    b = mol_block(1, 10, 4, 6)
    t = tree(b, memsize)
    print*,NINT(EXP(t%log_ncomb())), t%ncomb_frac(), t%ncomb_exp()
    print*,t%current_permutation()
!
!   call u%assert_equal(NINT(EXP(log_ncomb(q(:1)))),      4, 'a. log_ncomb [4]       ')
!   call u%assert_equal(NINT(EXP(log_ncomb(q(:2)))),     16, 'a. log_ncomb [4,3]     ')
!   call u%assert_equal(NINT(EXP(log_ncomb(q(:3)))),     40, 'a. log_ncomb [4,3,2]   ')
!   call u%assert_equal(NINT(EXP(log_ncomb(q(:4)))),     64, 'a. log_ncomb [4,3,2,1] ')
!   call u%assert_equal(memsize_queue(q),         [8,9,8,5], 'a. memsize_queue       ')
!
!   b = mol_block(2, 10, 6, 4)
!   q = queue([8, 6, 4, 2], [1, 1, 2, 4])
!   t = tree(b)
!             1         2         3
!    123456789012345678901234567890
!   [111111112222223-3-3-3-4---4---]
!
!   call u%assert_equal(NINT(EXP(log_ncomb(q(:1)))),      8, 'b. log_ncomb [8]       ')
!   call u%assert_equal(NINT(EXP(log_ncomb(q(:2)))),     56, 'b. log_ncomb [8,6]     ')
!   call u%assert_equal(NINT(EXP(log_ncomb(q(:3)))),    248, 'b. log_ncomb [8,6,4]   ')
!   call u%assert_equal(NINT(EXP(log_ncomb(q(:4)))),    632, 'b. log_ncomb [8,6,4,2] ')
!   call u%assert_equal(memsize_queue(q),         [8,6,8,8], 'b. memsize_queue       ')
!
!   allocate(W(SUM(memsize_queue(q))))
!   do i = 1, SIZE(W)
!     W(i) = -i
!   end do
!
!   call u%assert_equal(queue_pointer(q),         [1,1,1,1], '1. queue_pointer       ')
!   call setup_queue(t, q)
!   call u%assert_equal(queue_pointer(q),       [1,9,15,23], '1. queue_pointer(setup)')
!   call u%assert_equal(current_sequence(t, q),   [0,0,0,0], '1. current_sequence 0  ')
!   call u%assert_equal(current_permutation(t, q),[1,2,3,4], '1. current_permutation ')
!   call set_top_node(q(1), 999._RK, W, .true.)
!   call set_top_node(q(2), 999._RK, W, .true.)
!   call set_top_node(q(3), 999._RK, W, .true.)
!   call set_top_node(q(4), 999._RK, W, .true.)
!   call u%assert_equal(current_sequence(t, q),   [8,6,4,2], '2. current_sequence    ')
!   call u%assert_equal(current_permutation(t, q),[4,3,2,1], '2. current_permutation ')
!   call u%assert_equal(current_mapping(t, q),    [1,1,1,1], '2. current_mapping     ')
!   call u%assert_equal(current_pointer(q),    [8,14,21,27], '2. current_pointer     ')
!
!   call set_top_node(q(1), 999._RK, W, .false.)
!   call set_top_node(q(2), 999._RK, W, .false.)
!   call set_top_node(q(3), 999._RK, W, .false.)
!   call set_top_node(q(4), 999._RK, W, .false.)
!   call u%assert_equal(current_sequence(t, q),   [7,5,3,1], '3. current_sequence 2  ')
!   call u%assert_equal(current_permutation(t, q),[4,3,2,1], '3. current_permutation ')
!   call u%assert_equal(current_mapping(t, q),    [0,0,0,0], '3. current_mapping     ')
!   call u%assert_equal(current_pointer(q),    [7,13,19,23], '3. current_pointer     ')
!
!   call set_top_node(q(1), 999._RK, W, .false.)
!   call set_top_node(q(2), 999._RK, W, .false.)
!   call set_top_node(q(3), 999._RK, W, .false.)
!   call set_top_node(q(4), 999._RK, W, .false.)
!   call u%assert_equal(current_sequence(t, q),   [6,4,2,0], '4. current_sequence 3  ')
!   call u%assert_equal(current_permutation(t, q),[3,2,1,4], '4. current_permutation ')
!   call u%assert_equal(current_mapping(t, q),    [1,1,1,1], '4. current_mapping     ')
!   call u%assert_equal(current_pointer(q),    [6,12,17,19], '4. current_pointer     ')
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
  pure function memsize(b, p) result(res)
    type(mol_block), intent(in) :: b
    integer(IK), intent(in)     :: p
    integer(IK)                 :: res
    res = 1 + p
  end function memsize
!
end program main


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
    type(queue)     :: q(4)
    real(RK), allocatable :: W(:)
    integer(IK)           :: i
!
    b = mol_block(1, 10, 4, 6)
    q = queue([4, 3, 2, 1], [2, 3, 4, 5])
    t = tree(b)
    call u%assert_equal(NINT(EXP(log_ncomb(q(:1)))),      4, 'a. log_ncomb [4]       ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:2)))),     16, 'a. log_ncomb [4,3]     ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:3)))),     40, 'a. log_ncomb [4,3,2]   ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:4)))),     64, 'a. log_ncomb [4,3,2,1] ')
    call u%assert_equal(memsize_queue(q),         [8,9,8,5], 'a. memsize_queue       ')
!
    b = mol_block(2, 10, 6, 4)
    q = queue([8, 6, 4, 2], [1, 1, 2, 4])
    t = tree(b)
    call u%assert_equal(NINT(EXP(log_ncomb(q(:1)))),      8, 'b. log_ncomb [8]       ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:2)))),     56, 'b. log_ncomb [8,6]     ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:3)))),    248, 'b. log_ncomb [8,6,4]   ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:4)))),    632, 'b. log_ncomb [8,6,4,2] ')
    call u%assert_equal(memsize_queue(q),         [8,6,8,8], 'b. memsize_queue       ')
!
    allocate(W(SUM(memsize_queue(q))))
    do i = 1, SIZE(W)
      W(i) = -i
    end do
!
    call u%assert_equal(queue_pointer(q),         [1,1,1,1], '1. queue_pointer       ')
    call setup_queue(t, q)
    call u%assert_equal(queue_pointer(q),       [1,9,15,23], '1. queue_pointer(setup)')
    call u%assert_equal(current_sequence(t, q),   [0,0,0,0], '1. current_sequence 0  ')
    call u%assert_equal(current_permutation(t, q),[1,2,3,4], '1. current_permutation ')
    call set_top_node(q(1), 999._RK, W, .true.)
    call set_top_node(q(2), 999._RK, W, .true.)
    call set_top_node(q(3), 999._RK, W, .true.)
    call set_top_node(q(4), 999._RK, W, .true.)
    call u%assert_equal(current_sequence(t, q),   [8,6,4,2], '2. current_sequence    ')
    call u%assert_equal(current_permutation(t, q),[4,3,2,1], '2. current_permutation ')
    call u%assert_equal(current_mapping(t, q),    [1,1,1,1], '2. current_mapping     ')
!
    call set_top_node(q(1), 999._RK, W, .false.)
    call set_top_node(q(2), 999._RK, W, .false.)
    call set_top_node(q(3), 999._RK, W, .false.)
    call set_top_node(q(4), 999._RK, W, .false.)
    call u%assert_equal(current_sequence(t, q),   [7,5,3,1], '3. current_sequence 2  ')
    call u%assert_equal(current_permutation(t, q),[4,3,2,1], '3. current_permutation ')
    call u%assert_equal(current_mapping(t, q),    [0,0,0,0], '3. current_mapping     ')
!
    call set_top_node(q(1), 999._RK, W, .false.)
    call set_top_node(q(2), 999._RK, W, .false.)
    call set_top_node(q(3), 999._RK, W, .false.)
    call set_top_node(q(4), 999._RK, W, .false.)
    call u%assert_equal(current_sequence(t, q),   [6,4,2,0], '4. current_sequence 3  ')
    call u%assert_equal(current_permutation(t, q),[3,2,1,4], '4. current_permutation ')
    call u%assert_equal(current_mapping(t, q),    [1,1,1,1], '4. current_mapping     ')
!
  end subroutine test1
!
end program main

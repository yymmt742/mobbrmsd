program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
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
    type(tree)  :: t
    type(queue) :: q(4)
    real(RK), allocatable :: W(:)
    logical, parameter    :: T = .true., F = .false.
    integer(IK)           :: i
!
    q = queue([4, 3, 2, 1])
    t = tree(4, 1, 1)
    call u%assert_equal(NINT(EXP(log_ncomb(q(:1)))),  4, 'log_ncomb [4]       ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:2)))), 16, 'log_ncomb [4,3]     ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:3)))), 40, 'log_ncomb [4,3,2]   ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:4)))), 64, 'log_ncomb [4,3,2,1] ')
    call u%assert_equal(memsize_tree(t),             10, 'memsize_tree        ')
!
    q = queue([8, 6, 4, 2])
    t = tree(4, 2, 1)
    call u%assert_equal(NINT(EXP(log_ncomb(q(:1)))),  8, 'log_ncomb [8]       ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:2)))), 56, 'log_ncomb [8,6]     ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:3)))),248, 'log_ncomb [8,6,4]   ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:4)))),632, 'log_ncomb [8,6,4,2] ')
    call u%assert_equal(memsize_tree(t),             20, 'memsize_tree        ')
!
    allocate(W(memsize_tree(t)))
    do i = 1, SIZE(W)
      W(i) = i
    end do
!
    call u%assert_equal(queue_pointer(q),         [1,1,1,1], 'queue_pointer       ')
    call setup_queue(t, 1, q)
    call u%assert_equal(queue_pointer(q),       [1,9,15,19], 'queue_pointer(setup)')
!
    call u%assert_equal(current_sequence(q),      [0,0,0,0], 'current_sequence 0  ')
    call u%assert_equal(current_permutation(t, q),[1,2,3,4], 'current_permutation ')
    call set_top_node(t, q(1), 999._RK, W, .true.)
    call set_top_node(t, q(2), 999._RK, W, .true.)
    call set_top_node(t, q(3), 999._RK, W, .true.)
    call set_top_node(t, q(4), 999._RK, W, .true.)
    call u%assert_equal(current_sequence(q),      [1,1,1,1], 'current_sequence 1  ')
    call u%assert_equal(current_permutation(t, q),[1,2,3,4], 'current_permutation ')
    call u%assert_equal(current_mapping(t, q),    [0,0,0,0], 'current_mapping     ')
    call u%assert_equal(current_mapping(t, q),    [0,0,0,0], 'current_mapping     ')
!
    call set_top_node(t, q(1), 999._RK, W, .false.)
    call set_top_node(t, q(2), 999._RK, W, .false.)
    call set_top_node(t, q(3), 999._RK, W, .false.)
    call set_top_node(t, q(4), 999._RK, W, .false.)
    call u%assert_equal(current_sequence(q),      [2,2,2,2], 'current_sequence 2  ')
    call u%assert_equal(current_permutation(t, q),[1,2,3,4], 'current_permutation ')
    call u%assert_equal(current_mapping(t, q),    [1,1,1,1], 'current_mapping     ')
!
    call set_top_node(t, q(1), 999._RK, W, .false.)
    call set_top_node(t, q(2), 999._RK, W, .false.)
    call set_top_node(t, q(3), 999._RK, W, .false.)
    call set_top_node(t, q(4), 999._RK, W, .false.)
    call u%assert_equal(current_sequence(q),      [3,3,3,0], 'current_sequence 3  ')
    call u%assert_equal(current_permutation(t, q),[2,3,4,1], 'current_permutation ')
    call u%assert_equal(current_mapping(t, q),    [0,0,0,1], 'current_mapping     ')
!
  end subroutine test1
!
end program main

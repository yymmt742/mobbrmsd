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
    integer(IK)           :: i
!
    q = queue([1, 2, 3, 4])
    call setup_tree(t, q, 1, 1)
    call u%assert_equal(NINT(EXP(log_ncomb(q(:1)))),  1, 'log_ncomb [1]      ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:2)))),  3, 'log_ncomb [1,2]    ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:3)))),  9, 'log_ncomb [1,2,3]  ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:4)))), 33, 'log_ncomb [1,2,3,4]')
    call u%assert_equal(memsize_tree(t, q),          10, 'memsize_tree       ')
!
    q = queue([1, 5, 3, 1])
    call setup_tree(t, q, 1, 2)
    call u%assert_equal(NINT(EXP(log_ncomb(q(:1)))),  1, 'log_ncomb [1]      ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:2)))),  6, 'log_ncomb [1,5]    ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:3)))), 21, 'log_ncomb [1,5,3]  ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:4)))), 36, 'log_ncomb [1,5,3,1]')
    call u%assert_equal(memsize_tree(t, q),          20, 'memsize_tree       ')
!
    q = queue([2, 3, 3, 2])
    call setup_tree(t, q, 1, 3)
    call u%assert_equal(NINT(EXP(log_ncomb(q(:1)))),  2, 'log_ncomb [2]      ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:2)))),  8, 'log_ncomb [2,3]    ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:3)))), 26, 'log_ncomb [2,3,3]  ')
    call u%assert_equal(NINT(EXP(log_ncomb(q(:4)))), 62, 'log_ncomb [2,3,3,2]')
    call u%assert_equal(memsize_tree(t, q),          30, 'memsize_tree       ')
!
    allocate(W(memsize_tree(t, q)))
    do i = 1, SIZE(W)
      W(i) = COS(i**2 * 3.0_RK)
    end do
!
    call u%assert_equal(queue_pointer(q),   [1,7,16,25], 'queue_pointer      ')
!
    print*, current_sequence(q)
    call set_top_node(t, q(1), 999._RK, W, .true.)
    call set_top_node(t, q(2), 999._RK, W, .true.)
    call set_top_node(t, q(3), 999._RK, W, .true.)
    call set_top_node(t, q(4), 999._RK, W, .true.)
    print*, current_sequence(q)
    call set_top_node(t, q(1), 999._RK, W, .false.)
    call set_top_node(t, q(2), 999._RK, W, .false.)
    call set_top_node(t, q(3), 999._RK, W, .false.)
    call set_top_node(t, q(4), 999._RK, W, .false.)
    print*, current_sequence(q)
    call set_top_node(t, q(1), 999._RK, W, .false.)
    call set_top_node(t, q(2), 999._RK, W, .false.)
    call set_top_node(t, q(3), 999._RK, W, .false.)
    call set_top_node(t, q(4), 999._RK, W, .false.)
    print*, current_sequence(q)
    call set_top_node(t, q(1), 999._RK, W, .false.)
    call set_top_node(t, q(2), 999._RK, W, .false.)
    call set_top_node(t, q(3), 999._RK, W, .false.)
    call set_top_node(t, q(4), 999._RK, W, .false.)
    print*, current_sequence(q)
    call set_top_node(t, q(1), 999._RK, W, .false.)
    call set_top_node(t, q(2), 999._RK, W, .false.)
    call set_top_node(t, q(3), 999._RK, W, .false.)
    call set_top_node(t, q(4), 999._RK, W, .false.)
    print*, current_sequence(q)

!   call t%set_parent_node(W)
!   print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
!   call t%open_node()
!   call t%set_parent_node(W)
!   print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
!   call t%open_node()
!   call t%set_parent_node(W)
!   print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
!   call t%open_node()
!   call t%set_parent_node(W)
!   print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
!   W(t%upperbound) = W(t%parent_pointer())
!   call t%prune(W)
!   print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
!   call t%close_node()
!   print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
!   call t%prune(W)
!   print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
!
  end subroutine test1
!
end program main

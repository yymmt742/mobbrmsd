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
    type(tree) :: t
    real(RK), allocatable :: W(:)
    integer(IK)           :: i
!
!   t = tree(1, 1, 1, [1])
!   call u%assert_almost_equal(t%log_ncomb(),     LOG(ONE), 'log_ncomb [1]      ')
!   t = tree(1, 1, 1, [5])
!   call u%assert_almost_equal(t%log_ncomb(),  LOG(5.0_RK), 'log_ncomb [5]      ')
!   t = tree(1, 1, 3, [8,3,2])
!   call u%assert_almost_equal(t%log_ncomb(), LOG(80.0_RK), 'log_ncomb [8,3,2]  ')
!   t = tree(1, 1, 4, [1, 5, 3, 2])
!   call u%assert_almost_equal(t%log_ncomb(), LOG(51.0_RK), 'log_ncomb [1,5,3,2]')
!
!   allocate(W(t%memsize))
!   do i = 1, t%memsize
!     W(i) = COS(i**2*3D0)
!   end do
!
!   call t%open_node()
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

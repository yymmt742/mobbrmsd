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
    t = tree(1, 1, 4, [1,5,3,2])
    print*,t%memsize
    allocate(W(t%memsize))
    do i = 1, t%memsize
      W(i) = COS(i**2*3D0)
    end do
    print'(10f9.3)', W

    call t%open_node()
    call t%set_parent_node(W)
    print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
    call t%open_node()
    call t%set_parent_node(W)
    print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
    call t%open_node()
    call t%set_parent_node(W)
    print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
    call t%open_node()
    call t%set_parent_node(W)
    print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
    W(t%upperbound) = W(t%parent_pointer())
    call t%prune(W)
    print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
    call t%close_node()
    print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
    call t%prune(W)
    print *, t%nodes_pointer(), t%parent_pointer(), t%unfinished()
!
  end subroutine test1
!
end program main

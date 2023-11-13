module mod_branch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_lower_bound
  implicit none
!
  type node
!
    integer(IK), allocatable :: free_indices(:)
    real(RK)                 :: lower
!
  contains
!
    final             :: node_destroy
!
  end type node
!
contains
!
  pure function node_new(d, n, free_indices, x, y) result(res)
    integer(IK), intent(in) :: d
    integer(IK), intent(in) :: n
    integer(IK), intent(in) :: free_indices(:)
    real(RK), intent(in)    :: x(*)
    real(RK), intent(in)    :: y(*)
    real(RK)                :: w(lower_bound_worksize(d, n, free_indices))
    type(node)              :: res
    res%free_indices = free_indices
    call lower_bound(d, n, free_indices, x, y, w)
    res%lower = w(1)
  end function node_new
!
  pure elemental subroutine node_destroy(this)
    type(node), intent(inout) :: this
    if(ALLOCATED(this%free_indices)) deallocate(this%free_indices)
    this%lower = -HUGE(ZERO)
  end subroutine node_destroy
!
end module mod_branch

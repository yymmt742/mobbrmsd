module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: tree
!
  type node
    private
    sequence
    logical     :: alive
    integer(IK) :: p
  end type node
!
  type breadth
    private
    sequence
    integer(IK) :: inod, lowd, uppd
  end type breadth
!
  type tree
    integer(IK)                :: iscope = 0
    integer(IK)                :: memsize = 0
    integer(IK)                :: upperbound = 0
    type(node), allocatable    :: nodes(:)
    type(breadth), allocatable :: breadthes(:)
  contains
    procedure         :: n_depth          => tree_n_depth
    procedure         :: n_breadth        => tree_n_breadth
    procedure         :: prune            => tree_prune
!   procedure         :: lowerbounds      => tree_lowerbounds
    final             :: tree_destroy
  end type tree
!
  interface tree
    module procedure tree_new
  end interface tree
!
contains
!
  pure function tree_new(pw, memnode, ndepth, n_breadths) result(res)
    integer(IK), intent(in) :: pw, memnode, ndepth, n_breadths(ndepth)
    integer(IK)             :: n, i, j, k
    type(tree)              :: res
!
    n = SUM(n_breadths)
    res%upperbound = pw
    allocate (res%nodes(n))
    do concurrent(i=1:n)
      res%nodes(i) = node(.true., pw + 1 + (i - 1) * memnode)
    end do
!
    res%memsize = n * memnode + 1
!
    allocate (res%breadthes(ndepth))
    j = 0
    do i = 1, ndepth
      k = j + n_breadths(i)
      res%breadthes(i) = breadth(0, j, k)
      j = k
    end do
!
  end function tree_new
!
  pure elemental function tree_n_depth(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: res
    if (ALLOCATED(this%breadthes)) then
      res = SIZE(this%breadthes)
    else
      res = 0
    end if
  end function tree_n_depth
!
  pure elemental function tree_n_breadth(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: res
    res = 0
    if (ALLOCATED(this%breadthes)) then
      if (this%iscope < 1) return
      res = this%breadthes(this%iscope)%uppd - this%breadthes(this%iscope)%lowd
    end if
  end function tree_n_breadth
!
  pure subroutine tree_prune(this, W)
    class(tree), intent(inout) :: this
    real(RK), intent(in)       :: W(*)
    integer(IK)                :: i, l, u
    if (this%iscope < 1 .or. this%n_depth() < this%iscope) return
    l = this%breadthes(this%iscope)%lowd + 1
    u = this%breadthes(this%iscope)%uppd
    do concurrent(i=l:u)
      this%nodes(i)%alive = this%nodes(i)%alive .and. (W(this%nodes(i)%p) < W(this%upperbound))
    end do
  end subroutine tree_prune
!
  pure elemental subroutine tree_destroy(this)
    type(tree), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
    if (ALLOCATED(this%breadthes)) deallocate (this%breadthes)
  end subroutine tree_destroy
!
end module mod_tree

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
    private
    integer(IK)                :: iscope = 0
    integer(IK), public        :: memsize = 0
    integer(IK), public        :: upperbound = 0
    type(node), allocatable    :: nodes(:)
    type(breadth), allocatable :: breadthes(:)
  contains
    procedure         :: n_depth          => tree_n_depth
    procedure         :: n_breadth        => tree_n_breadth
    procedure         :: current_depth    => tree_current_depth
    procedure         :: reset            => tree_reset
    procedure         :: set_parent_node  => tree_set_parent_node
    procedure         :: prune            => tree_prune
    procedure         :: nodes_pointer    => tree_nodes_pointer
    procedure         :: parent_pointer   => tree_parent_pointer
    procedure         :: parent_index     => tree_parent_index
    procedure         :: current_pointer  => tree_current_pointer
    procedure         :: current_index    => tree_current_index
    procedure         :: open_node        => tree_open_node
    procedure         :: close_node       => tree_close_node
    procedure         :: finished         => tree_finished
    procedure         :: unfinished       => tree_unfinished
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
  pure elemental function tree_current_depth(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: res
    res = this%iscope
  end function tree_current_depth
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
  pure elemental function tree_nodes_pointer(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: res
    res = 0
    if (ALLOCATED(this%breadthes)) then
      if (this%iscope < 1) return
      res = this%nodes(this%breadthes(this%iscope)%lowd + 1)%p
    end if
  end function tree_nodes_pointer
!
  pure elemental function tree_parent_pointer(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: ip, res
    res = 0
    if (ALLOCATED(this%breadthes)) then
      if (this%iscope < 2) return
      ip = this%iscope - 1
      ip = this%breadthes(ip)%inod
      if (ip < 1) return
      res = this%nodes(ip)%p
    end if
  end function tree_parent_pointer
!
  pure elemental function tree_current_pointer(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: ip, res
    res = 0
    if (ALLOCATED(this%breadthes)) then
      if (this%iscope < 1) return
      ip = this%iscope
      ip = this%breadthes(ip)%inod
      if (ip < 1) return
      res = this%nodes(ip)%p
    end if
  end function tree_current_pointer
!
  pure elemental function tree_current_index(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: ip, res
    res = 0
    if (ALLOCATED(this%breadthes)) then
      if (this%iscope < 1) return
      ip = this%iscope
      ip = this%breadthes(ip)%inod - &
     &     this%breadthes(ip)%lowd
      if (ip > 0) res = ip
    end if
  end function tree_current_index
!
  pure elemental function tree_parent_index(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: ip, res
    res = 0
    if (ALLOCATED(this%breadthes)) then
      if (this%iscope < 2) return
      ip = this%iscope - 1
      ip = this%breadthes(ip)%inod - &
     &     this%breadthes(ip)%lowd
      if (ip > 0) res = ip
    end if
  end function tree_parent_index
!
  pure subroutine tree_reset(this)
    class(tree), intent(inout) :: this
    this%iscope = 0
    call this%open_node()
  end subroutine tree_reset
!
  pure subroutine tree_set_parent_node(this, W)
    class(tree), intent(inout) :: this
    real(RK), intent(in)       :: W(*)
    real(RK)                   :: lv
    integer(IK)                :: i, l, u
!
    if (this%iscope < 1 .or. this%n_depth() < this%iscope) return
    l = this%breadthes(this%iscope)%lowd + 1
    u = this%breadthes(this%iscope)%uppd
!
    this%breadthes(this%iscope)%inod = l
!
    lv = RHUGE
    do i = l, u
      if (this%nodes(i)%alive .and. W(this%nodes(i)%p) < lv) then
        this%breadthes(this%iscope)%inod = i
        lv = W(this%nodes(i)%p)
      end if
    end do
!
    if(this%iscope==this%n_depth()) this%nodes(l:u)%alive = .false.
!
  end subroutine tree_set_parent_node
!
  pure elemental subroutine tree_open_node(this)
    class(tree), intent(inout) :: this
    integer(IK)                :: i, l, u
!
    if (this%iscope < 0 .or. this%n_depth() <= this%iscope) return
    this%iscope = this%iscope + 1
    l = this%breadthes(this%iscope)%lowd + 1
    u = this%breadthes(this%iscope)%uppd
!
    do concurrent(i=l:u)
      this%nodes(i)%alive = .true.
    end do
!
    this%breadthes(this%iscope)%inod = l
!
  end subroutine tree_open_node
!
  pure elemental subroutine tree_close_node(this)
    class(tree), intent(inout) :: this
    integer(IK)                :: p
!
    if (this%iscope <= 1 .or. this%n_depth() < this%iscope) return
    this%iscope = this%iscope - 1
    p = this%breadthes(this%iscope)%inod
    this%nodes(p)%alive = .false.
!
  end subroutine tree_close_node
!
  pure elemental function tree_finished(this) result(res)
    class(tree), intent(in) :: this
    logical                 :: res
!
    res = .not. this%unfinished()
!
  end function tree_finished
!
  pure elemental function tree_unfinished(this) result(res)
    class(tree), intent(in) :: this
    logical                 :: res
    integer(IK)             :: i, l, u
!
    res = .true.
    if (this%iscope < 1 .or. this%n_depth() < this%iscope) return
    l = this%breadthes(this%iscope)%lowd + 1
    u = this%breadthes(this%iscope)%uppd
    res = ANY([(this%nodes(i)%alive, i=l, u)])
!
  end function tree_unfinished
!
  pure subroutine tree_prune(this, W)
    class(tree), intent(inout) :: this
    real(RK), intent(in)       :: W(*)
    integer(IK)                :: i, l, u
!
    if (this%iscope < 1 .or. this%n_depth() < this%iscope) return
    l = this%breadthes(this%iscope)%lowd + 1
    u = this%breadthes(this%iscope)%uppd
    do concurrent(i=l:u)
      this%nodes(i)%alive = this%nodes(i)%alive .and. (W(this%nodes(i)%p) < W(this%upperbound))
    end do
!
  end subroutine tree_prune
!
  pure elemental subroutine tree_destroy(this)
    type(tree), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
    if (ALLOCATED(this%breadthes)) deallocate (this%breadthes)
  end subroutine tree_destroy
!
end module mod_tree

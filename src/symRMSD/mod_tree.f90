module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: node
  public :: breadth
  public :: tree
  public :: memsize_tree
  public :: setup_tree
  public :: n_breadth
!
!| Node.
  type node
    private
    sequence
!|  p : pointer to memory. if p<0, the node has been explored.
    integer(IK) :: p
  end type node
!
!| A collection of nodes in a hierarchy.
  type breadth
    private
    sequence
!|  inod : current pointer
    integer(IK) :: inod
!|  lowd : (Pointer to the first node) - 1.
    integer(IK) :: lowd
!|  uppd : Pointer to the last node.
    integer(IK) :: uppd
  end type breadth
!
!| Tree.
  type tree
    private
    sequence
!|  number of hierarchy.
    integer(IK)       :: ndepth
!|  number of nodes.
    integer(IK)       :: nnodes
!|  current hierarchy.
    integer(IK)       :: idepth
!|  Memory size per node.
    integer(IK)       :: memnode
! contains
!   procedure         :: n_depth          => tree_n_depth
!   procedure         :: n_breadth        => tree_n_breadth
!   procedure         :: current_depth    => tree_current_depth
!   procedure         :: reset            => tree_reset
!   procedure         :: set_parent_node  => tree_set_parent_node
!   procedure         :: set_lowerbound   => tree_set_lowerbound
!   procedure         :: prune            => tree_prune
!   procedure         :: nodes_pointer    => tree_nodes_pointer
!   procedure         :: parent_pointer   => tree_parent_pointer
!   procedure         :: parent_index     => tree_parent_index
!   procedure         :: current_pointer  => tree_current_pointer
!   procedure         :: current_index    => tree_current_index
!   procedure         :: alive_nodes      => tree_alive_nodes
!   procedure         :: open_node        => tree_open_node
!   procedure         :: close_node       => tree_close_node
!   procedure         :: log_ncomb        => tree_log_ncomb
!   procedure         :: finished         => tree_finished
!   procedure         :: unfinished       => tree_unfinished
!   procedure         :: clear            => tree_clear
!   final             :: tree_destroy
  end type tree
!
  interface node
    module procedure node_new
  end interface node
!
  interface breadth
    module procedure breadth_new
  end interface breadth
!
  interface tree
    module procedure tree_new
  end interface tree
!
contains
!
  pure elemental function node_new() result(res)
    type(node) :: res
    res%p = 0
  end function node_new
!
  pure elemental function breadth_new(n_nodes) result(res)
!| n_nodes :: number of nodes in breadth, n_nodes>0.
    integer(IK), intent(in) :: n_nodes
    type(breadth)           :: res
    res%inod = 0
    res%lowd = 0
    res%uppd = MAX(1, n_nodes)
  end function breadth_new
!
  pure function tree_new(b, memnode) result(res)
!|  b :: breadth list, must be intiialized.
    type(breadth), intent(in) :: b(:)
!|  memnode :: memory size of each node.
    integer(IK), intent(in)   :: memnode
    type(tree)                :: res
!
    res%idepth  = 1
    res%nnodes  = SUM(b%uppd - b%lowd)
    res%ndepth  = SIZE(b)
    res%memnode = MAX(memnode, 1)
!
  end function tree_new
!
  pure elemental function memsize_tree(this) result(res)
!| n_nodes :: number of nodes in breadth, n_nodes>0.
    type(tree), intent(in) :: this
    integer(IK)            :: res
    res = this%nnodes * this%memnode
  end function memsize_tree
!
  pure subroutine setup_tree(this, b, n, p)
    type(tree), intent(in)       :: this
!|  b :: breadthes
    type(breadth), intent(inout) :: b(*)
!|  n :: nodes
    type(node), intent(inout)    :: n(*)
!|  p :: pointer
    integer(IK), intent(in)      :: p
    integer(IK)                  :: i, j, k
!
    do concurrent(i=1:this%nnodes)
      n(i)%p = p + (i - 1) * this%memnode
    end do
!
    j = 0
    do i = 1, this%nnodes
      k = j + b(i)%uppd - b(i)%lowd
      b(i)%lowd = j
      b(i)%uppd = k
      j = k
    end do
!
  end subroutine setup_tree
!
  pure function n_breadth(t, b) result(res)
!|  t :: tree
    type(tree), intent(in)    :: t
!|  b :: breadth list
    type(breadth), intent(in) :: b(*)
    integer(IK)               :: res
    res = 0
    if (t%idepth < 1) return
    res = b(t%idepth)%uppd - b(t%idepth)%lowd
  end function n_breadth
!
! pure elemental function tree_nodes_pointer(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: res
!   res = 0
!   if (ALLOCATED(this%breadthes)) then
!     if (this%iscope < 1) return
!     res = this%nodes(this%breadthes(this%iscope)%lowd + 1)%p
!   end if
! end function tree_nodes_pointer
!
! pure elemental function tree_parent_pointer(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: ip, res
!   res = 0
!   if (ALLOCATED(this%breadthes)) then
!     if (this%iscope < 2) return
!     ip = this%iscope - 1
!     ip = this%breadthes(ip)%inod
!     if (ip < 1) return
!     res = this%nodes(ip)%p
!   end if
! end function tree_parent_pointer
!
! pure elemental function tree_current_pointer(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: ip, res
!   res = 0
!   if (ALLOCATED(this%breadthes)) then
!     if (this%iscope < 1) return
!     ip = this%iscope
!     ip = this%breadthes(ip)%inod
!     if (ip < 1) return
!     res = this%nodes(ip)%p
!   end if
! end function tree_current_pointer
!
! pure elemental function tree_current_index(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: ip, res
!   res = 0
!   if (ALLOCATED(this%breadthes)) then
!     if (this%iscope < 1) return
!     ip = this%iscope
!     ip = this%breadthes(ip)%inod - &
!    &     this%breadthes(ip)%lowd
!     if (ip > 0) res = ip
!   end if
! end function tree_current_index
!
! pure function tree_alive_nodes(this) result(res)
!   class(tree), intent(in) :: this
!   logical, allocatable    :: res(:)
!   integer(IK)             :: l, u
!   allocate(res(0))
!   if (this%iscope < 1 .or. this%n_depth() < this%iscope) return
!   l = this%breadthes(this%iscope)%lowd + 1
!   u = this%breadthes(this%iscope)%uppd
!   res = this%nodes(l:u)%alive
! end function tree_alive_nodes
!
! pure elemental function tree_parent_index(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: ip, res
!   res = 0
!   if (ALLOCATED(this%breadthes)) then
!     if (this%iscope < 2) return
!     ip = this%iscope - 1
!     ip = this%breadthes(ip)%inod - &
!    &     this%breadthes(ip)%lowd
!     if (ip > 0) res = ip
!   end if
! end function tree_parent_index
!
! pure subroutine tree_reset(this)
!   class(tree), intent(inout) :: this
!   this%iscope = 0
!   if (ALLOCATED(this%nodes)) this%nodes%alive = .false.
!   call this%open_node()
! end subroutine tree_reset
!
! pure subroutine tree_set_parent_node(this, W)
!   class(tree), intent(inout) :: this
!   real(RK), intent(in)       :: W(*)
!   real(RK)                   :: lv
!   integer(IK)                :: i, l, u
!
!   if (this%iscope < 1 .or. this%n_depth() < this%iscope) return
!   l = this%breadthes(this%iscope)%lowd + 1
!   u = this%breadthes(this%iscope)%uppd
!
!   this%breadthes(this%iscope)%inod = l
!
!   lv = RHUGE
!   do i = l, u
!     if (this%nodes(i)%alive .and. W(this%nodes(i)%p) < lv) then
!       this%breadthes(this%iscope)%inod = i
!       lv = W(this%nodes(i)%p)
!     end if
!   end do
!
!   if(this%iscope==this%n_depth()) this%nodes(l:u)%alive = .false.
!
! end subroutine tree_set_parent_node
!
! pure elemental subroutine tree_open_node(this)
!   class(tree), intent(inout) :: this
!   integer(IK)                :: i, l, u
!
!   if (this%iscope < 0 .or. this%n_depth() <= this%iscope) return
!   if (this%iscope > 0) then
!     l = this%breadthes(this%iscope)%inod
!     if (l > 0) this%nodes(l)%alive = .false.
!   end if
!   this%iscope = this%iscope + 1
!   l = this%breadthes(this%iscope)%lowd + 1
!   u = this%breadthes(this%iscope)%uppd
!
!   do concurrent(i=l:u)
!     this%nodes(i)%alive = .true.
!   end do
!
!   this%breadthes(this%iscope)%inod = l
!
! end subroutine tree_open_node
!
! pure subroutine tree_set_lowerbound(this, W)
!   class(tree), intent(in) :: this
!   real(RK), intent(inout) :: W(*)
!   integer(IK)             :: i, nnod
!
!   w(this%lowerbound) = w(this%upperbound)
!   if (.not. ALLOCATED(this%nodes)) return
!   nnod = SIZE(this%nodes)
!   do i = 1, nnod
!     if (this%nodes(i)%alive) w(this%lowerbound) = MIN(w(this%lowerbound), w(this%nodes(i)%p))
!   end do
!
! end subroutine tree_set_lowerbound
!
! pure elemental subroutine tree_close_node(this)
!   class(tree), intent(inout) :: this
!   integer(IK)                :: p
!
!   if (this%iscope <= 1 .or. this%n_depth() < this%iscope) return
!   this%iscope = this%iscope - 1
!   p = this%breadthes(this%iscope)%inod
!
! end subroutine tree_close_node
!
! pure elemental function tree_log_ncomb(this) result(res)
!   class(tree), intent(in) :: this
!   real(RK)                :: tmp, res
!   integer(IK)             :: i, n
!   res = ZERO
!   if (.not. ALLOCATED(this%breadthes)) return
!   n = SIZE(this%breadthes)
!   if (n < 1) return
!   tmp = ZERO
!   do i = n, 2, -1
!     tmp = tmp - log_nnod(this%breadthes(i))
!     res = res + EXP(tmp)
!   end do
!   res = log_nnod(this%breadthes(1)) - tmp + LOG(ONE + res)
! end function tree_log_ncomb
!
! pure elemental function log_nnod(b) result(res)
!   type(breadth), intent(in) :: b
!   real(RK)                  :: res
!   res = LOG(REAL(b%uppd - b%lowd, RK))
! end function log_nnod
!
! pure elemental function tree_finished(this) result(res)
!   class(tree), intent(in) :: this
!   logical                 :: res
!
!   res = .not. this%unfinished()
!
! end function tree_finished
!
! pure elemental function tree_unfinished(this) result(res)
!   class(tree), intent(in) :: this
!   logical                 :: res
!   integer(IK)             :: i, l, u
!
!   res = .true.
!   if (this%iscope < 1 .or. this%n_depth() < this%iscope) return
!   l = this%breadthes(this%iscope)%lowd + 1
!   u = this%breadthes(this%iscope)%uppd
!   res = ANY([(this%nodes(i)%alive, i=l, u)])
!
! end function tree_unfinished
!
! pure subroutine tree_prune(this, W)
!   class(tree), intent(inout) :: this
!   real(RK), intent(in)       :: W(*)
!   integer(IK)                :: i, l, u
!
!   if (this%iscope < 1 .or. this%n_depth() < this%iscope) return
!   l = this%breadthes(this%iscope)%lowd + 1
!   u = this%breadthes(this%iscope)%uppd
!   do concurrent(i=l:u)
!     this%nodes(i)%alive = this%nodes(i)%alive .and. (W(this%nodes(i)%p) < W(this%upperbound))
!   end do
!
! end subroutine tree_prune
!
! pure elemental subroutine tree_clear(this)
!   class(tree), intent(inout) :: this
!   if (ALLOCATED(this%nodes)) deallocate (this%nodes)
!   if (ALLOCATED(this%breadthes)) deallocate (this%breadthes)
! end subroutine tree_clear
!
! pure elemental subroutine tree_destroy(this)
!   type(tree), intent(inout) :: this
!   call tree_clear(this)
! end subroutine tree_destroy
!
end module mod_tree


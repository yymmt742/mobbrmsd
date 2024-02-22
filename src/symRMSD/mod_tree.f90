module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: queue
  public :: nnode
  public :: nnodes
  public :: tree
  public :: setup_tree
  public :: memsize_tree
  public :: set_top_node
  public :: current_sequence
  public :: queue_pointer
  public :: node_pointer
  public :: is_empty
  public :: log_ncomb
!
!| A collection of nodes in a hierarchy.
  type queue
    private
    sequence
!|  i : Current pointer
    integer(IK) :: i
!|  p : pointer to memory.
    integer(IK) :: p
!|  n : nnodes.
    integer(IK) :: n
  end type queue
!
!| Tree.
  type tree
    private
    sequence
!|  number of hierarchy.
    integer(IK)       :: ndepth
!|  number of nodes.
    integer(IK)       :: nnodes
!|  Memory size per node.
    integer(IK)       :: nmnode
  end type tree
!
  interface queue
    module procedure queue_new
  end interface queue
!
contains
!
!| Constructer of queue
  pure elemental function queue_new(n_nodes) result(res)
!|  n_nodes :: number of nodes in queue, n_nodes>0.
    integer(IK), intent(in) :: n_nodes
    type(queue)             :: res
    res%i = -1
    res%p = 1
    res%n = MAX(1, n_nodes)
  end function queue_new
!
!| Count number of nodes in a queue.
  pure elemental function nnode(q) result(res)
!|  q :: queue
    type(queue), intent(in) :: q
    integer(IK)             :: res
    res = q%n
  end function nnode
!
!| Count number of nodes
  pure function nnodes(q) result(res)
!|  q :: queue list
    type(queue), intent(in) :: q(:)
    integer(IK)             :: res
    res = SUM(nnode(q))
  end function nnodes
!
!| Set up the tree.
  pure subroutine setup_tree(t, q, p, nmnode)
!|  t :: tree
    type(tree), intent(inout)  :: t
!|  q :: queue list, must be intiialized.
    type(queue), intent(inout) :: q(:)
!|  p :: pointer
    integer(IK), intent(in)    :: p
!|  nmnode :: memory size of each node.
    integer(IK), intent(in)    :: nmnode
    integer(IK)                :: i
!
    t%ndepth = SIZE(q)
    t%nnodes = nnodes(q)
    t%nmnode = MAX(nmnode, 1)
!
    if (t%nnodes < 1) return
    q(1)%p = p
    do i = 2, t%ndepth
      q(i)%p = t%nmnode * q(i - 1)%n + q(i - 1)%p
    end do
!
  end subroutine setup_tree
!
!| Inquire memsize of tree.
  pure function memsize_tree(t) result(res)
!|  t :: tree
    type(tree), intent(in)  :: t
    integer(IK)             :: res
    res = t%nnodes * t%nmnode
  end function memsize_tree
!
!| Returns a pointer to the current queue.
  pure elemental function queue_pointer(q) result(res)
!|  q :: queue
    type(queue), intent(in) :: q
    integer(IK)             :: res
    res = q%p
  end function queue_pointer
!
!| Returns a pointer to the current best node.
  pure elemental function node_pointer(t, q) result(res)
!|  t :: tree
    type(tree), intent(in)  :: t
!|  q :: queue
    type(queue), intent(in) :: q
    integer(IK)             :: res
    res = q%p + q%i * t%nmnode
  end function node_pointer
!
!| Returns a pointer to the current queue.
  pure function current_sequence(q) result(res)
!|  q :: queue
    type(queue), intent(in) :: q(:)
    integer(IK)             :: i, res(SIZE(q))
    do concurrent(i=1:SIZE(q))
      res(i) = q(i)%i + 1
    end do
  end function current_sequence
!
!| Set top node
  pure subroutine set_top_node(t, q, UB, W, reset)
!|  t :: tree
    type(tree), intent(in)     :: t
!|  q :: queue
    type(queue), intent(inout) :: q
!|  UB :: upperbound
    real(RK), intent(in)       :: UB
!|  W :: work array
    real(RK), intent(in)       :: W(*)
!|  reset :: If true, treat all nodes as unexplored.
    logical, intent(in)        :: reset
    real(RK)                   :: uv, lv
    integer(IK)                :: i, p
!
    if (q%i < 0 .or. reset) then
      lv = - RHUGE
    else
      lv = W(q%p + q%i * t%nmnode)
    end if
!
    q%i = -1
    uv = UB
!
    if (uv < lv) return
!
    do i = 0, q%n-1
      p = q%p + i * t%nmnode
      if (lv < W(p) .and. W(p) < uv) then
        q%i = i
        uv = W(p)
      end if
    end do
!
  end subroutine set_top_node
!
  pure function log_ncomb(q) result(res)
    type(queue), intent(in) :: q(:)
    real(RK)                :: tmp, res
    integer(IK)             :: i
!
    res = ZERO
    tmp = ZERO
    do i = SIZE(q), 1, -1
      tmp = tmp - LOG(REAL(q(i)%n, RK))
      res = res + EXP(tmp)
    end do
    res = ONE + res - EXP(tmp)
!
    if (res < 1.E-24_RK) then; res = -RHUGE
    else; res = LOG(res) - tmp
    end if
!
  end function log_ncomb
!
  pure elemental function is_empty(q) result(res)
    type(queue), intent(in) :: q
    logical                 :: res
    res = q%i < 0
  end function is_empty
!
! pure elemental function tree_parent_pointer(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: ip, res
!   res = 0
!   if (ALLOCATED(this%queuees)) then
!     if (this%iscope < 2) return
!     ip = this%iscope - 1
!     ip = this%queuees(ip)%inod
!     if (ip < 1) return
!     res = this%nodes(ip)%p
!   end if
! end function tree_parent_pointer
!
! pure elemental function tree_current_pointer(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: ip, res
!   res = 0
!   if (ALLOCATED(this%queuees)) then
!     if (this%iscope < 1) return
!     ip = this%iscope
!     ip = this%queuees(ip)%inod
!     if (ip < 1) return
!     res = this%nodes(ip)%p
!   end if
! end function tree_current_pointer
!
! pure elemental function tree_current_index(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: ip, res
!   res = 0
!   if (ALLOCATED(this%queuees)) then
!     if (this%iscope < 1) return
!     ip = this%iscope
!     ip = this%queuees(ip)%inod - &
!    &     this%queuees(ip)%lowd
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
!   l = this%queuees(this%iscope)%lowd + 1
!   u = this%queuees(this%iscope)%uppd
!   res = this%nodes(l:u)%alive
! end function tree_alive_nodes
!
! pure elemental function tree_parent_index(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: ip, res
!   res = 0
!   if (ALLOCATED(this%queuees)) then
!     if (this%iscope < 2) return
!     ip = this%iscope - 1
!     ip = this%queuees(ip)%inod - &
!    &     this%queuees(ip)%lowd
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
! pure subroutine climb(t, q, n)
!   type(tree), intent(inout)  :: t
!   type(queue), intent(inout) :: q(:)
!   type(node), intent(inout)  :: n(:)
!   integer(IK)                :: i
!
!   if (t%ndepth <= t%idepth) return
!
!   t%idepth = t%idepth + 1
!   q(t%idepth)%i = q(t%idepth)%l
!
!   do concurrent(i=q(t%idepth)%l:q(t%idepth)%u)
!     n(i)%p = ABS(n(i)%p)
!   end do
!
! end subroutine climb
!
! pure elemental subroutine go_down(t)
!   type(tree), intent(inout) :: t
!
!   if (t%idepth > 0) t%idepth = t%idepth - 1
!
! end subroutine go_down
!
end module mod_tree


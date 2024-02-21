module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: node
  public :: breadth
  public :: tree
  public :: memsize_tree
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
    sequence
    integer(IK)                :: p
    integer(IK)                :: iscope
    integer(IK), public        :: ubnode
    integer(IK), public        :: memnode
    integer(IK), public        :: memsize
    integer(IK), public        :: upperbound
    integer(IK), public        :: lowerbound
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
    res%alive = .TRUE.
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
  pure function tree_new(memnode, n_breadths) result(res)
!|  memnode :: memory size of each node.
    integer(IK), intent(in) :: memnode
!|  n_breadths :: number of nodes for each level.
    integer(IK), intent(in) :: n_breadths(:)
    type(tree)              :: res
    integer(IK)             :: pq, n
!
    n = SUM(n_breadths)
    res%p       = 1
    res%memnode = MAX(memnode, 1)
    res%memsize = n * res%memnode + 1
!
    pq = 1
    res%lowerbound = pq; pq = pq + 1
    res%ubnode     = pq; pq = pq + res%memnode
    res%upperbound = res%ubnode
!
!   allocate (res%nodes(n))
!   do concurrent(i=1:n)
!     res%nodes(i) = node(.true., pq + (i - 1) * res%memnode)
!   end do
!
!   allocate (res%breadthes(ndepth))
!   j = 0
!   do i = 1, ndepth
!     k = j + n_breadths(i)
!     res%breadthes(i) = breadth(0, j, k)
!     j = k
!   end do
!
  end function tree_new
!
  pure elemental function memsize_tree(this) result(res)
!| n_nodes :: number of nodes in breadth, n_nodes>0.
    type(tree), intent(in) :: this
    integer(IK)            :: res
    res = this%memsize
  end function memsize_tree
!
! pure subroutine tree_setup(this, b, n)
!   class(tree), intent(in)        :: this
!| b :: breadthes
!   type(breadthes), intent(inout) :: b(:)
!| n :: nodes
!   type(breadthes), intent(inout) :: n(:)
!   integer(IK)                    :: i, j, k
!   j = 0
!   do i = 1, SIZE(b)
!     k = j + b(i)%uppd - b(i)%lowd
!     b(i)%lowd = j
!     b(i)%uppd = k
!     j = k
!   end do
! end subroutine breadthes_setup
!
! pure subroutine nodes_setup(n, b)
!   type(nodes), intent(inout)     :: n(:)
!   type(breadthes), intent(inout) :: b(:)
!   integer(IK)                    :: i, j, k
!   j = 0
!   do i = 1, SIZE(b)
!     k = j + b(i)%uppd - b(i)%lowd
!     b(i)%lowd = j
!     b(i)%uppd = k
!     j = k
!   end do
! end subroutine nodes_setup
!
! pure elemental function tree_n_depth(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: res
!   if (ALLOCATED(this%breadthes)) then
!     res = SIZE(this%breadthes)
!   else
!     res = 0
!   end if
! end function tree_n_depth
!
! pure elemental function tree_current_depth(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: res
!   res = this%iscope
! end function tree_current_depth
!
! pure elemental function tree_n_breadth(this) result(res)
!   class(tree), intent(in) :: this
!   integer(IK)             :: res
!   res = 0
!   if (ALLOCATED(this%breadthes)) then
!     if (this%iscope < 1) return
!     res = this%breadthes(this%iscope)%uppd - this%breadthes(this%iscope)%lowd
!   end if
! end function tree_n_breadth
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

!| Tree structure for branch and bound.
module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, TEN => RTEN, ZERO => RZERO, RHUGE, LN_TO_L10
  use mod_mol_block
  implicit none
  private
  public :: tree
!
!|  queue, a collection of nodes in a hierarchy.
!   data is stored in heap by turning S first, [[a_1,...,a_S], [b_1,...,b_S],...]
  type queue
    private
    sequence
!|  i :: Current pointer.<br>
!        i in [0,1,...,n-1] :: current node<br>
!        i == -1 :: unexplored<br>
!        i <  -1 :: explored
    integer(IK) :: i
!|  p :: pointer to memory.
    integer(IK) :: p
!|  n :: nnodes.
    integer(IK) :: n
!|  x :: Memory size per node.
    integer(IK) :: x
  end type queue
!
!| Factorial tree.
  type tree
    private
!|  p :: pointer to work array.
    integer(IK), public :: p
!|  l :: current level.
    integer(IK)         :: l
!|  s :: scaling factor.
    integer(IK)         :: s
!|  c :: current queue.
    type(queue)         :: c
!|  q :: queue, a collection of nodes in a hierarchy.
    type(queue), allocatable :: q(:)
  contains
    procedure :: memsize             => tree_memsize
    procedure :: current_level       => tree_current_level
    procedure :: queue_pointer       => tree_queue_pointer
    procedure :: current_pointer     => tree_current_pointer
    procedure :: node_pointer        => tree_node_pointer
    procedure :: current_sequence    => tree_current_sequence
    procedure :: current_permutation => tree_current_permutation
    procedure :: current_mapping     => tree_current_mapping
    procedure :: expand              => tree_expand
    procedure :: leave               => tree_leave
    procedure :: select_top_node     => tree_select_top_node
    procedure :: n_perm              => tree_n_perm
    procedure :: n_depth             => tree_n_depth
    procedure :: log_ncomb           => tree_log_ncomb
    procedure :: ncomb_frac          => tree_ncomb_frac
    procedure :: ncomb_exp           => tree_ncomb_exp
    procedure :: queue_is_empty      => tree_queue_is_empty
    procedure :: queue_is_bottom     => tree_queue_is_bottom
    final     :: tree_destroy
  end type tree
!
  interface tree
    module procedure tree_new
  end interface tree
!
contains
!
!| Constructer of queue
  pure elemental function queue_new(n_nodes, memsize, p) result(res)
!|  n_nodes :: number of nodes in queue, n_nodes>0.
    integer(IK), intent(in) :: n_nodes
!|  memsize :: memory size of each node.
    integer(IK), intent(in) :: memsize
!|  p :: offset of pointer
    integer(IK), intent(in) :: p
    type(queue)             :: res
    res%i = -1
    res%p = p
    res%n = MAX(1, n_nodes)
    res%x = MAX(1, memsize)
  end function queue_new
!
!| Constructer of factorial tree.<br>
!  [s*m, s*(m-1),..., s*2, s]
  pure function tree_new(b, memsize) result(res)
!|  b :: mol_block.
    type(mol_block), intent(in) :: b
    interface
      pure function memsize(b, p) result(res)
        use mod_params, only: IK
        use mod_mol_block, only: mol_block
        type(mol_block), intent(in) :: b
        integer(IK), intent(in)     :: p
        integer(IK)                 :: res
      end function memsize
    end interface
    type(tree)                  :: res
    integer(IK)                 :: i
    res%p = 1
    res%s = MAX(b%s, 1)
    allocate (res%q(b%n2 + 1))
!
    res%q(1) = queue_new(1, MAX(memsize(b, 0), 1), 0) ! root node
    do i = 1, b%n2
      res%q(i + 1) = queue_new(res%s * (b%n1 - i + 1), MAX(memsize(b, i), 1),&
     &                         res%q(i)%x * res%q(i)%n + res%q(i)%p)
    end do
!
    res%l = 1
!
!   designate root as the current node.
    res%q(1)%i = 0
    res%c = res%q(1)
!
  end function tree_new
!
!| Inquire memsize of tree.
  pure elemental function memsize_queue(q) result(res)
!|  t :: tree
    type(queue), intent(in)  :: q
    integer(IK)              :: res
    res = q%x * q%n
  end function memsize_queue
!
!| Inquire memsize of tree.
  pure elemental function tree_memsize(this) result(res)
!|  t :: tree
    class(tree), intent(in)  :: this
    integer(IK)              :: res
    res = SUM(memsize_queue(this%q))
  end function tree_memsize
!
!| Returns a pointer to the queue.
  pure elemental function tree_current_level(this) result(res)
!|  this :: tree
    class(tree), intent(in)  :: this
    integer(IK)              :: res
    res = this%l - 1
  end function tree_current_level
!
!| Returns a pointer to the queue.
  pure elemental function tree_queue_pointer(this) result(res)
!|  this :: tree
    class(tree), intent(in)  :: this
    integer(IK)              :: res
    res = this%p + this%c%p
  end function tree_queue_pointer
!
!| Returns a pointer to the current queue.
  pure elemental function tree_current_pointer(this) result(res)
!|  this :: tree
    class(tree), intent(in)  :: this
    integer(IK)              :: res
    res = this%p + this%c%p + this%c%i * this%c%x
  end function tree_current_pointer
!
!| Returns a pointer to the current best node.
  pure elemental function tree_node_pointer(this, iper, imap) result(res)
!|  this :: tree
    class(tree), intent(in) :: this
!|  iper :: permutation index, must be [0,1,...,q%n/s-1].
    integer(IK), intent(in) :: iper
!|  imap :: mapping index, must be [0,1,...,s-1].
    integer(IK), intent(in) :: imap
    integer(IK)             :: res
    res = this%p + this%c%p + (iper * this%s + imap) * this%c%x
  end function tree_node_pointer
!
!| Returns current sequence.
  pure function tree_current_sequence(this) result(res)
!|  t :: tree
    class(tree), intent(in) :: this
    integer(IK)             :: i, res(SIZE(this%q)-1)
    do concurrent(i=2:SIZE(this%q))
      res(i - 1) = this%q(i)%i + 1
    end do
  end function tree_current_sequence
!
!| Returns current permutation.
  pure function tree_current_permutation(this) result(res)
!|  t :: tree
    class(tree), intent(in) :: this
    integer(IK)             :: i, p, res(SIZE(this%q) - 1)
    do concurrent(i=1:SIZE(res))
      res(i) = i
    end do
    do i = 1, SIZE(res)
      if (this%q(i + 1)%i < 0) return
      if(this%q(i + 1)%i<this%s) cycle
      p = i + this%q(i + 1)%i / this%s
      res(i:p) = [res(p), res(i:p - 1)]
    end do
  end function tree_current_permutation
!
!| Returns current mapping.
  pure function tree_current_mapping(this) result(res)
!|  this :: tree
    class(tree), intent(in) :: this
    integer(IK)             :: i, res(SIZE(this%q)-1)
    do concurrent(i=2:SIZE(this%q))
      res(i - 1) = MODULO(this%q(i)%i, this%s)
    end do
  end function tree_current_mapping
!
!| Expand current node
  pure elemental subroutine tree_expand(this)
!|  this :: tree
    class(tree), intent(inout) :: this
    if(this%queue_is_empty().or.this%queue_is_bottom()) return
    this%l = this%l + 1
    this%c = this%q(this%l)
    this%q(this%l)%i = 1
    this%c%i = -1
  end subroutine tree_expand
!
!| Leave current node
  pure elemental subroutine tree_leave(this)
!|  this :: tree
    class(tree), intent(inout) :: this
    if(this%l<2) return
    this%l = this%l - 1
    this%c = this%q(this%l)
  end subroutine tree_leave
!
!| Select top node
  pure subroutine tree_select_top_node(this, UB, W)
!|  this :: tree
    class(tree), intent(inout) :: this
!|  UB :: upperbound
    real(RK), intent(in)       :: UB
!|  W :: work array
    real(RK), intent(in)       :: W(*)
    real(RK)                   :: uv, lv
    integer(IK)                :: i, p
!
    if (this%c%i == -1) then
      lv = -RHUGE
    elseif (this%c%i < -1) then
      return
    else
      lv = W(this%current_pointer())
    end if
!
    this%c%i = -2
    uv = UB
!
    if (uv < lv) return
!
    p = this%queue_pointer()
    do i = 0, this%c%n - 1
      if (lv < W(p) .and. W(p) < uv) then
        this%c%i = i
        uv = W(p)
      end if
      p = p + this%c%x
    end do
!
    if (this%c%i < 0) return
    this%q(this%l)%i = this%c%i
!
  end subroutine tree_select_top_node
!
!|  returns number of permutation in q.
  pure elemental function tree_n_perm(this) result(res)
!|  this :: tree.
    class(tree), intent(in) :: this
    integer(IK)             :: res
    res = this%c%n / this%s
  end function tree_n_perm
!
  pure elemental function tree_n_depth(this) result(res)
!|  this :: tree.
    class(tree), intent(in) :: this
    integer(IK)             :: res
    res = SIZE(this%q) - 1
  end function tree_n_depth
!
!|  returns number of nodes in tree.
  pure elemental function tree_log_ncomb(this) result(res)
!|  this :: tree.
    class(tree), intent(in) :: this
    real(RK)                :: tmp, res
    integer(IK)             :: i
!
    res = ZERO
    tmp = ZERO
    do i = SIZE(this%q), 1, -1
      tmp = tmp - LOG(REAL(this%q(i)%n, RK))
      res = res + EXP(tmp)
    end do
    res = ONE + res - EXP(tmp)
!
    if (res < 1.E-24_RK) then; res = -RHUGE
    else; res = LOG(res) - tmp
    end if
!
  end function tree_log_ncomb
!
!| returns number of nodes in fraction.
  pure function tree_ncomb_frac(this) result(res)
!|  this :: tree.
    class(tree), intent(in) :: this
    real(RK)                :: tmp, res
    tmp = LN_TO_L10 * this%log_ncomb()
    res = TEN**(tmp - REAL(INT(tmp), RK))
  end function tree_ncomb_frac
!
!| returns number of nodes in exp.
  pure function tree_ncomb_exp(this) result(res)
!|  this :: tree.
    class(tree), intent(in) :: this
    integer(IK)             :: res
!
    res = INT(LN_TO_L10 * this%log_ncomb(), IK)
!
  end function tree_ncomb_exp
!
  pure elemental function tree_queue_is_empty(this) result(res)
!|  this :: tree.
    class(tree), intent(in) :: this
    logical                 :: res
    res = this%c%i < 0
  end function tree_queue_is_empty
!
  pure elemental function tree_queue_is_bottom(this) result(res)
!|  this :: tree.
    class(tree), intent(in) :: this
    logical                 :: res
    res = this%l == SIZE(this%q)
  end function tree_queue_is_bottom
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
  pure elemental subroutine tree_destroy(this)
    type(tree), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine tree_destroy
!
end module mod_tree


!| Tree structure for branch and bound. <br>
!   queue, a collection of nodes in a hierarchy. <br>
!   data is stored in heap by turning S first, [[a_1,...,a_S], [b_1,...,b_S],...] <br>
!   i :: Current pointer.<br>
!        i in [0,1,...,n-1] :: current node<br>
!        i == -1 :: unexplored<br>
!        i <  -1 :: explored
module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, TEN => RTEN, ZERO => RZERO, RHUGE, LN_TO_L10
  use mod_mol_block
  implicit none
  private
  public :: tree
  public :: tree_tuple
  public :: tree_memsize
  public :: tree_current_level
  public :: tree_root_pointer
  public :: tree_queue_pointer
  public :: tree_current_pointer
  public :: tree_node_pointer
  public :: tree_current_sequence
  public :: tree_current_permutation
  public :: tree_current_mapping
  public :: tree_expand
  public :: tree_leave
  public :: tree_select_top_node
  public :: tree_n_perm
  public :: tree_n_depth
  public :: tree_log_ncomb
  public :: tree_ncomb_frac
  public :: tree_ncomb_exp
  public :: tree_queue_is_empty
  public :: tree_queue_is_explored
  public :: tree_queue_is_unexplored
  public :: tree_queue_is_bottom
!
!| Factorial tree.
  type tree
    private
    sequence
!|  x :: pointer to work array.
    integer(IK), public :: x
!|  q :: pointer to queue array.
    integer(IK), public :: q
  end type tree
!
!| A set of tree and work arrays. <br>
!  This is mainly used for passing during initialization.
  type tree_tuple
    type(tree)               :: t
!!  tree.
    integer(IK), allocatable :: q(:)
!!  queue, a collection of nodes in a hierarchy.
    real(RK), allocatable    :: x(:)
!!  main memory.
  contains
    final     :: tree_tuple_destroy
  end type tree_tuple
!
  interface tree_tuple
    module procedure tree_tuple_new
  end interface tree_tuple
!
! actural pointer is calculated by t%q + q[sdlcr]
  integer(IK), parameter :: header_blocksize = 4
  integer(IK), parameter :: qs = 0 ! scaling.
  integer(IK), parameter :: qd = 1 ! tree depth.
  integer(IK), parameter :: ql = 2 ! current level.
  integer(IK), parameter :: qc = 3 ! current queue pointer.
  integer(IK), parameter :: qr = 4 ! pointer to root node
!
! actural pointer is calculated by q(t%q + q(t%q + qc))
  integer(IK), parameter :: queue_blocksize = 4
  integer(IK), parameter :: qi = 1
  integer(IK), parameter :: qp = 2
  integer(IK), parameter :: qn = 3
  integer(IK), parameter :: qx = 4
!
  integer(IK), parameter :: is_unexplored = -1
  integer(IK), parameter :: is_explored = -2
!
contains
!
  pure subroutine queue_set_state(q, i)
    integer(IK), intent(inout) :: q(*)
    integer(IK), intent(in)    :: i
    q(qi) = i
  end subroutine queue_set_state
!
  pure function queue_state(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = q(qi)
  end function queue_state
!
  pure function queue_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = q(qp)
  end function queue_pointer
!
  pure function queue_state_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = q(qp) + q(qi) * q(qx)
  end function queue_state_pointer
!
  pure function queue_node_pointer(q, i) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK), intent(in) :: i
    integer(IK)             :: res
    res = q(qp) + i * q(qx)
  end function queue_node_pointer
!
  pure function queue_nnodes(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = q(qn)
  end function queue_nnodes
!
  pure function queue_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = q(qn) * q(qx)
  end function queue_memsize
!
  pure function queue_memstride(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = q(qx)
  end function queue_memstride
!
!| Constructer of queue
  pure subroutine queue_init(n_nodes, memsize, p, q)
    integer(IK), intent(in) :: n_nodes
!!  n_nodes :: number of nodes in queue, n_nodes>0.
    integer(IK), intent(in) :: memsize
!!  memsize :: memory size of each node.
    integer(IK), intent(in) :: p
!!  p :: offset of pointer
    integer(IK), intent(inout) :: q(*)
!!  q :: queue array
    q(qi) = -1 ! current state -> unexplored
    q(qp) = p  ! pointer to work array.
    q(qn) = MAX(1, n_nodes)  ! max number of nodes in this queue.
    q(qx) = MAX(1, memsize)  ! memsize of a node.
  end subroutine queue_init
!
!| Constructer of factorial tree.<br>
!  [s*m, s*(m-1),..., s*2, s]
  pure function tree_tuple_new(b, memsize) result(res)
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
    type(tree_tuple)            :: res
    integer(IK)                 :: n, s, i, j, k, l
!
    res%t = tree(1, 1)
!
    n = mol_block_nmol(b)
    s = mol_block_nsym(b)
    allocate (res%q(header_blocksize + queue_blocksize * (n + 1)))
!
!   designate root as the current node.
    res%q(res%t%q + qs) = s
    res%q(res%t%q + qd) = n + 1
    res%q(res%t%q + ql) = 0
    res%q(res%t%q + qc) = header_blocksize ! pointer to root - 1
!
    j = 1
    k = 1
    l = res%t%q + res%q(res%t%q + qc) ! root_pointer
    call queue_init(j, memsize(b, 0), k, res%q(l))
!   set root state => 1
    res%q(l + qi - 1) = 0
!
    j = (n + 1) * s
!
    do i = 1, n
      j = j - s
      k = k + queue_memsize(res%q(l))
      l = l + queue_blocksize
      call queue_init(j, memsize(b, i), k, res%q(l))
    end do
!
    allocate (res%x(tree_memsize(res%t, res%q)))
!
  end function tree_tuple_new
!
!| Inquire memsize of tree.
  pure function tree_memsize(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res, i, l, u
    l = t%q + qr
    u = l + q(t%q + qd) * queue_blocksize
    res = 0
    do i = l, u, queue_blocksize
      res = res + queue_memsize(q(i))
    enddo
  end function tree_memsize
!
!| Returns current level.
  pure function tree_current_level(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res
    res = q(t%q + ql)
  end function tree_current_level
!
!| Returns a pointer to the root node.
  pure function tree_root_pointer(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res
    res = t%x + queue_pointer(q(t%q + qr)) - 1
  end function tree_root_pointer
!
!| Returns a pointer to the current queue.
  pure function tree_queue_pointer(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res
    res = t%q + q(t%q + qc)
    res = t%x + queue_pointer(q(res)) - 1
  end function tree_queue_pointer
!
!| Returns a pointer to the current queue.
  pure function tree_current_pointer(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res
    res = t%q + q(t%q + qc)
    res = t%x + queue_state_pointer(q(res)) - 1
  end function tree_current_pointer
!
!| Returns a pointer to the current best node.
  pure function tree_node_pointer(t, q, iper, imap) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK), intent(in) :: iper
!!  iper :: permutation index, must be [0,1,...,q%n/s-1].
    integer(IK), intent(in) :: imap
!!  imap :: mapping index, must be [0,1,...,s-1].
    integer(IK)             :: res
    res = t%q + q(t%q + qc)
    res = t%x + queue_node_pointer(q(res), (iper * q(t%q + qs) * imap))
  end function tree_node_pointer
!
!| Returns current sequence.
  pure function tree_current_sequence(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: i, j, res(tree_n_depth(t, q))
    j = t%q + qr
    do i = 1, SIZE(res)
      j = j + queue_blocksize
      res(i) = queue_state(q(j)) + 1
    enddo
  end function tree_current_sequence
!
!| Returns current permutation.
  pure function tree_current_permutation(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: i, j, p, s, res(tree_n_depth(t, q))
    do concurrent(i=1:SIZE(res))
      res(i) = i
    end do
    s = q(t%q + qs)
    j = t%q + qr
    do i = 1, SIZE(res)
      j = j + queue_blocksize
      p = queue_state(q(j))
      if (p < 0) return
      if (p < s) cycle
      p = i + p / s
      res(i:p) = [res(p), res(i:p - 1)]
    enddo
  end function tree_current_permutation
!
!| Returns current mapping.
  pure function tree_current_mapping(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: i, j, s, res(tree_n_depth(t, q))
    j = t%q + qr
    s = q(t%q + qs)
    do i = 1, SIZE(res)
      j = j + queue_blocksize
      res(i) = MODULO(queue_state(q(j)), s)
    end do
  end function tree_current_mapping
!
!| Expand current node
  pure subroutine tree_expand(t, q)
    type(tree), intent(in)     :: t
!!  t :: tree
    integer(IK), intent(inout) :: q(*)
!!  q :: queue
    integer(IK)                :: j, k
    if(tree_queue_is_empty(t, q).or.tree_queue_is_bottom(t, q)) return
    j = t%q + ql
    k = t%q + qc
    q(j) = q(j) + 1               ! l = l + 1
    q(k) = q(k) + queue_blocksize ! let current queue be l.
    call queue_set_state(q(t%q + q(k)), is_unexplored)
  end subroutine tree_expand
!
!| Leave current node
  pure subroutine tree_leave(t, q)
    type(tree), intent(in)     :: t
!!  t :: tree
    integer(IK), intent(inout) :: q(*)
!!  q :: queue
    integer(IK)                :: j, k
    if (tree_current_level(t, q) < 1) return
    j = t%q + ql
    k = t%q + qc
    q(j) = q(j) - 1               ! l = l - 1
    q(k) = q(k) - queue_blocksize ! let current queue be l.
  end subroutine tree_leave
!
!| Select top node
  pure subroutine tree_select_top_node(t, q, UB, W)
    type(tree), intent(in)     :: t
!!  t :: tree
    integer(IK), intent(inout) :: q(*)
!!  q :: queue
    real(RK), intent(in)       :: UB
!!  UB :: upperbound
    real(RK), intent(in)       :: W(*)
!!  W :: work array
    real(RK)                   :: uv, lv
    integer(IK)                :: i, j, p, n, b
!
    uv = UB
!
    if (tree_queue_is_explored(t, q)) then
      return
    elseif (tree_queue_is_unexplored(t, q)) then
      lv = -RHUGE
    else
      lv = W(tree_current_pointer(t, q))
    end if
!
    j = t%q + q(t%q + qc)
    call queue_set_state(q(j), is_explored)
!
    if (uv < lv) return
!
    p = tree_queue_pointer(t, q)
    n = queue_nnodes(q(j))
    b = queue_memstride(q(j))
!
    do i = 0, n - 1
      if (lv < W(p) .and. W(p) < uv) then
        call queue_set_state(q(j), i)
        uv = W(p)
      end if
      p = p + b
    end do
!
  end subroutine tree_select_top_node
!
!|  returns number of permutation in q.
  pure function tree_n_perm(t, q) result(res)
    type(tree), intent(in)  :: t
!!  this :: tree.
    integer(IK), intent(in) :: q(*)
!!  q :: queue.
    integer(IK)             :: res
    res = queue_nnodes(q(t%q + q(t%q + qc))) / q(t%q + qs)
  end function tree_n_perm
!
!|  returns number tree depth (without root node).
  pure function tree_n_depth(t, q) result(res)
    type(tree), intent(in)  :: t
!!  this :: tree.
    integer(IK), intent(in) :: q(*)
!!  q :: queue.
    integer(IK)             :: res
    res = q(t%q + qd) - 1
  end function tree_n_depth
!
!|  returns number of nodes in tree.
  pure function tree_log_ncomb(t, q) result(res)
    type(tree), intent(in)  :: t
!!  this :: tree.
    integer(IK), intent(in) :: q(*)
!!  q :: queue.
    real(RK)                :: tmp, res
    integer(IK)             :: i, l, u
!
    res = ZERO
    tmp = ZERO
    l = t%q + qr + qn - 1
    u = l + tree_n_depth(t, q) * queue_blocksize
    do i = u, l, -queue_blocksize
      tmp = tmp - LOG(REAL(q(i), RK))
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
  pure function tree_ncomb_frac(t, q) result(res)
    type(tree), intent(in)  :: t
!!  this :: tree.
    integer(IK), intent(in) :: q(*)
!!  q :: queue.
    real(RK)                :: tmp, res
    tmp = LN_TO_L10 * tree_log_ncomb(t, q)
    res = TEN**(tmp - REAL(INT(tmp), RK))
  end function tree_ncomb_frac
!
!| returns number of nodes in exp.
  pure function tree_ncomb_exp(t, q) result(res)
    type(tree), intent(in)  :: t
!!  this :: tree.
    integer(IK), intent(in) :: q(*)
!!  q :: queue.
    integer(IK)             :: res
!
    res = INT(LN_TO_L10 * tree_log_ncomb(t, q), IK)
!
  end function tree_ncomb_exp
!
  pure function tree_queue_is_unexplored(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    logical                 :: res
    integer(IK)             :: i
    i = t%q + q(t%q + qc)
    res = queue_state(q(i)) == -1
  end function tree_queue_is_unexplored
!
  pure function tree_queue_is_explored(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    logical                 :: res
    integer(IK)             :: i
    i = t%q + q(t%q + qc)
    res = queue_state(q(i)) < -1
  end function tree_queue_is_explored
!
  pure function tree_queue_is_empty(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    logical                 :: res
    integer(IK)             :: i
    i = t%q + q(t%q + qc)
    res = queue_state(q(i)) < 0
  end function tree_queue_is_empty
!
  pure function tree_queue_is_bottom(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    logical                 :: res
    res = q(t%q + ql) == tree_n_depth(t, q)
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
  pure elemental subroutine tree_tuple_destroy(this)
    type(tree_tuple), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%x)) deallocate (this%x)
  end subroutine tree_tuple_destroy
!
end module mod_tree


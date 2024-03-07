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
    integer(IK) :: s
    !! scaling.
    integer(IK) :: d
    !! tree depth.
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
! header pointer.
  integer(IK), parameter :: header_blocksize = 2
  integer(IK), parameter :: ql = 1 ! current level.
  integer(IK), parameter :: qc = 2 ! current queue pointer.
  integer(IK), parameter :: qr = 3 ! pointer to root node
!
! node block pointer is calculated by q(q(qc)) + q[ipnx] - 1
  integer(IK), parameter :: queue_blocksize = 4
  integer(IK), parameter :: qi = 1
  integer(IK), parameter :: qp = 2
  integer(IK), parameter :: qn = 3
  integer(IK), parameter :: qx = 4
!
! state parameter.
  integer(IK), parameter :: is_unexplored = -1
  integer(IK), parameter :: is_explored = -2
!
contains
!
!| queue functions.
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
!| initializer of queue
  pure subroutine queue_init(n_nodes, memsize, p, q)
    integer(IK), intent(in) :: n_nodes
!!  n_nodes :: number of nodes in queue, n_nodes>0.
    integer(IK), intent(in) :: memsize
!!  memsize :: memory size of each node.
    integer(IK), intent(in) :: p
!!  p :: offset of pointer
    integer(IK), intent(inout) :: q(*)
!!  q :: queue array
    q(qi) = is_unexplored    ! current state -> unexplored
    q(qp) = p                ! pointer to work array.
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
    integer(IK)                 :: i, j, k, l
!
    res%t = tree(mol_block_nsym(b), mol_block_nmol(b) + 1)
!
    allocate (res%q(header_blocksize + queue_blocksize * res%t%d))
!
    res%q(ql) = 0  ! designate root as the current node as 0.
    res%q(qc) = qr ! root_pointer
!
    l = res%q(qc) ! root_pointer
!
    call queue_init(1, memsize(b, 0), 1, res%q(qr)) ! construct root with one node.
    call queue_set_state(res%q(qr), 0)              ! set root state as 1.
!
    j = res%t%d * res%t%s
    k = 1
!
    do i = 1, res%t%d - 1
      j = j - res%t%s
      k = k + queue_memsize(res%q(l))
      l = l + queue_blocksize
      call queue_init(j, memsize(b, i), k, res%q(l))
    end do
!
    allocate (res%x(tree_memsize(res%t, res%q)))
!
  end function tree_tuple_new
!
!| Inquire total memsize of tree.
  pure function tree_memsize(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res, i, j
    j = qr
    res = 0
    do i = 1, t%d
      res = res + queue_memsize(q(j))
      j = j + queue_blocksize
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
    res = q(ql)
  end function tree_current_level
!
!| Returns a pointer to the root node.
  pure function tree_root_pointer(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res
    res = queue_pointer(q(qr))
  end function tree_root_pointer
!
!| Returns a pointer to the current queue.
  pure function tree_queue_pointer(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res
    res = queue_pointer(q(q(qc)))
  end function tree_queue_pointer
!
!| Returns a pointer to the current queue.
  pure function tree_current_pointer(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res
    res = queue_state_pointer(q(q(qc)))
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
    res = queue_node_pointer(q(q(qc)), (iper * t%s * imap))
  end function tree_node_pointer
!
!| Returns current sequence.
  pure function tree_current_sequence(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: i, j, res(t%d - 1)
    j = qr
    do i = 1, t%d - 1
      j = j + queue_blocksize
      res(i) = queue_state(q(j)) + 1
    end do
  end function tree_current_sequence
!
!| Returns current permutation.
  pure function tree_current_permutation(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: res(t%d - 1)
    integer(IK)             :: i, j, k, p, s
    do concurrent(i=1:SIZE(res))
      res(i) = i
    end do
    k = qr
    do i = 1, q(ql)
      k = k + queue_blocksize
      p = queue_state(q(k))
      if (p < 0) return
      if (p < t%s) cycle
      ! cyclic swap res(i:p)
      p = i + p / t%s
      s = res(p)
      do j = p, i + 1, -1
        res(j) = res(j - 1)
      end do
      res(i) = s
    enddo
  end function tree_current_permutation
!
!| Returns current mapping.
  pure function tree_current_mapping(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    integer(IK)             :: i, j, res(tree_n_depth(t, q))
    j = qr
    do i = 1, SIZE(res)
      j = j + queue_blocksize
      res(i) = MODULO(queue_state(q(j)), t%s)
    end do
  end function tree_current_mapping
!
!| Expand current node
  pure subroutine tree_expand(t, q)
    type(tree), intent(in)     :: t
!!  t :: tree
    integer(IK), intent(inout) :: q(*)
!!  q :: queue
    if(tree_queue_is_empty(t, q).or.tree_queue_is_bottom(t, q)) return
    q(ql) = q(ql) + 1               ! l = l + 1
    q(qc) = q(qc) + queue_blocksize ! let current queue be l.
    call queue_set_state(q(q(qc)), is_unexplored)
  end subroutine tree_expand
!
!| Leave current node
  pure subroutine tree_leave(t, q)
    type(tree), intent(in)     :: t
!!  t :: tree
    integer(IK), intent(inout) :: q(*)
!!  q :: queue
    if (tree_current_level(t, q) < 1) return
    q(ql) = q(ql) - 1               ! l = l - 1
    q(qc) = q(qc) - queue_blocksize ! let current queue be l.
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
    integer(IK)                :: i, c, p, n, b
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
    call queue_set_state(q(q(qc)), is_explored)
!
    if (uv < lv) return
!
    c = q(qc)
    p = tree_queue_pointer(t, q)
    n = queue_nnodes(q(c))
    b = queue_memstride(q(c))
!
    do i = 0, n - 1
      if (lv < W(p) .and. W(p) < uv) then
        call queue_set_state(q(c), i)
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
    res = queue_nnodes(q(q(qc))) / t%s
  end function tree_n_perm
!
!|  returns number tree depth (without root node).
  pure function tree_n_depth(t, q) result(res)
    type(tree), intent(in)  :: t
!!  this :: tree.
    integer(IK), intent(in) :: q(*)
!!  q :: queue.
    integer(IK)             :: res
    res = t%d - 1
  end function tree_n_depth
!
!|  returns number of nodes in tree.
  pure function tree_log_ncomb(t, q) result(res)
    type(tree), intent(in)  :: t
!!  this :: tree.
    integer(IK), intent(in) :: q(*)
!!  q :: queue.
    real(RK)                :: tmp, res
    integer(IK)             :: i, j
!
    res = ZERO
    tmp = ZERO
    j = qr + (t%d - 1) * queue_blocksize + qn - 1
    do i =  t%d, 1, -1
      tmp = tmp - LOG(REAL(q(j), RK))
      res = res + EXP(tmp)
      j = j - queue_blocksize
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
    res = queue_state(q(q(qc))) == -1
  end function tree_queue_is_unexplored
!
  pure function tree_queue_is_explored(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    logical                 :: res
    res = queue_state(q(q(qc))) < -1
  end function tree_queue_is_explored
!
  pure function tree_queue_is_empty(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    logical                 :: res
    res = queue_state(q(q(qc))) < 0
  end function tree_queue_is_empty
!
  pure function tree_queue_is_bottom(t, q) result(res)
    type(tree), intent(in)  :: t
!!  t :: tree
    integer(IK), intent(in) :: q(*)
!!  q :: queue
    logical                 :: res
    res = q(ql) == (t%d - 1)
  end function tree_queue_is_bottom
!
  pure elemental subroutine tree_tuple_destroy(this)
    type(tree_tuple), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%x)) deallocate (this%x)
  end subroutine tree_tuple_destroy
!
end module mod_tree


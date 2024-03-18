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
  public :: tree_memsize
  public :: tree_current_level
  public :: tree_current_state
  public :: tree_current_iper
  public :: tree_current_isym
  public :: tree_queue_pointer
  public :: tree_current_pointer
  public :: tree_node_pointer
  public :: tree_current_sequence
  public :: tree_current_permutation
  public :: tree_current_mapping
  public :: tree_expand
  public :: tree_leave
  public :: tree_select_top_node
  public :: tree_reset
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
!  This is mainly used for passing during initialization.
  type tree
    integer(IK), allocatable :: q(:)
!!  header list
    integer(IK), allocatable :: s(:)
!!  state list
  contains
    final     :: tree_destroy
  end type tree
!
  interface tree
    module procedure tree_new
  end interface tree
!
! header pointer.
  integer(IK), parameter :: queue_headersize = 2
  integer(IK), parameter :: qs = 1 ! scaling.
  integer(IK), parameter :: qd = 2 ! tree depth.
  integer(IK), parameter :: qr = 3 ! pointer to top node
!
! node block pointer is calculated by q(q(qc)) + q[ipnx] - 1
  integer(IK), parameter :: queue_blocksize = 3
  integer(IK), parameter :: qp = 1
  integer(IK), parameter :: qn = 2
  integer(IK), parameter :: qx = 3
!
  integer(IK), parameter :: state_headersize = 1
  integer(IK), parameter :: sl = 1 ! current level
  integer(IK), parameter :: sr = 2 ! top node state.
!
  integer(IK), parameter :: state_blocksize = 1
  integer(IK), parameter :: ss = 1 ! queue state
!
! state parameter.
  integer(IK), parameter :: is_unexplored = -1
  integer(IK), parameter :: is_explored = -2
!
contains
!
!| Constructer of factorial tree.<br>
!  [s*m, s*(m-1),..., s*2, s]
  pure function tree_new(b, memsize) result(res)
!|  b :: mol_block.
    integer(IK), intent(in) :: b(*)
    interface
      pure function memsize(b, p) result(res)
        use mod_params, only: IK
        use mod_mol_block, only: mol_block
        integer(IK), intent(in) :: b(*)
        integer(IK), intent(in) :: p
        integer(IK)             :: res
      end function memsize
    end interface
    type(tree)              :: res
    integer(IK)             :: i, j, k, l
!
    l = mol_block_nmol(b)
    allocate (res%q(queue_headersize + queue_blocksize * l))
    allocate (res%s(state_headersize + state_blocksize * l))
!
    res%q(qs) = mol_block_nsym(b)
    res%q(qd) = l
!
    j = res%q(qs) * res%q(qd)
    k = 1
    l = qr ! top node pointer
!
    do i = 1, res%q(qd)
      call queue_init(j, memsize(b, i), k, res%q(l))
      j = j - res%q(qs)
      k = k + queue_memsize(res%q(qr), i)
      l = l + queue_blocksize
    end do
!
    call tree_reset(res%q, res%s)
!
  end function tree_new
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
!   q(qi) = is_unexplored    ! current state -> unexplored
    q(qp) = p                ! pointer to work array.
    q(qn) = MAX(1, n_nodes)  ! max number of nodes in this queue.
    q(qx) = MAX(1, memsize)  ! memsize of a node.
  end subroutine queue_init
!
!| reset tree
  pure subroutine tree_reset(q, s)
    integer(IK), intent(in)    :: q(*)
!!  queue
    integer(IK), intent(inout) :: s(*)
!!  state
    integer(IK)                :: i
!
    s(sl) = 1 ! current state
    do concurrent(i=1:q(qd))
      s(sr + state_blocksize * (i - 1) + ss - 1) = is_unexplored ! queue state
    end do
!
  end subroutine tree_reset
!
!| Inquire total memsize of tree
  pure function tree_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK)             :: res, i
    res = 0
    do i = 1, q(qd)
      res = res + queue_memsize(q(qr), i)
    end do
  end function tree_memsize
!
!| Returns current level
  pure function tree_current_level(s) result(res)
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = s(sl)
  end function tree_current_level
!
!| Returns current level
  pure function tree_current_state(s) result(res)
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_state(s(sr), s(sl))
  end function tree_current_state
!
!| Returns current permutation
  pure function tree_current_iper(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_state(s(sr), s(sl)) / q(qs)
  end function tree_current_iper
!
!| Returns current mapping
  pure function tree_current_isym(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = MODULO(queue_state(s(sr), s(sl)), q(qs))
  end function tree_current_isym
!
!| Returns a pointer to the current queue.
  pure function tree_queue_pointer(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_pointer(q(qr), s(sl))
  end function tree_queue_pointer
!
!| Returns a pointer to the current queue.
  pure function tree_current_pointer(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_node_pointer(q(qr), s(sl), queue_state(s(sr), s(sl)))
  end function tree_current_pointer
!
!| Returns a pointer to the current best node.
  pure function tree_node_pointer(q, s, iper, isym) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK), intent(in) :: iper
!!  permutation index, must be [0,1,...,q%n/s-1].
    integer(IK), intent(in) :: isym
!!  mapping index, must be [0,1,...,s-1].
    integer(IK)             :: res
    res = queue_node_pointer(q(qr), s(sl), iper * q(qs) + isym)
  end function tree_node_pointer
!
!| Returns current sequence.
  pure function tree_current_sequence(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: i, res(q(qd))
    do concurrent(i=1:q(qd))
      res(i) = queue_state(s(sr), i) + 1
    end do
  end function tree_current_sequence
!
!| Returns current permutation.
  pure function tree_current_permutation(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res(q(qd))
    integer(IK)             :: i, j, p, t
    do concurrent(i=1:q(qd))
      res(i) = i
    end do
    do i = 1, q(qd)
      p = queue_state(s(sr), i)
      if (p < 0) return
      if (p < q(qs)) cycle
      ! cyclic swap res(i:p)
      p = i + p / q(qs)
      t = res(p)
      do j = p, i + 1, -1
        res(j) = res(j - 1)
      end do
      res(i) = t
    enddo
  end function tree_current_permutation
!
!| Returns current mapping.
  pure function tree_current_mapping(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: i, res(q(qd))
    do concurrent(i=1:q(qd))
      res(i) = MODULO(queue_state(s(sr), i), q(qs))
    end do
  end function tree_current_mapping
!
!| Expand current node
  pure subroutine tree_expand(q, s)
    integer(IK), intent(in)    :: q(*)
!!  queue
    integer(IK), intent(inout) :: s(*)
!!  state
     if(tree_queue_is_empty(q, s).or.tree_queue_is_bottom(q, s)) return
     s(sl) = s(sl) + 1
    call set_state(s, is_unexplored)
  end subroutine tree_expand
!
!| Leave current node
  pure subroutine tree_leave(q, s)
    integer(IK), intent(in)    :: q(*)
!!  queue
    integer(IK), intent(inout) :: s(*)
!!  state
    if (s(sl) < 2) return
    s(sl) = s(sl) - 1
  end subroutine tree_leave
!
!| Select top node
  pure subroutine tree_select_top_node(q, s, UB, W)
    integer(IK), intent(in)    :: q(*)
!!  queue
    integer(IK), intent(inout) :: s(*)
!!  state
    real(RK), intent(in)       :: UB
!!  upperbound
    real(RK), intent(in)       :: W(*)
!!  work array
    real(RK)                   :: uv, lv
    integer(IK)                :: i, p, n, b
!
    uv = UB
!
    if (tree_queue_is_explored(q, s)) then
      return
    elseif (tree_queue_is_unexplored(q, s)) then
      lv = -RHUGE
    else
      lv = W(tree_current_pointer(q, s))
    end if
!
    call set_state(s, is_explored)
!
    if (uv < lv) return
!
    p = tree_queue_pointer(q, s)
    n = queue_nnodes(q(qr), s) - 1
    b = queue_memstride(q(qr), s)
!
    do i = 0, n
      if (lv < W(p) .and. W(p) < uv) then
        call set_state(s, i)
        uv = W(p)
      end if
      p = p + b
    end do
!
  end subroutine tree_select_top_node
!
!|  returns number of permutation in q.
  pure function tree_n_perm(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_nnodes(q(qr), s(sl)) / q(qs)
  end function tree_n_perm
!
!|  returns number tree depth (without root node).
  pure function tree_n_depth(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK)             :: res
    res = q(qd)
  end function tree_n_depth
!
!|  returns number of nodes in tree.
  pure function tree_log_ncomb(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    real(RK)                :: tmp, res
    integer(IK)             :: i, j
!
    res = ZERO
    tmp = ZERO
    j = qr + (q(qd) - 1) * queue_blocksize + qn - 1
    do i =  q(qd), 1, -1
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
  pure function tree_ncomb_frac(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    real(RK)                :: tmp, res
    tmp = LN_TO_L10 * tree_log_ncomb(q)
    res = TEN**(tmp - REAL(INT(tmp), RK))
  end function tree_ncomb_frac
!
!| returns number of nodes in exp.
  pure function tree_ncomb_exp(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK)             :: res
!
    res = INT(LN_TO_L10 * tree_log_ncomb(q), IK)
!
  end function tree_ncomb_exp
!
  pure function tree_queue_is_unexplored(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = queue_state(s(sr), s(sl)) == -1
  end function tree_queue_is_unexplored
!
  pure function tree_queue_is_explored(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = queue_state(s(sr), s(sl)) < -1
  end function tree_queue_is_explored
!
  pure function tree_queue_is_empty(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = queue_state(s(sr), s(sl)) < 0
  end function tree_queue_is_empty
!
  pure function tree_queue_is_bottom(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  queue
    integer(IK), intent(in) :: s(*)
!!  current level
    logical                 :: res
    res = s(sl) == q(qd)
  end function tree_queue_is_bottom
!
  pure elemental subroutine tree_destroy(this)
    type(tree), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine tree_destroy
!
! module functions
!
  pure function queue_state(r, l) result(res)
    integer(IK), intent(in) :: r(state_blocksize, *)
    integer(IK), intent(in) :: l
    integer(IK)             :: res
    res = r(ss, l)
  end function queue_state
!
  pure subroutine set_state(s, t)
    integer(IK), intent(inout) :: s(*)! state vector
    integer(IK), intent(in)    :: t   ! state
    integer(IK)                :: i
    i = state_blocksize * (s(sl) - 1) + sr + ss - 1
    s(i) = t
  end subroutine set_state
!
  pure function queue_pointer(r, s) result(res)
    integer(IK), intent(in) :: r(queue_blocksize, *)
    integer(IK), intent(in) :: s(*)
    integer(IK)             :: res
    res = r(qp, s(sl))
  end function queue_pointer
!
  pure function queue_node_pointer(r, s, n) result(res)
    integer(IK), intent(in) :: r(queue_blocksize, *)
    integer(IK), intent(in) :: s(*)
    integer(IK), intent(in) :: n
    integer(IK)             :: res
    res = r(qp, s(sl)) + r(qx, s(sl)) * n
  end function queue_node_pointer
!
  pure function queue_nnodes(r, s) result(res)
    integer(IK), intent(in) :: r(queue_blocksize, *)
    integer(IK), intent(in) :: s(*)
    integer(IK)             :: res
    res = r(qn, s(sl))
  end function queue_nnodes
!
  pure function queue_memsize(r, l) result(res)
    integer(IK), intent(in) :: r(queue_blocksize, *)
    integer(IK), intent(in) :: l
    integer(IK)             :: res
    res = r(qn, l) * r(qx, l)
  end function queue_memsize
!
  pure function queue_memstride(r, s) result(res)
    integer(IK), intent(in) :: r(queue_blocksize, *)
    integer(IK), intent(in) :: s(*)
    integer(IK)             :: res
    res = r(qx, s(sl))
  end function queue_memstride
!
end module mod_tree


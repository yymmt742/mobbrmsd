!| Factorial tree for branch and bound. <br>
!   This structure handles factorial trees with nodes of constant memory size, ld. <br>
!
!   The factorial tree is completely defined by the number of level \(M\in\mathbb Z\) and a constant \(S\in\mathbb Z\).
!   A node in level \(p\in \{1,2,\dots,M\}\) has \((M-p)S\) childs,
!   and each level has \(S^p \prod_{i=1}^p (M-i+1)\) nodes.
!   Here, only \(M-p+1\) nodes for level \(p=1,2,\dots,M\) are kept as state vector.
!   All other information is discarded. <br>
!
!   The tree informations are kept by header q(\*) and state vector s(\*). <br>
!
!   @note
!   q(\*) constant and should be treated as an immutable variable. <br>
!   However, it is defined as an integer array
!   so that it can be treated as part of a data structure
!   at a higher level imprementation.
!   @endnote
!
!   The tree type has states \(p\) and \(\mathbf{q}^{(p)}\),
!   where \(p=1,2,\dots,M\) is the current level
!   and \(\mathbf{q}^{(p)}=\{q_1^{(p)},q_2^{(p)},\dots,q_p^{(p)}\}\) is the selected node indices.
!   Let \(q^*:=q_p^{(p)}\) be current node,
!   and a retained nodes belonging to \(p\) is defined as a \(p\)-queue. <br>
!
!   A \(p\)-queue has one of the following states,
!   corresponding to the state of the current node. <br>
!   - Node \(i\) is selected :: \(q_i^{(p)} \in [0,1,\dots,(M-p+1)S-1]\)<br>
!   - \(p\)-queue is unexplored :: \(q_i^{(p)} = -1\)<br>
!   - \(p\)-queue is explored :: \(q_i^{(p)} < -1\)<br>
!
!   The real data is kept in an external heap in the form W(ld, \*).
!   The value stored in W(1, \*) is treated as a node evaluation value.
!   (to be implemented by the user). <br>
!
!   @note
!   W can be updated dynamically,
!   but changes to the memory referenced by nodes
!   in the hierarchy below \(p\) will destroy the tree structure.
!   The user should control when to update W. <br>
!   @endnote
!
module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, TEN => RTEN, ZERO => RZERO, RHUGE, LN_TO_L10
  implicit none
  private
  public :: tree
  public :: tree_nnodes
  public :: tree_current_level
  public :: tree_current_state
  public :: tree_current_nnodes
  public :: tree_current_iper
  public :: tree_current_isym
  public :: tree_queue_pointer
  public :: tree_current_pointer
  public :: tree_node_pointer
  public :: tree_current_sequence
  public :: tree_current_permutation
  public :: tree_sequence_to_permutation
  public :: tree_current_mapping
  public :: tree_sequence_to_mapping
  public :: tree_expand
  public :: tree_leave
  public :: tree_select_top_node
  public :: tree_lowest_value
  public :: tree_reset
  public :: tree_n_sym
  public :: tree_n_perm
  public :: tree_n_depth
  public :: tree_log_ncomb
  public :: tree_ncomb_frac
  public :: tree_ncomb_exp
  public :: tree_is_empty
  public :: tree_is_unexplored
  public :: tree_queue_is_empty
  public :: tree_queue_is_left
  public :: tree_queue_is_explored
  public :: tree_queue_is_unexplored
  public :: tree_queue_is_root
  public :: tree_queue_is_bottom
!
!| Factorial tree.<br>
!   @note
!   This type is mainly used for passing during initialization.
!   @endnote
  type tree
    integer(IK), allocatable :: q(:)
!!  header
    integer(IK), allocatable :: s(:)
!!  state
  contains
    final     :: tree_destroy
  end type tree
!
! Constructer
  interface tree
    module procedure tree_new
  end interface tree
!
! --- pointers.
  integer(IK), parameter :: queue_headersize = 2
  integer(IK), parameter :: qs = 1 ! scaling.
  integer(IK), parameter :: qd = 2 ! tree depth.
  integer(IK), parameter :: qr = 3 ! pointer to top node
!
  integer(IK), parameter :: queue_blocksize = 2
  integer(IK), parameter :: qp = 1
  integer(IK), parameter :: qn = 2
!
  integer(IK), parameter :: state_headersize = 1
  integer(IK), parameter :: sl = 1 ! current level
  integer(IK), parameter :: sr = 2 ! top node state.
!
  integer(IK), parameter :: state_blocksize = 1
  integer(IK), parameter :: ss = 1 ! queue state
!
! ---state parameter.
  integer(IK), parameter :: is_unexplored = -1
  integer(IK), parameter :: is_explored = -2
!
contains
!
!| Constructer
  pure function tree_new(nmol, nsym) result(res)
    integer(IK), intent(in) :: nmol
    !! number of molecule
    integer(IK), intent(in) :: nsym
    !! number of molecular symmetry
    type(tree)              :: res
    integer(IK)             :: i, j, k, l
!
    allocate (res%q(queue_headersize + queue_blocksize * nmol))
    allocate (res%s(state_headersize + state_blocksize * nmol))
!
    res%q(qs) = nsym
    res%q(qd) = nmol
!
    j = res%q(qs) * res%q(qd)
    k = 1
    l = qr ! top node pointer
!
    do i = 1, res%q(qd)
      call queue_init(j, k, res%q(l))
      k = k + j
      j = j - res%q(qs)
      l = l + queue_blocksize
    end do
!
    call tree_reset(res%q, res%s)
!
  end function tree_new
!
!| Initializer of queue
  pure subroutine queue_init(n_nodes, p, q)
    integer(IK), intent(in) :: n_nodes
!!  n_nodes :: number of nodes in queue, n_nodes>0.
    integer(IK), intent(in) :: p
!!  p :: offset of pointer
    integer(IK), intent(inout) :: q(*)
!!  q :: queue array
!   q(qi) = is_unexplored    ! current state -> unexplored
    q(qp) = p                ! pointer to work array.
    q(qn) = MAX(1, n_nodes)  ! max number of nodes in this queue.
  end subroutine queue_init
!
!| Reset the state
  pure subroutine tree_reset(q, s)
    integer(IK), intent(in)    :: q(*)
!!  header
    integer(IK), intent(inout) :: s(*)
!!  state
    integer(IK)                :: i
!
    s(sl) = 1 ! current level
    do concurrent(i=1:q(qd))
      s(sr + state_blocksize * (i - 1) + ss - 1) = is_unexplored ! queue state
    end do
!
  end subroutine tree_reset
!
!| Inquire number of nodes,
!  defined by \(S\sum_{i=1}^{M}(M-i+1)=SM(M+1)/2\).
  pure function tree_nnodes(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK)             :: res, i
    res = 0
    do i = 1, q(qd)
      res = res + queue_nnodes(q(qr), i)
    end do
  end function tree_nnodes
!
!| Returns current level, \(p\).
  pure function tree_current_level(s) result(res)
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = s(sl)
  end function tree_current_level
!
!| Returns the state of current node, \(q^*\).
  pure function tree_current_state(s) result(res)
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_state(s(sr), s(sl))
  end function tree_current_state
!
!| Returns the number of nodes in \(p\)-queue.
  pure function tree_current_nnodes(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_nnodes(q(qr), s(sl))
  end function tree_current_nnodes
!
!| Returns current permutation indices of \(q^*\).
  pure function tree_current_iper(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_state(s(sr), s(sl)) / q(qs)
  end function tree_current_iper
!
!| Returns current mapping index of \(q^*\).
  pure function tree_current_isym(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = MODULO(queue_state(s(sr), s(sl)), q(qs))
  end function tree_current_isym
!
!| Returns a pointer to the current queue.
  pure function tree_queue_pointer(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_pointer(q(qr), s(sl))
  end function tree_queue_pointer
!
!| Returns a pointer to the current node, \(q^*\).
  pure function tree_current_pointer(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_node_pointer(q(qr), s, queue_state(s(sr), s(sl)))
  end function tree_current_pointer
!
!| Returns a pointer to a node specified by iper, isym.
  pure function tree_node_pointer(q, s, iper, isym) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK), intent(in) :: iper
!!  permutation index, must be [0,1,...,q%n/s-1].
    integer(IK), intent(in) :: isym
!!  mapping index, must be [0,1,...,s-1].
    integer(IK)             :: res
    res = queue_node_pointer(q(qr), s, iper * q(qs) + isym)
  end function tree_node_pointer
!
!| Returns current sequence indices.
  pure function tree_current_sequence(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: i, res(q(qd))
    do concurrent(i=1:q(qd))
      res(i) = queue_state(s(sr), i)
    end do
  end function tree_current_sequence
!
!| Returns current permutation indices.
  pure function tree_current_permutation(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res(q(qd))
    integer(IK)             :: i, j, p, t
    do concurrent(i=1:q(qd))
      res(i) = i
    end do
    do i = 1, s(sl)
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
    end do
  end function tree_current_permutation
!
!| Convert sequence to permutation.
  pure function tree_sequence_to_permutation(q, z) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: z(*)
!!  state
    integer(IK)             :: res(q(qd))
    integer(IK)             :: i, j, p, t
    do concurrent(i=1:q(qd))
      res(i) = i
    end do
    do i = 1, q(qd)
      if (z(i) < 0) return
      if (z(i) < q(qs)) cycle
      ! cyclic swap res(i:p)
      p = i + z(i) / q(qs)
      t = res(p)
      do j = p, i + 1, -1
        res(j) = res(j - 1)
      end do
      res(i) = t
    end do
  end function tree_sequence_to_permutation
!
!| Returns current mapping indices.
  pure function tree_current_mapping(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: i, res(q(qd))
    do concurrent(i=1:q(qd))
      res(i) = MODULO(queue_state(s(sr), i), q(qs))
    end do
  end function tree_current_mapping
!
!| Convert sequence to mapping.
  pure function tree_sequence_to_mapping(q, z) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: z(*)
!!  state
    integer(IK)             :: i, res(q(qd))
    do concurrent(i=1:q(qd))
      res(i) = MODULO(z(i), q(qs))
    end do
  end function tree_sequence_to_mapping
!
!| Expand current node.
  pure subroutine tree_expand(q, s)
    integer(IK), intent(in)    :: q(*)
!!  header
    integer(IK), intent(inout) :: s(*)
!!  state
    if (tree_queue_is_bottom(q, s)) return
    s(sl) = s(sl) + 1
    call set_state(s, is_unexplored)
  end subroutine tree_expand
!
!| Leave current node.
  pure subroutine tree_leave(q, s)
    integer(IK), intent(in)    :: q(*)
!!  header
    integer(IK), intent(inout) :: s(*)
!!  state
    if (s(sl) > 1) s(sl) = s(sl) - 1
  end subroutine tree_leave
!
!| Select top node, using W(1, \*).
  pure subroutine tree_select_top_node(q, s, ld, UB, W)
    integer(IK), intent(in)    :: q(*)
!!  header
    integer(IK), intent(inout) :: s(*)
!!  state
    integer(IK), intent(in)    :: ld
!!  leading dimension
    real(RK), intent(in)       :: UB
!!  upperbound
    real(RK), intent(in)       :: W(ld, *)
!!  work array
    real(RK)                   :: uv, lv
    integer(IK)                :: i, p1, pn, cp
!
    uv = UB
!
    if (tree_queue_is_explored(q, s)) then
      return
    elseif (tree_queue_is_unexplored(q, s)) then
      cp = tree_queue_pointer(q, s) - 1
      lv = -RHUGE
    else
      cp = tree_current_pointer(q, s)
      lv = W(1, cp)
    end if
!
    call set_state(s, is_explored)
!
    if (uv < lv) return
!
    p1 = tree_queue_pointer(q, s)
    pn = p1 + queue_nnodes(q(qr), s(sl)) - 1
!
    do i = p1, cp - 1
      if (lv < W(1, i) .and. W(1, i) < uv) then
        call set_state(s, i - p1)
        uv = W(1, i)
      end if
    end do
!
    do i = cp + 1, pn
      if (lv < W(1, i) .and. W(1, i) < uv) then
        call set_state(s, i - p1)
        uv = W(1, i)
      end if
    end do
!
  end subroutine tree_select_top_node
!
!| Returns the minimum value of the surviving nodes, excluding the current value.
!  If tree is empty, Returns -infty.
  pure function tree_lowest_value(q, s, ld, W) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK), intent(in) :: ld
!!  leading dimension
    real(RK), intent(in)    :: W(ld, *)
!!  work array
    real(RK)                :: res
    integer(IK)             :: i
    if (tree_is_unexplored(q, s)) then
      res = -RHUGE
      return
    end if
    associate (p => s(sl))
      res = RHUGE
      do i = 1, p
        call queue_second_value(q(qr), s(sr), i, ld, W, res)
      end do
    end associate
  end function tree_lowest_value
!
  pure subroutine queue_second_value(r, s, i, ld, W, res)
    integer(IK), intent(in) :: r(queue_blocksize, *), s(*), i, ld
    real(RK), intent(in)    :: W(ld, *)
    real(RK), intent(inout) :: res
    integer(IK)             :: p, c, u
    associate (l => r(qp, i))
      u = l + r(qn, i) - 1
      c = l + s(i)
      do p = l, c - 1
        if (W(1, p) < W(1, c)) cycle
        res = MIN(W(1, p), res)
      end do
      do p = c + 1, u
        if (W(1, p) < W(1, c)) cycle
        res = MIN(W(1, p), res)
      end do
    end associate
  end subroutine queue_second_value
!
!| Returns number of symmetry in queue.
  pure function tree_n_sym(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK)             :: res
    res = q(qs)
  end function tree_n_sym
!
!| Returns number of permutation in queue.
  pure function tree_n_perm(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    integer(IK)             :: res
    res = queue_nnodes(q(qr), s(sl)) / q(qs)
  end function tree_n_perm
!
!| Returns number tree depth (without root node).
  pure function tree_n_depth(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK)             :: res
    res = q(qd)
  end function tree_n_depth
!
!| Returns number of nodes in tree.
  pure function tree_log_ncomb(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    real(RK)                :: tmp, res
    integer(IK)             :: i, j
!
    res = ZERO
    tmp = ZERO
    j = qr + (q(qd) - 1) * queue_blocksize + qn - 1
    do i = q(qd), 1, -1
      tmp = tmp - LOG(real(q(j), RK))
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
!| Returns number of nodes in fraction.
  pure function tree_ncomb_frac(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    real(RK)                :: tmp, res
    tmp = LN_TO_L10 * tree_log_ncomb(q)
    res = TEN**(tmp - real(INT(tmp), RK))
  end function tree_ncomb_frac
!
!| Returns number of nodes in exp.
  pure function tree_ncomb_exp(q) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK)             :: res
!
    res = INT(LN_TO_L10 * tree_log_ncomb(q), IK)
!
  end function tree_ncomb_exp
!
!| Returns true if current node is unexplored.
  pure function tree_queue_is_unexplored(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = queue_state(s(sr), s(sl)) == is_unexplored
  end function tree_queue_is_unexplored
!
!| Returns true if current node is explored.
  pure function tree_queue_is_explored(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = queue_state(s(sr), s(sl)) < is_unexplored
  end function tree_queue_is_explored
!
!| Returns true if \(p\)-queue is empty.
  pure function tree_is_empty(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = tree_queue_is_empty(q, s) .and. tree_queue_is_bottom(q, s)
  end function tree_is_empty
!
!| Returns true if tree is unexplored.
  pure function tree_is_unexplored(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = tree_queue_is_root(q, s) .and. tree_queue_is_root(q, s)
  end function tree_is_unexplored
!
!| Returns true if \(p\)-queue is explored.
  pure function tree_queue_is_empty(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = queue_state(s(sr), s(sl)) < 0
  end function tree_queue_is_empty
!
!| Returns true if \(p\)-queue has current node.
  pure function tree_queue_is_left(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = queue_state(s(sr), s(sl)) >= 0
  end function tree_queue_is_left
!
!| Returns true if \(p=1\).
  pure function tree_queue_is_root(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = s(sl) == 1
  end function tree_queue_is_root
!
!| Returns true if \(p=M\).
  pure function tree_queue_is_bottom(q, s) result(res)
    integer(IK), intent(in) :: q(*)
!!  header
    integer(IK), intent(in) :: s(*)
!!  state
    logical                 :: res
    res = s(sl) == q(qd)
  end function tree_queue_is_bottom
!
!| Destoructer
  pure elemental subroutine tree_destroy(this)
    type(tree), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine tree_destroy
!
! module functions
!
  pure function queue_state(s, l) result(res)
    integer(IK), intent(in) :: s(state_blocksize, *)
    integer(IK), intent(in) :: l
    integer(IK)             :: res
    res = s(ss, l)
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
    res = r(qp, s(sl)) + n
  end function queue_node_pointer
!
  pure function queue_nnodes(r, i) result(res)
    integer(IK), intent(in) :: r(queue_blocksize, *)
    integer(IK), intent(in) :: i
    integer(IK)             :: res
    res = r(qn, i)
  end function queue_nnodes
!
end module mod_tree


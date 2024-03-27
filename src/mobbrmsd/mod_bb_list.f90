!| mod_bb_list
module mod_bb_list
  use mod_params, only: DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_bb_block
  use mod_rotation
  implicit none
  private
  public :: bb_list
  public :: bb_list_memsize
  public :: bb_list_n_block
  public :: bb_list_n_atoms
  public :: bb_list_log_n_nodes
  public :: bb_list_setup
  public :: bb_list_run
  public :: bb_list_swap_y
  public :: bb_list_rotation_matrix
!
  integer(IK), parameter :: header_size = 1
  integer(IK), parameter :: nb = 1 ! number of block
!
  integer(IK), parameter :: header_sttsize = 1
  integer(IK), parameter :: sb = 1
  integer(IK), parameter :: ss = 2 ! pointer to best state vector
!
  integer(IK), parameter :: header_memsize = 4
  integer(IK), parameter :: ub = 1
  integer(IK), parameter :: lb = 2
  integer(IK), parameter :: nc = 3
  integer(IK), parameter :: rt = 4
!
!| bb_list<br>
!  This is mainly used for passing during initialization.
  type bb_list
    integer(IK), allocatable :: q(:)
    !! integer array
    integer(IK), allocatable :: s(:)
    !! work integer array
  contains
    final           :: bb_list_destroy
  end type bb_list
!
  interface bb_list
    module procedure bb_list_new
  end interface bb_list
!
contains
!
  pure function n_block(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = q(nb)
  end function n_block
!
  pure function s_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = header_size + q(nb) + 1
  end function s_pointer
!
  pure function w_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = header_size + 2 * q(nb) + 1
  end function w_pointer
!
  pure function x_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = header_size + 3 * q(nb) + 1
  end function x_pointer
!
  pure function q_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = header_size + 1
  end function q_pointer
!
!| generate node instance
  pure function bb_list_new(blk) result(res)
    type(bb_block), intent(in) :: blk(:)
    type(bb_list)     :: res
    integer(IK)                :: q(header_size), s(header_sttsize)
    integer(IK)                :: pq(SIZE(blk)), px(SIZE(blk)), ps(SIZE(blk)), pw(SIZE(blk))
    integer(IK)                :: i, j, nstat
!
    q(nb) = SIZE(blk)
    nstat = SUM([(bb_block_nmol(blk(i)%q), i=1, q(nb))])
    s(sb) = 1
!
    j = header_size + 4 * q(nb) + 1
    do i = 1, q(nb)
      pq(i) = j
      j = j + SIZE(blk(i)%q)
    end do
!
    j = header_sttsize + nstat + 1
    do i = 1, q(nb)
      ps(i) = j
      j = j + SIZE(blk(i)%s)
    end do
!
    j = header_memsize + 1
    do i = 1, q(nb)
      pw(i) = j
      j = j + bb_block_memsize(blk(i)%q) + bb_block_worksize(blk(i)%q)
    end do
!
    j = 1
    do i = 1, q(nb)
      px(i) = j
      j = j + bb_block_molsize(blk(i)%q)
    end do
!
    allocate (res%q, source=[q, pq, ps, pw, px, [(blk(i)%q, i=1, SIZE(blk))]])
    allocate (res%s, source=[s, [(-1, i=1, nstat)], [(blk(i)%s, i=1, SIZE(blk))]])
!
  end function bb_list_new
!
!| Inquire worksize of f_matrix.
  pure function bb_list_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block.
    integer(IK)             :: res, i, j, n
!
    res = header_memsize
    n = n_block(q)
    j = q_pointer(q)
    do i = 1, n
      res = res + bb_block_memsize(q(q(j))) + bb_block_worksize(q(q(j)))
      j = j + 1
    end do
!
  end function bb_list_memsize
!
!| Returns number of molecular blocks.
  pure function bb_list_n_block(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block.
    integer(IK)             :: res
!
    res = q(nb)
!
  end function bb_list_n_block
!
!| Returns number of total atoms.
  pure function bb_list_n_atoms(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block.
    integer(IK)             :: res, i, j, n
!
    res = 0
    n = n_block(q)
    j = q_pointer(q)
    do i = 1, n
      res = res + bb_block_natm(q(q(j)))
      j = j + 1
    end do
!
  end function bb_list_n_atoms
!
!| Returns the logarithm of the total number of nodes.
  pure function bb_list_log_n_nodes(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block.
    integer(IK)             :: i, j, n
    real(RK)                :: res
!
    res = ZERO
    n = n_block(q)
    j = q_pointer(q)
    do i = 1, n
      res = res + bb_block_log_ncomb(q(q(j)))
      j = j + 1
    end do
!
  end function bb_list_log_n_nodes
!
!| Setup
  pure subroutine bb_list_setup(q, s, X, Y, W)
    integer(IK), intent(in)    :: q(*)
    !! header
    integer(IK), intent(inout) :: s(*)
    !! state
    real(RK), intent(in)       :: X(*)
    !! reference coordinate
    real(RK), intent(in)       :: Y(*)
    !! target coordinate
    real(RK), intent(inout)    :: W(*)
    !! work array
    integer(IK)                :: i, n, ps, pq, px, pw
!
    s(sb) = 1
    W(ub) = RHUGE
    W(lb) = -RHUGE
    W(nc) = ZERO
    W(rt) = ZERO
!
    ps = s_pointer(q)
    px = x_pointer(q)
    pw = w_pointer(q)
    pq = q_pointer(q)
!
    n = n_block(q)
!
    do concurrent(i=0:n - 1)
      call bb_block_setup(q(q(pq + i)), X(q(px + i)), Y(q(px + i)), s(q(ps + i)), W(q(pw + i)), zfill=(i == 0))
    end do
!
    call save_state(n, q(pq), q(ps), q, s)
!
  end subroutine bb_list_setup
!
!| run branch and bound
  pure subroutine bb_list_run(q, s, W, cutoff, difflim, maxeval)
    integer(IK), intent(in)           :: q(*)
    !! header
    integer(IK), intent(inout)        :: s(*)
    !! state
    real(RK), intent(inout)           :: W(*)
    !! work array
    real(RK), intent(in), optional    :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional    :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional :: maxeval
    !! The search ends when ncount exceeds maxiter.
    real(RK)                          :: coff, diff, nlim
    integer(IK)                       :: pq, ps, pw
    integer(IK)                       :: j, n
!
    n = n_block(q)
    pq = q_pointer(q)
    ps = s_pointer(q)
    pw = w_pointer(q)
!
    coff = RHUGE
    diff = ZERO
    nlim = RHUGE
!
    if (PRESENT(cutoff)) coff = MIN(coff, cutoff)
    if (PRESENT(difflim)) diff = MAX(diff, difflim)
    if (PRESENT(maxeval)) nlim = MERGE(real(maxeval, RK), nlim, maxeval >= 0)
!
    call run_bb(n, q(pq), q(ps), q(pw), q, coff, diff, nlim, s(sb), s, W)
!
    W(nc) = ZERO
    do j = 0, n - 1
      W(nc) = W(nc) + bb_block_evaluation_count(W(q(pw + j)))
    end do
!
    W(rt) = LOG(W(nc))
    do j = 0, n - 1
      W(rt) = W(rt) - bb_block_log_ncomb(q(q(pq + j)))
    end do
!
  end subroutine bb_list_run
!
  pure subroutine run_bb(n, pq, ps, pw, q, coff, diff, nlim, b, s, W)
  integer(IK), intent(in)    :: n, pq(n), ps(n), pw(n), q(*)
  real(RK), intent(in)       :: coff, diff, nlim
  integer(IK), intent(inout) :: b, s(*)
  real(RK), intent(inout)    :: W(*)
!
    do
      do
        call bb_block_expand(W(ub), q(pq(b)), s(ps(b)), W(pw(b)))
        if (b == n .or. bb_block_queue_is_empty(q(pq(b)), s(ps(b)))) exit
        b = b + 1
        call bb_block_inheritance(W(ub), q(pq(b)), s(ps(b)), W(pw(b)), &
       &                          q(pq(b - 1)), s(ps(b - 1)), W(pw(b - 1)))
      enddo
!
      if (b == n &
        & .and. .not. bb_block_queue_is_empty(q(pq(b)), s(ps(b))) &
        & .and. bb_block_queue_is_bottom(q(pq(b)), s(ps(b)))) then
        call lowerbound(n, pq, ps, pw, q, s, W)
        block
          real(RK) :: cv
          cv = bb_block_current_value(q(pq(b)), s(ps(b)), W(pw(b)))
          if (W(ub) > cv) then
            W(ub) = cv
            call save_state(n, pq, ps, q, s)
          end if
        end block
      end if
!
      do
        call bb_block_leave(W(ub), q(pq(b)), s(ps(b)), W(pw(b)))
        if (b == 1 .or. .not. bb_block_queue_is_empty(q(pq(b)), s(ps(b)))) exit
        b = b - 1
      enddo
!
      if (b == 1 .and. bb_block_queue_is_empty(q(pq(b)), s(ps(b)))) return
      if (nlim <= W(nc)) return
      if (W(lb) > coff) return
      if (W(ub) - W(lb) < diff) return
    end do
!
  end subroutine run_bb
!
  pure subroutine lowerbound(n, pq, ps, pw, q, s, W)
  integer(IK), intent(in)    :: n, pq(n), ps(n), pw(n)
  integer(IK), intent(in)    :: q(*), s(*)
  real(RK), intent(inout)    :: W(*)
  real(RK)                   :: lv
  integer(IK)                :: b
    lv = RHUGE
    do b = 1, n
      lv = MIN(lv, bb_block_lowest_value(q(pq(b)), s(ps(b)), W(pw(b))))
    end do
    W(lb) = MAX(W(lb), lv)
  end subroutine lowerbound
!
  pure subroutine save_state(n, pq, ps, q, s)
  integer(IK), intent(in)    :: n, pq(n), ps(n)
  integer(IK), intent(in)    :: q(*)
  integer(IK), intent(inout) :: s(*)
  integer(IK)                :: b, p
    p = ss
    do b = 1, n
      call bb_block_save_state(q(pq(b)), s(ps(b)), s(p))
      p = p + bb_block_nmol(q(pq(b)))
    end do
  end subroutine save_state
!
!| Swap target coordinate.
  pure subroutine bb_list_swap_y(q, s, Y)
    integer(IK), intent(in) :: q(*)
    !! header
    integer(IK), intent(in) :: s(*)
    !! state
    real(RK), intent(inout) :: Y(*)
    !! target coordinate
    integer(IK)             :: i, n, pb, pq, px
!
    px = x_pointer(q)
    pq = q_pointer(q)
    pb = ss
!
    n = n_block(q)
!
    do i = 0, n - 1
      call bb_block_swap_y(q(q(pq + i)), s(pb), Y(q(px + i)))
      pb = pb + bb_block_nmol(q(q(pq + i)))
    end do
!
  end subroutine bb_list_swap_y
!
!| Sum covariance matrix by saved state z.
  pure subroutine bb_list_rotation_matrix(q, s, W, R)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! state
    real(RK), intent(in)    :: W(*)
    !! main memory
    real(RK), intent(inout) :: R(*)
    !! rotation matrix
    real(RK)                :: G, C(DD), V(rotation_worksize())
    integer(IK)             :: i, n, pb, pq, pw
!
    n = n_block(q)
    pb = ss
    pw = w_pointer(q)
    pq = q_pointer(q)
!
    G = ZERO
    C = ZERO
!
    do i = 0, n - 1
      call bb_block_covmat_add(q(q(pq + i)), s(pb), W(q(pw + i)), G, C)
      pb = pb + bb_block_nmol(q(q(pq + i)))
    end do
!
    call estimate_rotation(G, C, R, V)
!
  end subroutine bb_list_rotation_matrix
!
!| destractor
  pure elemental subroutine bb_list_destroy(this)
    type(bb_list), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine bb_list_destroy
!
end module mod_bb_list


!| mod_branch_and_bound
module mod_branch_and_bound
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_bb_block
  implicit none
  private
  public :: branch_and_bound
  public :: branch_and_bound_memsize
  public :: branch_and_bound_setup
  public :: branch_and_bound_run
  public :: DEF_maxeval
  public :: DEF_cutoff
!
  integer(IK), parameter :: DEF_maxeval = -1
  real(RK), parameter    :: DEF_cutoff  = RHUGE
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
!| branch_and_bound<br>
!  This is mainly used for passing during initialization.
  type branch_and_bound
    integer(IK), allocatable :: q(:)
    !! integer array
    integer(IK), allocatable :: s(:)
    !! work integer array
  contains
    final           :: branch_and_bound_destroy
  end type branch_and_bound
!
  interface branch_and_bound
    module procedure branch_and_bound_new
  end interface branch_and_bound
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
  pure function branch_and_bound_new(blk) result(res)
    type(bb_block), intent(in) :: blk(:)
    type(branch_and_bound)     :: res
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
  end function branch_and_bound_new
!
!| Inquire worksize of f_matrix.
  pure function branch_and_bound_memsize(q) result(res)
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
  end function branch_and_bound_memsize
!
!| Setup
  pure subroutine branch_and_bound_setup(q, s, X, Y, W)
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
  end subroutine branch_and_bound_setup
!
  pure subroutine branch_and_bound_run(q, s, W, cutoff, difflim, maxiter)
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
    integer(IK), intent(in), optional :: maxiter
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
    if (PRESENT(cutoff)) coff = MAX(coff, cutoff)
    if (PRESENT(difflim)) diff = MAX(diff, difflim)
    if (PRESENT(maxiter)) nlim = maxiter
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
  end subroutine branch_and_bound_run
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
  pure elemental subroutine branch_and_bound_destroy(this)
    type(branch_and_bound), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine branch_and_bound_destroy
!
end module mod_branch_and_bound


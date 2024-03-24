!| mod_branch_and_bound
module mod_branch_and_bound
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_bb_block
  implicit none
  private
  public :: branch_and_bound
  public :: branch_and_bound_memsize
  public :: branch_and_bound_worksize
  public :: branch_and_bound_setup
  public :: DEF_maxeval
  public :: DEF_cutoff
!
  integer(IK), parameter :: DEF_maxeval = -1
  real(RK), parameter    :: DEF_cutoff  = RHUGE
!
  integer(IK), parameter :: header_size = 1
!
  integer(IK), parameter :: nb = 1
  !! number of block
!
  integer(IK), parameter :: header_memsize = 1
  integer(IK), parameter :: bx = header_memsize + 1
!   ratio = pi; pi = pi + 1
!   nsrch = pi; pi = pi + 1
!   lncmb = pi; pi = pi + 1
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
  pure function s_pointer(q, i) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK), intent(in) :: i
    integer(IK)             :: res
    res = q(q(header_size + i))
  end function s_pointer
!
  pure function w_pointer(q, i) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK), intent(in) :: i
    integer(IK)             :: res
    res = q(q(header_size + i) + 1)
  end function w_pointer
!
  pure function x_pointer(q, i) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK), intent(in) :: i
    integer(IK)             :: res
    res = q(q(header_size + i) + 2)
  end function x_pointer
!
  pure function q_pointer(q, i) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK), intent(in) :: i
    integer(IK)             :: res
    res = q(header_size + i) + 3
  end function q_pointer
!
!| generate node instance
  pure function branch_and_bound_new(blk) result(res)
    type(bb_block), intent(in) :: blk(:)
    type(branch_and_bound)     :: res
    integer(IK)                :: q(header_size)
    integer(IK)                :: p(SIZE(blk))
    integer(IK)                :: r(3, SIZE(blk))
    integer(IK)                :: i, j
!
    q(nb) = SIZE(blk)
!
    j = header_size + SIZE(p) + 1
    do i = 1, SIZE(p)
      p(i) = j
      j = j + 3 + SIZE(blk(i)%q)
    end do
!
    j = 1
    do i = 1, SIZE(r, 2)
      r(1, i) = j
      j = j + SIZE(blk(i)%s)
    end do
!
    j = 1
    do i = 1, SIZE(r, 2)
      r(2, i) = j
      j = j + bb_block_memsize(blk(i)%q)
    end do
!
    j = 1
    do i = 1, SIZE(r, 2)
      r(3, i) = j
      j = j + bb_block_molsize(blk(i)%q)
    end do
!
    allocate (res%q, source=[q, p, [([r(:,i), blk(i)%q], i=1, SIZE(blk))]])
    allocate (res%s, source=[[(blk(i)%s, i=1, SIZE(blk))]])
!
  end function branch_and_bound_new
!
!| Inquire worksize of f_matrix.
  pure function branch_and_bound_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block.
    integer(IK)             :: res, i, p, n
!
    res = 0
    n = n_block(q)
    do i = 1, n
      p = q_pointer(q, i)
      res = res + bb_block_memsize(q(p))
    end do
!
  end function branch_and_bound_memsize
!
!| Inquire worksize of bb_block.
  pure function branch_and_bound_worksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array.
    integer(IK)             :: i, p, n, res, mem
!
    res = 0
    mem = 0
!
    n = n_block(q)
    do i = 1, n
      p = q_pointer(q, i)
      mem = mem + bb_block_memsize(q(p))
      res = MAX(res, mem + bb_block_worksize(q(p)))
    end do
    res = res - branch_and_bound_memsize(q)
!
  end function branch_and_bound_worksize
!
  subroutine branch_and_bound_setup(q, s, X, Y, W)
    integer(IK), intent(in)    :: q(*)
    integer(IK), intent(inout) :: s(*)
    real(RK), intent(in)       :: X(*)
    real(RK), intent(in)       :: Y(*)
    real(RK), intent(inout)    :: W(*)
    integer(IK)                :: i, b, n
!
    b = 1
!
    n = n_block(q)
    do i = 1, n
      block
        integer(IK) :: ps, pq, px, pw
        ps = s_pointer(q, i)
        px = x_pointer(q, i)
        pw = w_pointer(q, i)
        pq = q_pointer(q, i)
        call bb_block_setup(q(pq), X(px), Y(px), s(ps), W(pw))
      end block
    end do
!
  end subroutine branch_and_bound_setup
!
  pure subroutine branch_and_bound_run(q, s, W)
    integer(IK), intent(in)    :: q(*)
    integer(IK), intent(inout) :: s(*)
    real(RK), intent(inout)    :: W(*)
!   logical, intent(in)        :: swap_y
!   integer(IK)                :: cur, pp, cix, ncount
!
!   tr = this%tr
!   bi = this%bi
!   ncount = 0
!
!   call tr%set_parent_node(W)
!
!   do
!     do
!       call tr%open_node()
!       cur = tr%current_depth() - 1
!       call set_hc(this%dx, tr, bi, cur, this%nd, this%bs, W)
!       call tr%prune(W)
!       ncount = ncount + tr%n_breadth()
!
!       if (tr%finished()) exit
!
!       call tr%set_parent_node(W)
!       cix = tr%current_index()
!       bi(cur)%isym = (cix - 1) / bi(cur)%nper
!       call swap_iper(this%nd, cur, cix, bi)
!       if (cur == this%nd) then
!         pp = tr%current_pointer()
!         if (W(tr%upperbound) > W(pp)) then
!           call breadth_indicator_save(bi)
!           call DCOPY(tr%memnode, W(pp), 1, W(tr%ubnode), 1)
!         end if
!         exit
!       end if
!     end do
!
!     call tr%set_lowerbound(W)
!
!     do while (tr%finished() .and. cur > 0)
!       call tr%close_node()
!       call paws_iper(this%nd, cur, bi)
!       call tr%prune(W)
!       cur = cur - 1
!     end do
!
!     if (cur == 0) exit
!     if (this%maxeval > 0 .and. this%maxeval < ncount) exit
!     if (this%cutoff < W(this%lowerbound)) exit
!
!     call tr%set_parent_node(W)
!     block
!       integer(IK) :: cix
!       cix = tr%current_index()
!       bi(cur)%isym = (cix - 1) / bi(cur)%nper
!       call swap_iper(this%nd, cur, cix, bi)
!     end block
!   end do
!
!   W(this%nsrch) = REAL(ncount, RK)
!   if (ncount < 1) then
!     W(this%ratio) = -RHUGE
!   else
!     W(this%ratio) = LOG(W(this%nsrch)) - W(this%lncmb)
!   end if
!
!   if(.not.swap_y) return
!
!   block
!     integer(IK) :: ig, ic
!     ig = tr%ubnode + 1
!     ic = tr%ubnode + 2
!     call rotation(this%mn, this%dmn, W(ig), W(ic), W(this%yp))
!   end block
!
!   block
!     integer(IK) :: i
!     do concurrent(i = 1:this%dx%l)
!       call swap(this%dx%m(i)%m, this%dx%m(i)%g, &
!      &          bi(this%p(i):this%p(i)+this%dx%m(i)%g-1)%jper, &
!      &          bi(this%p(i):this%p(i)+this%dx%m(i)%g-1)%jsym, &
!      &          this%ms(i), W(this%q(i)))
!     end do
!   end block
!
! contains
!
!   pure subroutine rotation(mn, dmn, H, C, Y)
!     integer(IK), intent(in) :: mn, dmn
!     real(RK), intent(in)    :: H, C(d, d)
!     real(RK), intent(inout) :: Y(d, mn)
!     integer(IK)             :: nw
!
!     nw = dd + MAX(dmn, worksize_rotation_matrix())
!
!     block
!       real(RK) :: WL(nw)
!       call estimate_rotation_matrix(H, C, WL(1), WL(dd + 1))
!       call DGEMM('T', 'N', D, mn, D, ONE, WL, D, Y, D, ZERO, WL(dd + 1), D)
!       call DCOPY(dmn, WL(DD + 1), 1, Y, 1)
!     end block
!
!   end subroutine rotation
!
  end subroutine branch_and_bound_run
!
! pure elemental subroutine branch_and_bound_clear(this)
!   class(branch_and_bound), intent(inout) :: this
!   call this%dx%clear()
!   call this%tr%clear()
!   if (ALLOCATED(this%p)) deallocate (this%p)
!   if (ALLOCATED(this%q)) deallocate (this%q)
!   if (ALLOCATED(this%ms)) deallocate (this%ms)
! end subroutine branch_and_bound_clear
!
  pure elemental subroutine branch_and_bound_destroy(this)
    type(branch_and_bound), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine branch_and_bound_destroy
!
!!!
!
! pure elemental subroutine breadth_indicator_save(this)
!   type(breadth_indicator), intent(inout) :: this
!   this%jper = this%iper
!   this%jsym = this%isym
! end subroutine breadth_indicator_save
!
end module mod_branch_and_bound


module mod_branch_and_bound
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_mol_block
  use mod_group_permutation
  use mod_mol_symmetry
  use mod_d_matrix
  use mod_estimate_rotation_matrix
  use mod_tree
  implicit none
  private
  public :: branch_and_bound, DEF_maxeval
!
  integer(IK), parameter :: DEF_maxeval = -1
!
  type breadth_indicator
    sequence
    private
    integer(IK) :: ispc
    integer(IK) :: iper
    integer(IK) :: jper
    integer(IK) :: isym
    integer(IK) :: jsym
    integer(IK) :: nper
    integer(IK) :: nsym
    integer(IK) :: nnod
  end type breadth_indicator
!
  type branch_and_bound
    private
    integer(IK)                :: bs, nd
    integer(IK), public        :: mn, dmn, memsize, maxeval
    integer(IK), public        :: ratio, nsrch, lncmb, xp, yp
    integer(IK), public        :: upperbound, lowerbound
    integer(IK), allocatable   :: p(:), q(:)
    type(d_matrix_list)        :: dx
    type(tree)                 :: tr
    type(breadth_indicator), allocatable :: bi(:)
    type(mol_symmetry), allocatable      :: ms(:)
  contains
    procedure :: setup      => branch_and_bound_setup
    procedure :: run        => branch_and_bound_run
    procedure :: clear      => branch_and_bound_clear
    final     :: branch_and_bound_destroy
  end type branch_and_bound
!
  interface branch_and_bound
    module procedure branch_and_bound_new
  end interface branch_and_bound
!
  interface
    include 'dgemm.h'
    include 'dcopy.h'
  end interface
!
contains
!
!| generate node instance
  pure function branch_and_bound_new(blk, ms, maxeval) result(res)
    type(mol_block_list), intent(in)         :: blk
    type(mol_symmetry), intent(in), optional :: ms(*)
    integer(IK), intent(in), optional        :: maxeval
    type(branch_and_bound)                   :: res
    integer(IK)                              :: i, j, pi
!
    res%mn = blk%mn
    res%dmn = blk%d * blk%mn
!
    pi = 1
    res%ratio = pi; pi = pi + 1
    res%nsrch = pi; pi = pi + 1
    res%lncmb = pi; pi = pi + 1
    res%xp = pi;    pi = pi + res%dmn
    res%yp = pi;    pi = pi + res%dmn
!
    res%dx = d_matrix_list(blk, pi); pi = pi + res%dx%memsize()
!
    res%nd = res%dx%n_depth()
    res%bs = res%dx%dd + 2
!
    allocate (res%bi(res%nd))
!
    do concurrent(j=1:res%dx%l)
      block
        integer(IK) :: k
        k = SUM(res%dx%m(:j - 1)%g)
        do concurrent(i=1:res%dx%m(j)%g)
          block
            integer(IK) :: ib, nper, nnod
            nper = res%dx%m(j)%g - i + 1
            nnod = nper * res%dx%m(j)%s
            ib = k + i
            res%bi(ib) = breadth_indicator(j, i, i, 0, 0, nper, res%dx%m(j)%s, nnod)
          end block
        end do
      end block
    end do
!
    res%tr = tree(pi, res%bs, res%nd + 1, [1, res%bi%nnod]); pi = pi + res%tr%memsize
!
    res%upperbound = res%tr%upperbound
    res%lowerbound = res%tr%lowerbound
!
    if (PRESENT(maxeval)) then
      res%maxeval = maxeval
    else
      res%maxeval = DEF_maxeval
    end if
!
    allocate (res%p(res%dx%l))
    res%p(1) = 1
    do i = 2, res%dx%l
      res%p(i) = res%p(i - 1) + res%dx%m(i - 1)%g
    end do
!
    allocate (res%q(res%dx%l))
    res%q(1) = res%yp
    do i = 2, res%dx%l
      res%q(i) = res%q(i - 1) + res%dx%d * res%dx%m(i - 1)%m * res%dx%m(i - 1)%n
    end do
!
    allocate (res%ms(res%dx%l))
    if (PRESENT(ms)) then
      do concurrent(i=1:res%dx%l)
        res%ms(i) = ms(i)
      end do
    end if
!
    call res%tr%reset()
!
    res%memsize = pi
!
  end function branch_and_bound_new
!
  pure subroutine branch_and_bound_setup(this, X, Y, W)
    class(branch_and_bound), intent(in) :: this
    real(RK), intent(in)                :: X(*)
    real(RK), intent(in)                :: Y(*)
    real(RK), intent(inout)             :: W(*)
    integer(IK)                         :: p
!
    call DCOPY(this%dmn, X, 1, W(this%xp), 1)
    call DCOPY(this%dmn, Y, 1, W(this%yp), 1)
    call this%dx%eval(this%ms, X, Y, W)
!
    p = this%tr%nodes_pointer()
    W(p) = W(this%dx%o)
    p = p + 1
    W(p) = W(this%dx%h)
    p = p + 1
    call DCOPY(this%dx%dd, W(this%dx%c), 1, W(p), 1)
!
    W(this%tr%upperbound) = RHUGE
    W(this%tr%lowerbound) = W(this%dx%o)
    W(this%lncmb) = this%tr%log_ncomb()
!
  end subroutine branch_and_bound_setup
!
  pure subroutine branch_and_bound_run(this, W, swap_y)
    class(branch_and_bound), intent(in)  :: this
    real(RK), intent(inout)              :: W(*)
    logical, intent(in)                  :: swap_y
    type(tree)                           :: tr
    type(breadth_indicator), allocatable :: bi(:)
    integer(IK)                          :: cur, pp, cix, ncount
!
    tr = this%tr
    bi = this%bi
    ncount = 0
!
    call tr%set_parent_node(W)
!
    do
      do
        call tr%open_node()
        cur = tr%current_depth() - 1
        call set_hc(this%dx, tr, bi, cur, this%nd, this%bs, W)
        call tr%prune(W)
        ncount = ncount + tr%n_breadth()
!
        if (tr%finished()) exit
!
        call tr%set_parent_node(W)
        cix = tr%current_index()
        bi(cur)%isym = (cix - 1) / bi(cur)%nper
        call swap_iper(this%nd, cur, cix, bi)
        if (cur == this%nd) then
          pp = tr%current_pointer()
          if (W(tr%upperbound) > W(pp)) then
            call breadth_indicator_save(bi)
            call DCOPY(tr%memnode, W(pp), 1, W(tr%ubnode), 1)
          end if
          exit
        end if
      end do
!
      call tr%set_lowerbound(W)
!
      do while (tr%finished() .and. cur > 0)
        call tr%close_node()
        call paws_iper(this%nd, cur, bi)
        call tr%prune(W)
        cur = cur - 1
      end do
!
      if (this%maxeval > 0 .and. this%maxeval < ncount) exit
      if (cur == 0) exit
!
      call tr%set_parent_node(W)
      block
        integer(IK) :: cix
        cix = tr%current_index()
        bi(cur)%isym = (cix - 1) / bi(cur)%nper
        call swap_iper(this%nd, cur, cix, bi)
      end block
    end do
!
    W(this%nsrch) = REAL(ncount, RK)
    if (ncount < 1) then
      W(this%ratio) = -RHUGE
    else
      W(this%ratio) = LOG(W(this%nsrch)) - W(this%lncmb)
    end if
!
    if(.not.swap_y) return
!
    block
      integer(IK) :: ig, ic
      ig = tr%ubnode + 1
      ic = tr%ubnode + 2
      call rotation(this%dx%d, this%dx%dd, this%mn, this%dmn, W(ig), W(ic), W(this%yp))
    end block
!
    block
      integer(IK) :: i
      do concurrent(i = 1:this%dx%l)
        call swap(this%dx%d, this%dx%m(i)%m, this%dx%m(i)%g, &
       &          bi(this%p(i):this%p(i)+this%dx%m(i)%g-1)%jper, &
       &          bi(this%p(i):this%p(i)+this%dx%m(i)%g-1)%jsym, &
       &          this%ms(i), W(this%q(i)))
      end do
    end block
!
  contains
!
    pure subroutine swap_iper(nd, cur, inod, b)
      integer(IK), intent(in)                :: nd, cur, inod
      type(breadth_indicator), intent(inout) :: b(nd)
      integer(IK)                            :: ip, sw, i
      ip = MODULO(inod - 1, b(cur)%nper)
      sw = b(cur + ip)%iper
      do i = cur + ip - 1, cur, -1
        b(i + 1)%iper = b(i)%iper
      end do
      b(cur)%iper = sw
    end subroutine swap_iper
!
    pure subroutine paws_iper(nd, cur, bi)
      integer(IK), intent(in)                :: nd, cur
      type(breadth_indicator), intent(inout) :: bi(nd)
      integer(IK)                            :: sw, i
      if (cur < 2) return
      do i = cur, nd
        if (bi(i - 1)%ispc /= bi(i)%ispc .or. bi(i - 1)%iper < bi(i)%iper) return
        sw = bi(i - 1)%iper
        bi(i - 1)%iper = bi(i)%iper
        bi(i)%iper = sw
      end do
    end subroutine paws_iper
!
    pure subroutine set_hc(dm, tr, bi, cur, nd, bs, W)
      integer(IK), intent(in)             :: cur, nd, bs
      type(d_matrix_list), intent(in)     :: dm
      type(tree), intent(in)              :: tr
      type(breadth_indicator), intent(in) :: bi(nd)
      real(RK), intent(inout)             :: W(*)
      integer(IK)                         :: iper(nd)
      integer(IK)                         :: i, j, p, q, ph, nx
!
      p = tr%parent_pointer()
      ph = p + 1
      q = tr%nodes_pointer() - bs
      iper = bi%iper
      nx = dm%dd + 1
!
      do concurrent(i=1:bi(cur)%nper, j=0:bi(cur)%nsym - 1)
        block
          integer(IK) :: t, h, c
          t = q + bs * (i + bi(cur)%nper * j)
          h = t + 1
          c = t + 2
          call DCOPY(nx, W(ph), 1, W(h), 1)
          call dm%partial_eval(cur, iper, i, j, W, W(t), W(h), W(c))
        end block
      end do
!
    end subroutine set_hc
!
    pure subroutine swap(d, m, g, iper, isym, ms, X)
      integer(IK), intent(in)             :: d, m, g
      integer(IK), intent(in)             :: iper(g), isym(g)
      type(mol_symmetry), intent(in)      :: ms
      real(RK), intent(inout)             :: X(d, m, g)
      type(group_permutation)             :: gp
      integer(IK)                         :: i
      gp = group_permutation(iper)
      do concurrent(i=1:g)
        call ms%swap(d, X(1, 1, i), isym(i))
      end do
      call gp%reverse(d * m, X)
    end subroutine swap
!
    pure subroutine rotation(d, dd, mn, dmn, H, C, Y)
      integer(IK), intent(in) :: d, dd, mn, dmn
      real(RK), intent(in)    :: H, C(d, d)
      real(RK), intent(inout) :: Y(d, mn)
      real(RK)                :: V(1)
      integer(IK)             :: nw
!
      call estimate_rotation_matrix(-d, V(1), V, V, V)
      nw = dd + MAX(dmn, NINT(V(1)))
!
      block
        real(RK) :: WL(nw)
        call estimate_rotation_matrix(d, H, C, WL(1), WL(dd + 1))
        call DGEMM('T', 'N', d, mn, d, ONE, WL, d, Y, d, ZERO, WL(dd + 1), d)
        call DCOPY(dmn, WL(dd + 1), 1, Y, 1)
      end block
!
    end subroutine rotation
!
  end subroutine branch_and_bound_run
!
  pure elemental subroutine branch_and_bound_clear(this)
    class(branch_and_bound), intent(inout) :: this
    call this%dx%clear()
    call this%tr%clear()
    if (ALLOCATED(this%p)) deallocate (this%p)
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%ms)) deallocate (this%ms)
  end subroutine branch_and_bound_clear
!
  pure elemental subroutine branch_and_bound_destroy(this)
    type(branch_and_bound), intent(inout) :: this
    call branch_and_bound_clear(this)
  end subroutine branch_and_bound_destroy
!
!!!
!
  pure elemental subroutine breadth_indicator_save(this)
    type(breadth_indicator), intent(inout) :: this
    this%jper = this%iper
    this%jsym = this%isym
  end subroutine breadth_indicator_save
!
end module mod_branch_and_bound

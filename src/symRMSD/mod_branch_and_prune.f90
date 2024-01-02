module mod_branch_and_prune
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_mol_block
  use mod_group_permutation
  use mod_mol_symmetry
  use mod_d_matrix
  use mod_tree
  implicit none
  private
  public :: branch_and_prune
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
  type branch_and_prune
    private
    integer(IK)                :: bs, nd, mem
    integer(IK), public        :: dmn
    integer(IK), allocatable   :: p(:), q(:)
    type(d_matrix_list)        :: dm
    type(tree)                 :: tr
    type(breadth_indicator), allocatable :: bi(:)
    type(mol_symmetry), allocatable      :: ms(:)
  contains
    procedure :: memsize    => branch_and_prune_memsize
    procedure :: upperbound => branch_and_prune_upperbound
    procedure :: lowerbound => branch_and_prune_lowerbound
    procedure :: setup      => branch_and_prune_setup
    procedure :: run        => branch_and_prune_run
    procedure :: clear      => branch_and_prune_clear
    final     :: branch_and_prune_destroy
  end type branch_and_prune
!
  interface branch_and_prune
    module procedure branch_and_prune_new
  end interface branch_and_prune
!
  interface
    include 'dcopy.h'
  end interface
!
contains
!
!| generate node instance
  pure function branch_and_prune_new(blk, ms) result(res)
    type(mol_block_list), intent(in)         :: blk
    type(mol_symmetry), intent(in), optional :: ms(*)
    type(branch_and_prune)                   :: res
    integer(IK)                              :: i, j, yp, pi
!
    yp = 1
    pi = yp + blk%d * blk%mn
!
    res%dm = d_matrix_list(blk, pi)
    res%nd = res%dm%n_depth()
    allocate (res%bi(res%nd))
!
    do concurrent(j=1:res%dm%l)
      block
        integer(IK) :: k
        k = SUM(res%dm%m(:j - 1)%g)
        do concurrent(i=1:res%dm%m(j)%g)
          block
            integer(IK) :: ib, nper, nnod
            nper = res%dm%m(j)%g - i + 1
            nnod = nper * res%dm%m(j)%s
            ib = k + i
            res%bi(ib) = breadth_indicator(j, i, i, 0, 0, nper, res%dm%m(j)%s, nnod)
          end block
        end do
      end block
    end do
!
    pi = pi + res%dm%memsize()
!
    res%bs = res%dm%dd + 2
    res%dmn = blk%d * blk%mn
    res%mem = res%dmn
!
    res%tr = tree(pi, res%bs, res%nd + 1, [1, res%bi%nnod])
!
    allocate (res%p(res%dm%l))
    res%p(1) = 1
    do i = 2, res%dm%l
      res%p(i) = res%p(i - 1) + res%dm%m(i - 1)%g
    end do
!
    allocate (res%q(res%dm%l))
    res%q(1) = yp
    do i = 2, res%dm%l
      res%q(i) = res%q(i - 1) + res%dm%d * res%dm%m(i - 1)%m * res%dm%m(i - 1)%n
    end do
!
    allocate (res%ms(res%dm%l))
    if (PRESENT(ms)) then
      do concurrent(i=1:res%dm%l)
        res%ms(i) = ms(i)
      end do
    end if
!
    call res%tr%reset()
!
  end function branch_and_prune_new
!
  pure elemental function branch_and_prune_memsize(this) result(res)
    class(branch_and_prune), intent(in) :: this
    integer(IK)                         :: res
    res = this%dm%memsize() + this%tr%memsize + this%mem
  end function branch_and_prune_memsize
!
  pure subroutine branch_and_prune_setup(this, X, Y, W)
    class(branch_and_prune), intent(in) :: this
    real(RK), intent(in)                :: X(*)
    real(RK), intent(in)                :: Y(*)
    real(RK), intent(inout)             :: W(*)
    integer(IK)                         :: p
!
    call this%dm%eval(this%ms, X, Y, W)
!
    p = this%tr%nodes_pointer()
    W(p) = W(this%dm%o)
    p = p + 1
    W(p) = W(this%dm%h)
    p = p + 1
    call dcopy(this%dm%dd, W(this%dm%c), 1, W(p), 1)
    call dcopy(this%dmn, Y, 1, W(this%q(1)), 1)
!
    W(this%tr%upperbound) = RHUGE
    W(this%tr%lowerbound) = W(this%dm%o)
!
  end subroutine branch_and_prune_setup
!
  pure subroutine branch_and_prune_run(this, W, swap_y)
    class(branch_and_prune), intent(in)  :: this
    real(RK), intent(inout)              :: W(*)
    logical, intent(in)                  :: swap_y
    type(tree)                           :: tr
    type(breadth_indicator), allocatable :: bi(:)
    integer(IK)                          :: cur, pp, cix
!
    tr = this%tr
    bi = this%bi
!
    call tr%set_parent_node(W)
!
    do
      do
        call tr%open_node()
        cur = tr%current_depth() - 1
        call set_hc(this%dm, tr, bi, cur, this%nd, this%bs, W)
        call tr%prune(W)
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
            W(tr%upperbound) = W(pp)
            call breadth_indicator_save(bi)
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
    if(.not.swap_y) return
!
    block
      integer(IK) :: i
      do concurrent(i = 1:this%dm%l)
        call swap(this%dm%d, this%dm%m(i)%m, this%dm%m(i)%g, &
       &          bi(this%p(i):this%p(i)+this%dm%m(i)%g-1)%jper, &
       &          bi(this%p(i):this%p(i)+this%dm%m(i)%g-1)%jsym, &
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
      real(RK)                            :: lb, lf
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
          call dcopy(nx, W(ph), 1, W(h), 1)
          call dm%partial_eval(cur, iper, i, j, W, W(t), W(h), W(c), LB=LB, LF=LF)
          w(t) = LF
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
  end subroutine branch_and_prune_run
!
  pure function branch_and_prune_upperbound(this, W) result(res)
    class(branch_and_prune), intent(in) :: this
    real(RK), intent(in)                :: W(*)
    real(RK)                            :: res
    res = W(this%tr%upperbound)
  end function branch_and_prune_upperbound
!
  pure function branch_and_prune_lowerbound(this, W) result(res)
    class(branch_and_prune), intent(in) :: this
    real(RK), intent(in)                :: W(*)
    real(RK)                            :: res
    res = W(this%tr%lowerbound)
  end function branch_and_prune_lowerbound
!
  pure elemental subroutine branch_and_prune_clear(this)
    class(branch_and_prune), intent(inout) :: this
    call this%dm%clear()
    call this%tr%clear()
    if (ALLOCATED(this%p)) deallocate (this%p)
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%ms)) deallocate (this%ms)
  end subroutine branch_and_prune_clear
!
  pure elemental subroutine branch_and_prune_destroy(this)
    type(branch_and_prune), intent(inout) :: this
    call branch_and_prune_clear(this)
  end subroutine branch_and_prune_destroy
!
!!!
!
  pure elemental subroutine breadth_indicator_save(this)
    type(breadth_indicator), intent(inout) :: this
    this%jper = this%iper
    this%jsym = this%isym
  end subroutine breadth_indicator_save
!
end module mod_branch_and_prune

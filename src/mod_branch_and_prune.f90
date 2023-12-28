module mod_branch_and_prune
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_mol_block
  use mod_group_permutation
  use mod_molecular_rotation
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
    integer(IK)                :: bs, nd
    integer(IK), allocatable   :: p(:), q(:)
    type(d_matrix_list)        :: dm
    type(tree)                 :: tr
    type(breadth_indicator), allocatable  :: bi(:)
    type(molecular_rotation), allocatable :: mr(:)
    type(group_permutation), allocatable  :: gp(:)
  contains
    procedure :: memsize    => branch_and_prune_memsize
    procedure :: upperbound => branch_and_prune_upperbound
    procedure :: setup      => branch_and_prune_setup
    procedure :: run        => branch_and_prune_run
    procedure :: swap       => branch_and_prune_swap
    procedure :: clear      => branch_and_prune_clear
    final     :: branch_and_prune_destroy
  end type branch_and_prune
!
  interface branch_and_prune
    module procedure branch_and_prune_new
  end interface branch_and_prune
!
contains
!
!| generate node instance
  pure function branch_and_prune_new(blk, p, rot) result(res)
    type(mol_block_list), intent(in)     :: blk
    integer(IK), intent(in)              :: p
    type(molecular_rotation), intent(in), optional :: rot(*)
    type(branch_and_prune)               :: res
    integer(IK)                          :: i, j
!
    res%dm = d_matrix_list(blk, p)
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
            res%bi(ib) = breadth_indicator(j, 0, 0, 0, 0, nper, res%dm%m(j)%s, nnod)
          end block
        end do
      end block
    end do
!
    res%bs = res%dm%dd + 2
    res%tr = tree(p + res%dm%memsize(), res%bs, res%nd + 1, [1, res%bi%nnod])
!
    allocate (res%p(res%dm%l))
    res%p(1) = 1
    do i = 2, res%dm%l
      res%p(i) = res%p(i - 1) + res%dm%m(i - 1)%g
    end do
!
    allocate (res%q(res%dm%l))
    res%q(1) = 1
    do i = 2, res%dm%l
      res%q(i) = res%q(i - 1) + res%dm%d * res%dm%m(i - 1)%m * res%dm%m(i - 1)%n
    end do
!
    allocate (res%gp(res%dm%l))
    allocate (res%mr(res%dm%l))
    if (PRESENT(rot)) then
      do concurrent(i=1:res%dm%l)
        res%mr(i) = rot(i)
      end do
    end if
!
  end function branch_and_prune_new
!
  pure elemental function branch_and_prune_memsize(this) result(res)
    class(branch_and_prune), intent(in) :: this
    integer(IK)                         :: res
    res = this%dm%memsize() + this%tr%memsize
  end function branch_and_prune_memsize
!
  pure subroutine branch_and_prune_setup(this, X, Y, W)
    class(branch_and_prune), intent(inout) :: this
    real(RK), intent(in)                   :: X(*)
    real(RK), intent(in)                   :: Y(*)
    real(RK), intent(inout)                :: W(*)
    integer(IK)                            :: p, i, j, k
!
    call this%dm%eval(this%mr, X, Y, W)
    call this%tr%reset()
    W(this%tr%upperbound) = RHUGE
!
    p = this%tr%nodes_pointer()
    W(p) = W(this%dm%o)
    p = p + 1
    W(p) = W(this%dm%h)
    p = p + 1
    call copy(this%dm%dd, W(this%dm%c), W(p))
!
    call this%tr%set_parent_node(W)
!
    k = 0
    do j = 1, this%dm%l
      do i = 1, this%dm%m(j)%g
        k = k + 1
        this%bi(k)%iper = i
        this%bi(k)%jper = i
        this%bi(k)%isym = 0
        this%bi(k)%jsym = 0
      end do
    end do
!
  end subroutine branch_and_prune_setup
!
  subroutine branch_and_prune_run(this, W)
    class(branch_and_prune), intent(inout) :: this
    real(RK), intent(inout)                :: W(*)
    integer(IK)                            :: cur
!
    print*, this%tr%nodes_pointer()
!
    do
      do
        call this%tr%open_node()
        cur = this%tr%current_depth() - 1
        call set_hc(this%dm, this%tr, this%bi, cur, this%nd, this%bs, W)
        !call this%tr%prune(W)
!print'(11f9.3)',W(195:337)
        if (this%tr%finished()) exit
        call this%tr%set_parent_node(W)
        block
          integer(IK) :: cix
          cix = this%tr%current_index()
          this%bi(cur)%isym = (cix - 1) / this%bi(cur)%nper
          call swap_iper(this%nd, cur, cix, this%bi)
!print '(A,6i4,2f9.3)', ' cisn', this%bi%iper, this%bi%isym
        end block
        if (cur == this%nd) then
          block
            integer(IK) :: i, pp
            pp = this%tr%current_pointer()
print '(A,6i4,2f9.3)', ' open', this%bi%iper, this%bi%isym, W(this%tr%upperbound), W(pp)
!
            if (W(this%tr%upperbound) > W(pp)) then
              W(this%tr%upperbound) = W(pp)
              do concurrent(i=1:this%nd)
                this%bi(i)%jper = this%bi(i)%iper
                this%bi(i)%jsym = this%bi(i)%isym
              end do
              print *, 'save', W(this%tr%upperbound)
            end if
          end block
          exit
        end if
      end do
!
      do while (this%tr%finished() .and. cur > 0)
        call this%tr%close_node()
        call paws_iper(this%nd, cur, this%bi)
        !call this%tr%prune(W)
        cur = cur - 1
      end do
!
      if (cur == 0) exit
!
      call this%tr%set_parent_node(W)
      block
        integer(IK) :: cix
        cix = this%tr%current_index()
        this%bi(cur)%isym = (cix - 1) / this%bi(cur)%nper
        call swap_iper(this%nd, cur, cix, this%bi)
!print '(A,6i4,2f9.3)', ' cisn', this%bi%iper, this%bi%isym
      end block
!
    end do
!
    call set_gp(this%dm%l, this%dm%m%g, this%p, this%bi, this%gp)
!
  contains
!
    pure subroutine set_gp(l, g, p, bi, gp)
      integer(IK), intent(in)                :: l, g(l), p(l)
      type(breadth_indicator), intent(in)    :: bi(*)
      type(group_permutation), intent(inout) :: gp(l)
      integer(IK)                            :: i
      do concurrent(i=1:l)
        block
          integer(IK) :: iper(g(i))
          iper = bi(p(i):p(i) + g(i) - 1)%jper
          gp(i) = group_permutation(iper)
        end block
      end do
    end subroutine set_gp
!
    pure subroutine swap_iper(nd, cur, inod, b)
      integer(IK), intent(in)                :: nd, cur, inod
      type(breadth_indicator), intent(inout) :: b(nd)
      integer(IK)                            :: ip, sw, i
!
      ip = MODULO(inod - 1, b(cur)%nper)
      if (ip < 1) return
!
      sw = b(cur + ip)%iper
      do i = cur + ip - 1, cur, -1
        b(i + 1)%iper = b(i)%iper
      end do
      b(cur)%iper = sw
!
    end subroutine swap_iper
!
    pure subroutine paws_iper(nd, cur, bi)
      integer(IK), intent(in)                :: nd, cur
      type(breadth_indicator), intent(inout) :: bi(nd)
      integer(IK)                            :: sw, i
!
      if (cur < 2) return
      do i = cur, nd - 1
        if (bi(i - 1)%ispc /= bi(i)%ispc .or. bi(i - 1)%iper < bi(i)%iper) return
        sw = bi(i - 1)%iper
        bi(i - 1)%iper = bi(i)%iper
        bi(i)%iper = sw
      end do
!
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
          W(t) = ZERO
          call copy(nx, W(ph), W(h))
          call dm%partial_eval(cur, iper, i, j, W, W(t), W(h), W(c))
        end block
      end do
!
    end subroutine set_hc
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
  pure subroutine branch_and_prune_swap(this, X)
    class(branch_and_prune), intent(in) :: this
    real(RK), intent(inout)             :: X(*)
    integer(IK)                         :: i
!
    do i = 1, this%dm%l
      call swap(this%dm%d, this%dm%m(i)%m, this%dm%m(i)%g, &
     &          this%bi(this%p(i)), this%mr(i), this%gp(i), X(this%q(i)))
    end do
!
  contains
!
    pure subroutine swap(d, m, g, bi, mr, gp, X)
      integer(IK), intent(in)              :: d, m, g
      type(breadth_indicator), intent(in)  :: bi(g)
      type(molecular_rotation), intent(in) :: mr
      type(group_permutation), intent(in)  :: gp
      real(RK), intent(inout)              :: X(d, m, g)
      integer(IK)                          :: i
      call gp%swap(d * m, X)
      do concurrent(i=1:g)
        call mr%swap(d, X(1, 1, i), bi(i)%jsym)
      end do
    end subroutine swap
!
  end subroutine branch_and_prune_swap
!
  pure subroutine copy(d, source, dest)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: source(*)
    real(RK), intent(inout) :: dest(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dest(i) = source(i)
    end do
  end subroutine copy
!
  pure elemental subroutine branch_and_prune_clear(this)
    class(branch_and_prune), intent(inout) :: this
    call this%dm%clear()
    call this%tr%clear()
    if (ALLOCATED(this%p)) deallocate (this%p)
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%mr)) deallocate (this%mr)
    if (ALLOCATED(this%gp)) deallocate (this%gp)
  end subroutine branch_and_prune_clear
!
  pure elemental subroutine branch_and_prune_destroy(this)
    type(branch_and_prune), intent(inout) :: this
    call branch_and_prune_clear(this)
  end subroutine branch_and_prune_destroy
!
end module mod_branch_and_prune

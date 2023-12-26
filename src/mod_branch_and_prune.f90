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
    integer(IK)                :: bs
    integer(IK), allocatable   :: p(:), q(:)
    type(d_matrix_list)        :: dmat
    type(tree)                 :: t
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
    integer(IK)                          :: n, i, j, k
!
    res%dmat = d_matrix_list(blk, p)
    n = res%dmat%n_depth()
    allocate (res%bi(n))
!
    k = 0
    do j = 1, res%dmat%l
      do i = res%dmat%m(j)%g, 1, -1
        k = k + 1
        res%bi(k) = breadth_indicator(j, 0, 0, 0, 0, i, res%dmat%m(j)%s, i * res%dmat%m(j)%s)
      end do
    end do
!
    res%bs = res%dmat%dd + 2
    res%t = tree(p + res%dmat%memsize(), res%bs, n + 1, [1, res%bi%nnod])
!
    allocate (res%p(res%dmat%l))
    res%p(1) = 1
    do i=2,res%dmat%l
      res%p(i) = res%p(i - 1) + res%dmat%m(i - 1)%g
    enddo
!
    allocate (res%q(res%dmat%l))
    res%q(1) = 1
    do i=2,res%dmat%l
      res%q(i) = res%q(i - 1) + res%dmat%d * res%dmat%m(i - 1)%m * res%dmat%m(i - 1)%n
    enddo
!
    allocate (res%gp(res%dmat%l))
    allocate (res%mr(res%dmat%l))
    if (PRESENT(rot)) then
      do concurrent(i=1:res%dmat%l)
        res%mr(i) = rot(i)
      end do
    end if
!
  end function branch_and_prune_new
!
  pure elemental function branch_and_prune_memsize(this) result(res)
    class(branch_and_prune), intent(in) :: this
    integer(IK)                         :: res
    res = this%dmat%memsize() + this%t%memsize
  end function branch_and_prune_memsize
!
  pure subroutine branch_and_prune_setup(this, X, Y, W)
    class(branch_and_prune), intent(inout) :: this
    real(RK), intent(in)                   :: X(*)
    real(RK), intent(in)                   :: Y(*)
    real(RK), intent(inout)                :: W(*)
    integer(IK)                            :: p, i, j, k
!
    call this%dmat%eval(this%mr, X, Y, W)
    call this%t%reset()
    call this%t%open_node()
!
    W(this%t%upperbound) = RHUGE
!
    p = this%t%parent_pointer()
    W(p) = W(this%dmat%o)
    p = p + 1
    W(p) = W(this%dmat%h)
    p = p + 1
    call copy(this%dmat%dd, W(this%dmat%c), W(p))
    call this%t%set_parent_node(W)
!
    k = 0
    do j = 1, this%dmat%l
      do i = 1, this%dmat%m(j)%g
        k = k + 1
        this%bi(k)%iper = i
        this%bi(k)%jper = i
        this%bi(k)%isym = 1
        this%bi(k)%jsym = 1
      end do
    end do
!
  end subroutine branch_and_prune_setup
!
   subroutine branch_and_prune_run(this, W)
    class(branch_and_prune), intent(inout) :: this
    real(RK), intent(inout)                :: W(*)
    integer(IK)                            :: pp, cur, nd
!
    nd = this%t%n_depth() - 1
!
    do
!
      cur = this%t%current_depth() - 1
!
      call set_hc(this%dmat, this%t, this%bi, cur, this%bs, W)
      call this%t%prune(W)
      call this%t%set_parent_node(W)
!
      if (cur == nd) then
        pp = this%t%current_pointer()
        if (W(this%t%upperbound) > W(pp)) then
          W(this%t%upperbound) = W(pp)
          block
            integer(IK) :: i
            do concurrent(i=1:nd)
              this%bi(i)%jper = this%bi(i)%iper
              this%bi(i)%jsym = this%bi(i)%isym
            end do
            print*,this%bi%jper
            print*,this%bi%jsym
          end block
        end if
      end if
!
      if (this%t%finished()) then
        if (cur == 1) exit
        call this%t%close_node()
        call paws_iper(nd, cur, this%bi)
      else
        call swap_iper(nd, cur, this%t%current_index(), this%bi)
        call this%t%open_node()
      end if
!
    end do
!
    call set_gp(this%dmat%l, this%dmat%m%g, this%p, this%bi, this%gp)
!
    print*, w(this%dmat%o), W(this%t%upperbound)
    print'(8i4)',this%bi
!
  contains
!
     subroutine set_gp(l, g, p, bi, gp)
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
      b(cur)%isym = (inod - 1) / b(cur)%nper + 1
!
      sw = b(cur + ip)%iper
      do i = cur + ip - 1, cur, -1
        b(i + 1)%iper = b(i)%iper
      enddo
      b(cur)%iper = sw
!
    end subroutine swap_iper
!
    pure subroutine paws_iper(nd, cur, b)
      integer(IK), intent(in)                :: nd, cur
      type(breadth_indicator), intent(inout) :: b(nd)
      integer(IK)                            :: sw, i
!
      sw = b(cur)%iper
      do i = cur + 1, nd
        if (i == nd .or. b(i)%ispc /= b(cur)%ispc .or. sw < b(i)%iper) then
          b(i - 1)%iper = sw
          exit
        end if
        b(i - 1)%iper = b(i)%iper
      enddo
!
    end subroutine paws_iper
!
    subroutine set_hc(dmat, tr, b, cur, bs, W)
      type(d_matrix_list), intent(in)     :: dmat
      type(tree), intent(in)              :: tr
      type(breadth_indicator), intent(in) :: b(:)
      integer(IK), intent(in)             :: cur, bs
      real(RK), intent(inout)             :: W(*)
      integer(IK)                         :: i, j, p, q, ph, pc
!
      p = tr%parent_pointer()
      ph = p + 1
      pc = p + 2
      q = tr%nodes_pointer()
!
      do concurrent(i=1:b(cur)%nper, j=1:b(cur)%nsym)
        block
          integer(IK) :: k, t, h, c
          k = i + b(cur)%nper * (j - 1)
          t = q + bs * (k - 1)
          h = t + 1
          c = t + 2
          w(h) = W(ph)
          call copy(dmat%dd, W(pc), W(c))
          call dmat%partial_eval(cur, b%iper, i, j, W, W(t), W(h), W(c))
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
    res = W(this%t%upperbound)
  end function branch_and_prune_upperbound
!
  subroutine branch_and_prune_swap(this, X)
    class(branch_and_prune), intent(in) :: this
    real(RK), intent(inout)             :: X(*)
    integer(IK)                         :: i
!
print*,this%q/3
    do i = 1, this%dmat%l
      call swap(this%dmat%d, this%dmat%m(i)%m, this%dmat%m(i)%g, &
     &          this%bi(this%p(i)), this%mr(i), this%gp(i), X(this%q(i)))
    end do
!
  contains
!
    subroutine swap(d, m, g, bi, mr, gp, X)
      integer(IK), intent(in)              :: d, m, g
      type(breadth_indicator), intent(in)  :: bi(g)
      type(molecular_rotation), intent(in) :: mr
      type(group_permutation), intent(in)  :: gp
      real(RK), intent(inout)              :: X(d, m, g)
      integer(IK)                          :: i
      call gp%reverse(d * m, X)
      do concurrent(i=1:g)
        call mr%reverse(d, X(1, 1, i), bi(i)%jsym - 1)
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
    call this%dmat%clear()
!   call this%t%clear()
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

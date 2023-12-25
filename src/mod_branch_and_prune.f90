module mod_branch_and_prune
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_mol_block
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
    integer(IK) :: isym
    integer(IK) :: nper
    integer(IK) :: nsym
    integer(IK) :: nnod
  end type breadth_indicator
!
  type branch_and_prune
    integer(IK)                :: bs
    type(d_matrix_list)        :: dmat
    type(tree)                 :: t
    type(breadth_indicator), allocatable  :: b(:)
    type(molecular_rotation), allocatable :: r(:)
  contains
    procedure :: memsize    => branch_and_prune_memsize
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
    allocate (res%b(n))
!
    k = 0
    do j = 1, res%dmat%l
      do i = res%dmat%m(j)%g, 1, -1
        k = k + 1
        res%b(k) = breadth_indicator(j, 0, 0, i, res%dmat%m(j)%s, i * res%dmat%m(j)%s)
      end do
    end do
!
    res%bs = res%dmat%dd + 2
    res%t = tree(p + res%dmat%memsize(), res%bs, n + 1, [1, res%b%nnod])
    call res%t%open_node()
!
    allocate (res%r(res%dmat%l))
    if (PRESENT(rot)) then
      do concurrent(i=1:res%dmat%l)
        res%r(i) = rot(i)
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
    call this%dmat%eval(this%r, X, Y, W)
    call this%t%reset()
!
    p = this%t%parent_pointer()
    W(p) = W(this%dmat%o)
    p = p + 1
    W(p) = W(this%dmat%h)
    p = p + 1
    call copy(this%dmat%dd, W(this%dmat%c), W(p))
!
    W(this%t%upperbound) = RHUGE
!
    k = 0
    do j = 1, this%dmat%l
      do i = 1, this%dmat%m(j)%g
        k = k + 1
        this%b(k)%iper = i
        this%b(k)%isym = 1
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
    print*,this%t%parent_pointer()
    print'(3f9.3)',w(this%t%parent_pointer():this%t%parent_pointer()+1)
    print'(3f9.3)',w(this%t%parent_pointer()+2:this%t%parent_pointer()+10)
    print*, w(this%dmat%o), W(this%t%upperbound)
    print'(6i4)',this%b
    nd = this%t%n_depth()

    do
!
      cur = this%t%current_depth()
!
      if (this%t%finished())then
        call paws_iper(nd, cur, this%b)
        call this%t%close_node()
        if (cur==1) exit
        cycle
      endif
!
      call set_hc(this%dmat, this%t, this%b, cur, W)
      call this%t%prune(W)
      call this%t%set_parent_node(W)
print*,cur, this%t%parent_pointer(), W(this%t%parent_pointer()), W(this%t%upperbound), this%t%parent_index()
!
      call swap_iper(nd, cur, this%t%parent_index(), this%b)
      if (cur==nd) then
        pp = this%t%parent_pointer()
        W(this%t%upperbound) = W(pp)
        print*,pp, W(this%t%upperbound)
        print'(L4,*(i4))',this%t%finished(), cur, nd, this%b%iper
        call paws_iper(nd, cur, this%b)
        call this%t%close_node()
      else
        call this%t%open_node()
      endif
!
    end do
!
  contains
!
    subroutine swap_iper(nd, cur, inod, b)
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
    subroutine paws_iper(nd, cur, b)
      integer(IK), intent(in)                :: nd, cur
      type(breadth_indicator), intent(inout) :: b(nd)
      integer(IK)                            :: sw, i
!
      sw = b(cur)%iper
      do i = cur + 1, nd
        if (b(i)%ispc /= b(cur)%ispc.or.sw < b(i)%iper) then
          b(i - 1)%iper = sw
          exit
        end if
        b(i - 1)%iper = b(i)%iper
      enddo
!
    end subroutine paws_iper
!
    subroutine set_hc(dmat, tr, b, cur, W)
      type(d_matrix_list), intent(in)     :: dmat
      type(tree), intent(in)              :: tr
      type(breadth_indicator), intent(in) :: b(:)
      integer(IK), intent(in)             :: cur
      real(RK), intent(inout)             :: W(*)
      integer(IK)                         :: i, j, p, q, ph, pc, bs
!
      bs = dmat%dd + 2
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
! recursive subroutine breadth_search(idep, ndep, dmat, W, childs, upperbound, bstper)
!   integer(IK), intent(in)         :: idep, ndep
!   type(d_matrix_list), intent(in) :: dmat
!   real(RK), intent(in)            :: W(*)
!   type(breadth), intent(inout)    :: childs(*)
!   real(RK), intent(inout)         :: upperbound
!   integer(IK), intent(inout)      :: bstper(*)
!   integer(IK)                     :: i
!
!   if (idep == ndep) then
!     block
!       real(RK) :: lb
!       call childs(idep)%set_node_minloc()
!       lb = childs(idep)%lowerbound()
!       if (upperbound > lb) then
!         upperbound = lb
!         do i = 1, ndep
!           print *, childs(i)%iper(), childs(i)%isym()
!         end do
!       end if
!     end block
!     return
!   end if
!
!   do while (childs(idep)%not_finished())
!     call childs(idep)%set_node_index()
!     childs(idep + 1) = childs(idep)%generate_breadth(dmat, W)
!     call breadth_search(idep + 1, ndep, dmat, W, childs, upperbound, bstper)
!     call childs(idep)%prune(upperbound)
!     print*,childs(idep)%nodes%alive
!   end do
!
! end subroutine breadth_search
!
  pure elemental subroutine branch_and_prune_clear(this)
    class(branch_and_prune), intent(inout) :: this
    call this%dmat%clear()
!   call this%t%clear()
    if (ALLOCATED(this%r)) deallocate (this%r)
  end subroutine branch_and_prune_clear
!
  pure elemental subroutine branch_and_prune_destroy(this)
    type(branch_and_prune), intent(inout) :: this
    call branch_and_prune_clear(this)
  end subroutine branch_and_prune_destroy
!
end module mod_branch_and_prune

module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_molecular_rotation
  use mod_molecular_permutation
  use mod_mol_block
  use mod_block_lower_bound
  implicit none
  private
  public :: node, breadth
!
  type node
    private
    type(mol_block_list)                     :: blk
    !! molecular permutation block
    type(molecular_permutation), allocatable :: per(:)
    !! molecular permutation
    real(RK), public                         :: lowerbound = RHUGE
    !! the lower bound
  contains
    procedure         :: generate_breadth => node_generate_breadth
    final             :: node_destroy
  end type node
!
  type breadth
    type(node), allocatable :: nodes(:)
  contains
    final             :: breadth_destroy
  end type breadth
!
  interface node
    module procedure node_new
  end interface node
!
contains
!
!| generate node instance
  pure function node_new(blk, rot, per, x, y) result(res)
    class(mol_block_list), intent(in)        :: blk
    class(molecular_rotation), intent(in)    :: rot(:)
    !! molecular rotation
    class(molecular_permutation), intent(in) :: per(:)
    !! molecular permutation
    real(RK), intent(in)                     :: x(*)
    !! reference molecular coordinate, x(d,*)
    real(RK), intent(in)                     :: y(*)
    !! target molecular coordinate, y(d,*)
    real(RK)                                 :: x_(blk%d * blk%n_atom())
    real(RK)                                 :: y_(blk%d * blk%n_atom())
    real(RK)                                 :: w(block_lower_bound_worksize(blk))
    integer(IK)                              :: i, s, mn, dmn
    type(node)                               :: res
!
    s = blk%n_spc()
    mn = blk%n_atom()
    dmn = blk%d * mn
!
    res%blk = blk
    allocate(res%per(s))
    do concurrent(i=1:s)
      res%per(i) = per(i)
    end do
!
    call copy(dmn, x, x_)
    call centering(blk%d, mn, x_)
    call copy(dmn, y, y_)
    do concurrent(i=1:s)
      call per(i)%mol_swap(blk%d, blk%b(i)%m, rot(i), y_(blk%b(i)%p))
    end do
    call centering(blk%d, mn, y_)
    call block_lower_bound(blk, x_, y_, w)
    res%lowerbound = w(1)
!
  end function node_new
!
  pure subroutine copy(n, x, w)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: x(*)
    real(RK), intent(inout) :: w(*)
    integer(IK)             :: i
    do concurrent(i=1:n)
      w(i) = x(i)
    end do
  end subroutine copy
!
  pure subroutine centering(d, n, x)
    integer(IK), intent(in) :: d, n
    real(RK), intent(inout) :: x(d, n)
    real(RK)                :: c(d)
    integer(IK)             :: i, j
    c = SUM(X, 2) / n
    do concurrent(i=1:d, j=1:n)
      x(i, j) = x(i, j) - c(i)
    end do
  end subroutine centering
!
!| generate childe nodes instance
  pure function node_generate_breadth(this, rot, x, y) result(res)
    class(node), intent(in)               :: this
    !! molecular template
    class(molecular_rotation), intent(in) :: rot(:)
    !! molecular template
    real(RK), intent(in)                  :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)                  :: y(*)
    !! target molecular coordinate, y(d,m,n)
    type(breadth)                         :: res
    type(mol_block_list)                  :: b_child
    integer(IK)                           :: nper, nsym, nspc
    integer(IK)                           :: iper, isym, ispc
!
    b_child = this%blk%child()
    ispc = b_child%ispecies()
    nspc = b_child%nspecies()
    if (ispc < 1) ispc = nspc
    nper = b_child%b(ispc)%g + 1
    nsym = rot(ispc)%n_sym()
!
    allocate (res%nodes(nper * (nsym + 1)))
!
    do concurrent(iper=1:nper, isym=0:nsym)
      block
        type(molecular_permutation) :: p_child(nspc)
        integer(IK)                 :: i
        do concurrent(i=1:ispc - 1)
          p_child(i) = this%per(i)
        end do
        p_child(ispc) = this%per(ispc)%child(iper, isym)
        do concurrent(i=ispc + 1:nspc)
          p_child(i) = this%per(i)
        end do
        res%nodes(isym * nper + iper) = node(b_child, rot, p_child, x, y)
      end block
    end do
!
  end function node_generate_breadth
!
  pure elemental subroutine node_destroy(this)
    type(node), intent(inout) :: this
    this%lowerbound = RHUGE
  end subroutine node_destroy
!
  pure elemental subroutine breadth_destroy(this)
    type(breadth), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
  end subroutine breadth_destroy
!
end module mod_tree

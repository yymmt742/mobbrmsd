module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_molecule
  use mod_molecular_permutation
  use mod_lower_bound
  implicit none
  private
  public :: node, childs
!
  type node
!
    private
!
    type(molecular_permutation) :: prm
    !! Molecular permutation indicator
    real(RK), public            :: lower = -HUGE(ZERO)
    !! the lower bound
!
  contains
!
    procedure         :: generate_childs => node_generate_childs
    final             :: node_destroy
!
  end type node
!
  type childs
    type(node), allocatable :: nodes(:)
  contains
    final             :: childs_destroy
  end type childs
!
  interface node
    module procedure node_new
  end interface node
!
contains
!
!| generate node instance
  pure function node_new(mol, prm, x, y) result(res)
    class(molecule), intent(in)              :: mol
    !! molecular template
    class(molecular_permutation), intent(in) :: prm
    !! molecular permutation
    real(RK), intent(in)                     :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)                     :: y(*)
    !! target molecular coordinate, y(d,m,n)
    real(RK)                                 :: w(node_worksize(mol, prm))
    type(node)                               :: res
    integer(IK)                              :: d, m, n, dm, dmn, ix, iy, iw
!
    res%prm = prm
!
    d = prm%d
    m = mol%natom()
    n = prm%n
    dm  = d * m
    dmn = dm * n
!
    ix = 1
    iy = ix + dmn
    iw = iy + dmn
!
    w(ix:ix + dmn - 1) = x(:dmn)
    call prm%sort_matrix(mol, y, w(iy))
!
    call lower_bound(dm, n, prm%free_indices(), w(ix), w(iy), w(iw))
    res%lower = w(iw)
!
  end function node_new
!
!| generate childe nodes instance
  pure function node_generate_childs(this, mol, x, y) result(res)
    class(node), intent(in) :: this
    class(molecule), intent(in)              :: mol
    !! molecular template
    real(RK), intent(in)    :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)    :: y(*)
    !! target molecular coordinate, y(d,m,n)
    type(childs)            :: res
    integer(IK)             :: f, g, i
!
!   f = SIZE(this%free_indices)
!   g = SIZE(this%fixed_indices)
!
!   ALLOCATE(res%nodes(f))
!
!   do concurrent(i=1:f)
!     block
!       integer(IK) :: fixed(g + 1), free(f - 1)
!       fixed = [this%fixed_indices, this%free_indices(i)]
!       free = [this%free_indices(:i - 1), this%free_indices(i + 1:)]
!       res%nodes(i) = node(d, m, n, fixed, free, x, y)
!     end block
!   end do
!
  end function node_generate_childs
!
  pure function node_worksize(mol, prm) result(res)
    class(molecule), intent(in)              :: mol
    class(molecular_permutation), intent(in) :: prm
    integer(IK)                              :: d, m, n, res
!
    res = 2 * (prm%d * mol%natom() * prm%n) &
   &    + lower_bound_worksize(prm%d * mol%natom(), prm%n, prm%free_indices())
!
  end function node_worksize
!
! pure function nlist(m, free) result(res)
! integer(IK), intent(in) :: m, free(:)
! integer(IK)             :: res(SIZE(free) * m)
! integer(IK)             :: i, j
!
!   do concurrent(j=1:SIZE(free), i=1:m)
!     res((j - 1) * m + i) = (free(j) - 1) * m + i
!   end do
!
! end function nlist
!
  pure elemental subroutine node_destroy(this)
    type(node), intent(inout) :: this
    this%lower = -HUGE(ZERO)
  end subroutine node_destroy
!
  pure elemental subroutine childs_destroy(this)
    type(childs), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
  end subroutine childs_destroy
!
end module mod_tree

module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_molecule
  use mod_molecular_rotation
  use mod_molecular_permutation
  use mod_block_lower_bound
  implicit none
  private
  public :: node, childs
!
  type node
    sequence
    private
    integer(IK)              :: depth = 0
    !! node depth
    integer(IK)              :: iprm  = 0
    !! molecular permutation indicator
    integer(IK)              :: isym  = 0
    !! molecular symmetry indicator
    real(RK), public         :: lower = -HUGE(ZERO)
    !! the lower bound
  contains
    procedure         :: generate_childs => node_generate_childs
    final             :: node_destroy
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
  pure function node_new(rot, d, m, n, k, iper, isym, x, yfix, yfree) result(res)
    class(molecular_rotation), intent(in) :: rot
    !! molecular rotation
    integer(IK), intent(in)               :: d
    !! spatial dimension
    integer(IK), intent(in)               :: m
    !! molecular rotation
    integer(IK), intent(in)               :: n
    !! molecular rotation
    integer(IK), intent(in)               :: k
    !! node depth
    integer(IK), intent(in)               :: iper
    !! permutation index
    integer(IK), intent(in)               :: isym
    !! rotation index
    real(RK), intent(in)                  :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)                  :: yfix(*)
    !! target molecular coordinate, y(d,m,depth)
    real(RK), intent(in)                  :: yfree(*)
    !! target molecular coordinate, y(d,m,n-depth)
    real(RK)                              :: w(node_worksize(mol, prm))
    type(node)                            :: res
    integer(IK)                           :: i, d, m, n, g, mn, dmn
!
    dmk = d * m * k
    do concurrent(i=1:dmk)
      w(i) = yfix(i)
    end do
    call block_lower_bound(d, m, n, mol%free_index(), [(i, i=g + 1, n)], x, w, w(dmn + 1))
    res%lower = w(dmn + 1)
!
  end function node_new
!
  pure subroutine sort_matrix(rot, d, m, n, depth, iprm, isym, yfix, yfree, w)
    class(molecular_rotation), intent(in) :: rot
    !! molecular rotation indicator
    integer(IK), intent(in)               :: d, n, m, depth, iprm, isym
    !! molecular rotation
    real(RK), intent(in)                  :: yfix(d, m, depth)
    !! target molecular coordinate, y(d,m,depth)
    real(RK), intent(in)                  :: yfree(d, m, n - depth)
    !! target molecular coordinate, y(d,m,n-depth)
    real(RK), intent(inout)               :: w(*)
    integer(IK)                           :: s
!
    if (ALLOCATED(prm%sym)) then
      call sort(mol, prm%d, mol%natom(), prm%n, prm%nfree(), prm%sym, prm%swap_indices(), source, dest)
    else
      call sort(mol, prm%d, mol%natom(), prm%n, prm%nfree(), [(1, s=1, prm%n)], prm%swap_indices(), source, dest)
    end if
!
  contains
!
    pure subroutine sort(mol, d, m, n, f, sym, swp, source, dest)
      class(molecule), intent(in) :: mol
      integer(IK), intent(in)     :: d, m, n, f, sym(f), swp(n)
      real(RK), intent(in)        :: source(d, m, n)
      real(RK), intent(inout)     :: dest(d, m, n)
      integer(IK)                 :: i, j, k
!
      do concurrent(k=1:f)
        block
          integer(IK) :: w(m)
          w = mol%sym_index(sym(k))
          do concurrent(i=1:d, j=1:m)
            dest(i, j, k) = source(i, w(j), swp(k))
          end do
        end block
      end do
!
      do concurrent(i=1:d, j=1:m, k=f + 1:n)
        dest(i, j, k) = source(i, j, swp(k))
      end do
!
    end subroutine sort
!
  end subroutine sort_matrix
!
!| generate childe nodes instance
  pure function node_generate_childs(this, mol, x, y) result(res)
    class(node), intent(in)     :: this
    class(molecule), intent(in) :: mol
    !! molecular template
    real(RK), intent(in)        :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)        :: y(*)
    !! target molecular coordinate, y(d,m,n)
    type(childs)                :: res
    integer(IK)                 :: f, g, i
!
    f = this%prm%nfree()
    g = this%prm%nfix()
!
    ALLOCATE(res%nodes(f))
!
    do concurrent(i=1:f)
      block
        type(molecular_permutation) :: prm
        integer(IK)                 :: fix(g + 1)
        fix = [this%prm%fix, this%prm%free(i)]
        prm = molecular_permutation(mol, f + g, fix=fix, d=prm%d)
        res%nodes(i) = node(mol, prm, x, y)
      end block
    end do
!
  end function node_generate_childs
!
  pure function node_worksize(mol, prm) result(res)
    class(molecule), intent(in)              :: mol
    class(molecular_permutation), intent(in) :: prm
    integer(IK)                              :: n, m, g, i, res
!
    n = prm%n
    m = mol%natom()
    g = prm%nfix()
    res = prm%d * m * n &
   &    + block_lower_bound_worksize(prm%d, m, n, mol%free_index(), [(i, i=g + 1, n)])
!
  end function node_worksize
!
  pure function free_indices(mol, prm) result(res)
    class(molecule), intent(in)              :: mol
    class(molecular_permutation), intent(in) :: prm
    integer(IK)                              :: res(mol%natom() * prm%nfree())
    integer(IK)                              :: n, m, f, i, j, nfix
!
    n = prm%n
    m = mol%natom()
    f = prm%nfree()
    nfix = m * (n - f)
!
    do concurrent(i=1:m, j=1:f)
      res(i + (j - 1) * m) =  nfix + i + m * (j - 1)
    end do
!
  end function free_indices
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

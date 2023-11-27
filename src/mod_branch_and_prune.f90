module mod_branch_and_prune
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_molecule
  use mod_molecular_permutation
  use mod_tree
  implicit none
  private
  public :: node, childs
!
contains
!
!| generate node instance
  pure subroutine branch_and_prune(mol, x, y) result(res)
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
    integer(IK)                              :: d, m, n, mn, dmn
!
    res%prm = prm
!
    d = prm%d
    m = mol%natom()
    n = prm%n
    mn  = m * n
    dmn = d * mn
!
    call sort_matrix(mol, res%prm, y, w)
    call lower_bound(d, mn, free_indices(mol, prm), x, w, w(dmn + 1))
    res%lower = w(dmn + 1)
!
  end function branch_and_prune
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
    integer(IK)                              :: nm, i, res
!
    nm = prm%n * mol%natom()
    res = prm%d * mol%natom() * prm%n
    if (ALLOCATED(prm%free)) then
      res = lower_bound_worksize(prm%d, nm, free_indices(mol, prm))
    else
      res = lower_bound_worksize(prm%d, nm, [(i, i=1, nm)])
    end if
!
  end function node_worksize
!
  pure subroutine sort_matrix(mol, prm, source, dest)
    class(molecular_permutation), intent(in) :: prm
    class(molecule), intent(in)              :: mol
    real(RK), intent(in)                     :: source(*)
    real(RK), intent(inout)                  :: dest(*)
    integer(IK)                              :: s
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
          w = mol%sym_index(sym(swp(k)))
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
  pure elemental subroutine branch_and_prune_destroy(this)
    type(childs), intent(inout) :: this
    call this%mol%clear()
    call this%prm%clear()
  end subroutine branch_and_prune_destroy
!
end module mod_branch_and_prune

program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_molecule
  use mod_molecular_permutation
  use mod_tree
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test node')
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter          :: d = 3
    integer, parameter          :: m = 12
    integer, parameter          :: n = 8
    integer, parameter          :: sym(m) = [1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8]
    real(RK)                    :: X(d * m * n), Y(d * m * n)
    type(molecule)              :: mol
    type(molecular_permutation) :: prm
    type(node)                  :: w
    type(childs)                :: c, p, q, r
!
    mol = molecule(m, d=1, sym=sym)
    prm = molecular_permutation(mol, n, sym=[1, 1, 1, 1], fix=[1, 2, 3, 4])
!   print*,prm%swap_indices()
!
    X = [sample(d, m * n)]
    Y = X + 0.00*[sample(d, m * n)]
!
    w = node(mol, prm, X, Y)
    print *, w%lower
    c = w%generate_childs(mol, X, Y)
    print *, c%nodes(:)%lower
    p = c%nodes(1)%generate_childs(mol, X, Y)
    print *, p%nodes(:)%lower
    q = p%nodes(1)%generate_childs(mol, X, Y)
    print *, q%nodes(:)%lower
    r = q%nodes(1)%generate_childs(mol, X, Y)
    print *, r%nodes(:)%lower
!
    p = c%nodes(2)%generate_childs(mol, X, Y)
    print *, p%nodes(:)%lower
    p = c%nodes(3)%generate_childs(mol, X, Y)
    print *, p%nodes(:)%lower
    p = c%nodes(4)%generate_childs(mol, X, Y)
    print *, p%nodes(:)%lower
    print*
!
    mol = molecule(m, d=1, sym=[1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8])
    prm = molecular_permutation(mol, n, sym=[1, 1, 1, 1], fix=[5, 6, 7, 8])
!
    w = node(mol, prm, X, Y)
    print *, w%lower
    c = w%generate_childs(mol, X, Y)
    print *, c%nodes(:)%lower
    p = c%nodes(1)%generate_childs(mol, X, Y)
    print *, p%nodes(:)%lower
    p = c%nodes(2)%generate_childs(mol, X, Y)
    print *, p%nodes(:)%lower
    p = c%nodes(3)%generate_childs(mol, X, Y)
    print *, p%nodes(:)%lower
    p = c%nodes(4)%generate_childs(mol, X, Y)
    print *, p%nodes(:)%lower
!
  end subroutine test1
!
  function sample(d, n) result(res)
    integer, intent(in)  :: d, n
    real(RK)             :: cnt(d)
    real(RK)             :: res(d, n)
    integer              :: i
    call RANDOM_NUMBER(res)
    cnt = SUM(res, 2) / n
    do concurrent(i=1:n)
      res(:, i) = res(:, i) - cnt
    enddo
  end function sample
!
end program main

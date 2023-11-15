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
    real(RK)                    :: X(d * m * n), Y(d * m * n)
    type(molecule)              :: mol
    type(molecular_permutation) :: prm
    type(node)                  :: w
    type(childs)                :: c, p
!
    mol = molecule(m, RESHAPE([1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8], [12, 1]))
    prm = molecular_permutation(mol, n, sym=[1, 2, 1, 1, 2, 2, 1, 2], fix=[1, 2, 3, 4])
!
    X = [sample(d, m * n)]
    Y = X + 0.5*[sample(d, m * n)]
!
    w = node(mol, prm, X, Y)
    print *, w%lower
!   c = w%generate_childs(d, m, n, X, Y)
!   print *, c%nodes(:)%lower
!   p = c%nodes(1)%generate_childs(d, m, n, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(2)%generate_childs(d, m, n, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(3)%generate_childs(d, m, n, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(4)%generate_childs(d, m, n, X, Y)
!   print *, p%nodes(:)%lower
!
!   w = node(d, m, n, [7, 4, 3, 1], [2, 5, 6, 8], X, Y)
!   print *, w%lower
!   c = w%generate_childs(d, m, n, X, Y)
!   print *, c%nodes(:)%lower
!   p = c%nodes(1)%generate_childs(d, m, n, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(2)%generate_childs(d, m, n, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(3)%generate_childs(d, m, n, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(4)%generate_childs(d, m, n, X, Y)
!   print *, p%nodes(:)%lower
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

program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_block_lower_bound
  use mod_molecular_rotation
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
    integer, parameter       :: d = 3
    integer, parameter       :: s = 3
    integer, parameter       :: m1 = 5, n1 = 3, f1 = 3, g1 = 0
    integer, parameter       :: m2 = 3, n2 = 4, f2 = 2, g2 = 3
    integer, parameter       :: m3 = 7, n3 = 5, f3 = 7, g3 = 5
    integer, parameter       :: mn = m1 * n1 + m2 * n2 + m3 * n3
    type(molecular_rotation) :: r(s)
    type(mol_block_list)     :: b
    real(RK)                 :: X(d, mn), Y(d, mn)
    type(node)               :: a
    integer                  :: i
!
    b = mol_block_list(d, s, [m1,m2,m3], [n1,n2,n3], [f1,f2,f3])
    b%b(1)%g = g1
    b%b(2)%g = g2
    b%b(3)%g = g3
    r(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m1, 2]))
    r(2) = molecular_rotation(RESHAPE([2, 1, 3], [m2, 1]))
    r(3) = molecular_rotation(RESHAPE([7, 6, 5, 4, 3, 2, 1], [m3, 1]))
!
    do concurrent(i=1:mn)
      X(:, i) = i
      Y(:, i) = i
    enddo
!
!   a = node(b, x, y)
!
  end subroutine test1
!
! subroutine test1()
!   integer, parameter          :: d = 3
!   integer, parameter          :: m = 12
!   integer, parameter          :: n = 8
!   integer, parameter          :: sym(m) = [1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8]
!   real(RK)                    :: X(d * m * n), Y(d * m * n)
!   type(molecule)              :: mol
!   type(molecular_permutation) :: prm
!   type(node)                  :: w
!   type(childs)                :: c, p, q, r
!
!   mol = molecule(m, d=1, sym=sym)
!   prm = molecular_permutation(mol, n, sym=[1, 1, 1, 1], fix=[1, 2, 3, 4])
!   print*,prm%swap_indices()
!
!   X = [sample(d, m * n)]
!   Y = X + 0.00*[sample(d, m * n)]
!
!   w = node(mol, prm, X, Y)
!   print *, w%lower
!   c = w%generate_childs(mol, X, Y)
!   print *, c%nodes(:)%lower
!   p = c%nodes(1)%generate_childs(mol, X, Y)
!   print *, p%nodes(:)%lower
!   q = p%nodes(1)%generate_childs(mol, X, Y)
!   print *, q%nodes(:)%lower
!   r = q%nodes(1)%generate_childs(mol, X, Y)
!   print *, r%nodes(:)%lower
!
!   p = c%nodes(2)%generate_childs(mol, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(3)%generate_childs(mol, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(4)%generate_childs(mol, X, Y)
!   print *, p%nodes(:)%lower
!   print*
!
!   mol = molecule(m, d=1, sym=[1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8])
!   prm = molecular_permutation(mol, n, sym=[1, 1, 1, 1], fix=[5, 6, 7, 8])
!
!   w = node(mol, prm, X, Y)
!   print *, w%lower
!   c = w%generate_childs(mol, X, Y)
!   print *, c%nodes(:)%lower
!   p = c%nodes(1)%generate_childs(mol, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(2)%generate_childs(mol, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(3)%generate_childs(mol, X, Y)
!   print *, p%nodes(:)%lower
!   p = c%nodes(4)%generate_childs(mol, X, Y)
!   print *, p%nodes(:)%lower
!
! end subroutine test1
!
! function sample(d, n) result(res)
!   integer, intent(in)  :: d, n
!   real(RK)             :: cnt(d)
!   real(RK)             :: res(d, n)
!   integer              :: i
!   call RANDOM_NUMBER(res)
!   cnt = SUM(res, 2) / n
!   do concurrent(i=1:n)
!     res(:, i) = res(:, i) - cnt
!   enddo
! end function sample
!
end program main

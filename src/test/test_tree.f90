program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_block_lower_bound
  use mod_mol_block
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
    integer, parameter          :: d = 3
    integer, parameter          :: s = 3
    integer, parameter          :: m1 = 5, n1 = 3, f1 = 3
    integer, parameter          :: m2 = 3, n2 = 4, f2 = 2
    integer, parameter          :: m3 = 7, n3 = 5, f3 = 7
    integer, parameter          :: mn = m1 * n1 + m2 * n2 + m3 * n3
    type(mol_block_list)        :: blk
    type(molecular_rotation)    :: rot(s)
    type(molecular_permutation) :: per(s)
    real(RK)                    :: X(d, mn), Y(d, mn)
    type(node)                  :: a
    type(breadth)               :: z
    integer                     :: i
!
    blk = mol_block_list(d, s, [m1,m2,m3], [n1,n2,n3], [f1,f2,f3])
    rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m1, 2]))
    rot(2) = molecular_rotation(RESHAPE([2, 1, 3], [m2, 1]))
    rot(3) = molecular_rotation(RESHAPE([7, 6, 5, 4, 3, 2, 1], [m3, 1]))
    per(1) = molecular_permutation(n1)
    per(2) = molecular_permutation(n2)
    per(3) = molecular_permutation(n3)
!
    do concurrent(i=1:mn)
      X(:, i) = i
      Y(:, i) = i
    enddo
!
    a = node(blk, rot, per, [x], [y])
    z = a%generate_breadth(rot, [x], [y])
    print'(*(f9.3))',z%nodes%lower
!
  end subroutine test1
!
! subroutine test2()
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
! end subroutine test2
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

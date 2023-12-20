program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_Procrustes
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
! call test2()
! call test3()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter          :: d = 3
    integer, parameter          :: s = 3
    integer, parameter          :: m1 = 5, n1 = 3, f1 = 3, g1 = 2
    integer, parameter          :: m2 = 3, n2 = 4, f2 = 2, g2 = 2
    integer, parameter          :: m3 = 7, n3 = 5, f3 = 7, g3 = 3
    integer, parameter          :: mn = m1 * n1 + m2 * n2 + m3 * n3
    type(mol_block_list)        :: blk
    type(molecular_rotation)    :: rot(s)
    type(molecular_permutation) :: per(s)
    real(RK)                    :: X(d, mn), Y(d, mn)
    type(node)                  :: a
    type(breadth)               :: z1, z2, z3, z4, z5, z6, z7, z8
    integer                     :: i, ml
!
    blk = mol_block_list(d, s, [m1,m2,m3], [n1,n2,n3], [f1,f2,f3], [g1,g2,g3])
    rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m1, 2]))
    rot(2) = molecular_rotation(RESHAPE([(i, i=1,0)], [0, 1]))
    rot(3) = molecular_rotation(RESHAPE([7, 6, 5, 4, 3, 2, 1], [m3, 1]))
    per(1) = molecular_permutation(n1)
    per(2) = molecular_permutation(n2)
    per(3) = molecular_permutation(n3)
!
    X = sample(d,mn)
    Y = 0.5D0 * X + 0.5D0 * sample(d, mn)
!
    a = node(blk, [x], [y])
    print'(*(f16.3))',a%lowerbound
    z1 = a%generate_breadth(rot, [x], [y])
    print'(*(f16.3))',z1%nodes%lowerbound
    ml = MINLOC(z1%nodes%lowerbound, 1)
    z2 = z1%nodes(ml)%generate_breadth(rot, [x], [y])
    print'(*(f16.3))',z2%nodes%lowerbound
    ml = MINLOC(z2%nodes%lowerbound, 1)
    z3 = z2%nodes(ml)%generate_breadth(rot, [x], [y])
    print'(*(f16.3))',z3%nodes%lowerbound
    ml = MINLOC(z3%nodes%lowerbound, 1)
    z4 = z3%nodes(ml)%generate_breadth(rot, [x], [y])
    print'(*(f16.3))',z4%nodes%lowerbound
!
    ml = MINLOC(z4%nodes%lowerbound, 1)
    z5 = z4%nodes(ml)%generate_breadth(rot, [x], [y])
    print'(*(f16.3))',z5%nodes%lowerbound
!
    ml = MINLOC(z5%nodes%lowerbound, 1)
    z6 = z5%nodes(ml)%generate_breadth(rot, [x], [y])
    print'(*(f16.3))',z6%nodes%lowerbound
!
    ml = MINLOC(z6%nodes%lowerbound, 1)
    z7 = z6%nodes(ml)%generate_breadth(rot, [x], [y])
    print'(*(f16.3))',z7%nodes%lowerbound
!
    ml = MINLOC(z7%nodes%lowerbound, 1)
    z8 = z7%nodes(ml)%generate_breadth(rot, [x], [y])
    print'(*(f16.3))',z8%nodes%lowerbound
!
  end subroutine test1
!
! subroutine test2()
!   integer, parameter          :: d = 3
!   integer, parameter          :: s = 2
!   integer, parameter          :: m1 = 5, n1 = 3, f1 = 3
!   integer, parameter          :: m2 = 3, n2 = 4, f2 = 2
!   integer, parameter          :: mn = m1 * n1 + m2 * n2
!   type(mol_block_list)        :: blk
!   type(molecular_rotation)    :: rot(s)
!   type(molecular_permutation) :: per(s)
!   real(RK)                    :: X(d, mn), Y(d, mn)
!   type(node)                  :: a
!   type(breadth)               :: z(6)
!   integer                     :: i, m
!
!   blk = mol_block_list(d, s, [m1,m2], [n1,n2], [f1,f2])
!   rot(1) = molecular_rotation(RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m1, 2]))
!   rot(2) = molecular_rotation(RESHAPE([2, 1, 3], [m2, 1]))
!   per = molecular_permutation([n1, n2])
!
!   do concurrent(i=1:mn)
!     X(:, i) = [1, i, i*i]
!   enddo
!   Y = X
!
!   a = node(blk, rot, per, [x], [y])
!   z(1) = a%generate_breadth(rot, [x], [y])
!   m    = MINLOC(z(1)%nodes%lowerbound, 1)
!   print'(I4,*(f9.2))',m, z(1)%nodes%lowerbound
!   z(2) = z(1)%nodes(m)%generate_breadth(rot, [x], [y])
!   m    = MINLOC(z(2)%nodes%lowerbound, 1)
!   print'(I4,*(f9.2))',m, z(2)%nodes%lowerbound
!   z(3) = z(2)%nodes(m)%generate_breadth(rot, [x], [y])
!   m    = MINLOC(z(3)%nodes%lowerbound, 1)
!   print'(I4,*(f9.2))',m, z(3)%nodes%lowerbound
!   z(4) = z(3)%nodes(m)%generate_breadth(rot, [x], [y])
!   m    = MINLOC(z(4)%nodes%lowerbound, 1)
!   print'(I4,*(f9.2))',m, z(4)%nodes%lowerbound
!   z(5) = z(4)%nodes(m)%generate_breadth(rot, [x], [y])
!   m    = MINLOC(z(5)%nodes%lowerbound, 1)
!   print'(I4,*(f9.2))',m, z(5)%nodes%lowerbound
!   z(6) = z(5)%nodes(m)%generate_breadth(rot, [x], [y])
!   m    = MINLOC(z(6)%nodes%lowerbound, 1)
!   print'(I4,*(f9.2))',m, z(6)%nodes%lowerbound
!
! end subroutine test2
!
! subroutine test3()
!   integer, parameter          :: d = 3
!   integer, parameter          :: s = 1
!   integer, parameter          :: m = 1, n = 6, f = 0
!   integer, parameter          :: swap(n) = [2,1,3,4,5,6]
!   !integer, parameter          :: swap(n) = [1,2,3,4,5,6]
!   !integer, parameter          :: swap(n) = [2,5,1,6,3,4]
!   integer, parameter          :: flip(m, 2) = RESHAPE([1, 1], [m, 2])
!   !integer, parameter          :: flip(m, 2) = RESHAPE([2, 3, 1, 4, 5, 3, 1, 2, 4, 5], [m, 2])
!   integer, parameter          :: iflp(n) = [0,0,0,0,0,0]
!   !integer, parameter          :: iflp(n) = [2,0,1,2,0,1]
!   !integer, parameter          :: swap(n) = [5,9,2,7,6,4,1,8,3,10]
!   type(mol_block_list)        :: blk
!   type(molecular_rotation)    :: rot(1)
!   type(molecular_permutation) :: per(1)
!   real(RK)                    :: X(d, m, n), Y(d, m, n)
!   real(RK)                    :: P(m*n,m*n)
!   real(RK), allocatable       :: w(:)
!   type(node)                  :: a
!   type(breadth)               :: z(9)
!   integer                     :: i, j, ml
!
!   blk = mol_block_list(d, s, [m], [n], [f])
!   rot(1) = molecular_rotation(flip)
!   per = molecular_permutation([n])
!
!all random_number(X)
!   do concurrent(i=1:m, j=1:n)
!     X(:, i, j) = [1, i, 10*j]
!   enddo
!   do concurrent(i=1:m, j=1:n)
!     if (iflp(j) < 1) then
!       Y(:, i, swap(j)) = X(:, i, j)
!     else
!       Y(:, flip(i, iflp(j)), swap(j)) = X(:, i, j)
!     end if
!   enddo
!
!llocate(w(block_lower_bound_worksize(blk)))
!rint'(3f9.3)', X
!rint*
!rint'(3f9.3)', Y
!rint*
!all Procrustes(6, MATMUL(TRANSPOSE(X(:, 1, :)), Y(:, 1, :)), P, w)
!rint'(6f9.3)', X(:, 1, :) - MATMUL(Y(:, 1, :), TRANSPOSE(P))
!rint'(6f9.3)', P
!rint *
!eturn
!   !blk%b(1)%g = 4
!   print*, blk%d, blk%b
!   call block_lower_bound(blk, x, y, w, nrand=0)
!   print*,w(1)
!   print'(6f9.3)', MATMUL( TRANSPOSE(reshape(X,[3,6])),reshape(Y,[3,6]))
!rint*
!   a = node(blk, rot, per, [x], [y])
!   print*, a%lowerbound
!   return
!   do i = 1, 9
!     z(i) = a%generate_breadth(rot, [x], [y])
!     ml   = MINLOC(z(i)%nodes%lowerbound, 1)
!     print'(I4,*(f9.2))',ml, z(i)%nodes%lowerbound
!   end do
!
! end subroutine test3
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

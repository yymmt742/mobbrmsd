module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_molecular_rotation
  use mod_molecular_permutation
  use mod_block_lower_bound
  implicit none
  private
  public :: node, breadth
!
  type node
    private
    real(RK), public         :: lower = -HUGE(ZERO)
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
  pure function node_new(b, x, y) result(res)
    type(mol_block_list), intent(in)      :: b
    !! molecular block
    real(RK), intent(in)                  :: x(*)
    !! reference molecular coordinate, x(d,*)
    real(RK), intent(in)                  :: y(*)
    !! target molecular coordinate, y(d,*)
    real(RK)                              :: w(block_lower_bound_worksize(b))
    type(node)                            :: res
!
    call block_lower_bound(b, x, y, w)
    res%lower = w(1)
!
  end function node_new
!
!| generate childe nodes instance
  pure function node_generate_breadth(this, rot, b, x, y) result(res)
    class(node), intent(in)               :: this
    !! molecular template
    class(molecular_rotation), intent(in) :: rot(:)
    !! molecular template
    type(mol_block_list), intent(in)      :: b
    !! mol_block_list
    real(RK), intent(in)                  :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)                  :: y(*)
    !! target molecular coordinate, y(d,m,n)
    type(breadth)                         :: res
    type(mol_block_list)                  :: b_child
    integer(IK)                           :: mn, dmn, ispc, iper, isym, nper, nsym
!
    mn  = SUM(b%b%n * b%b%m)
    dmn = b%d * mn
    b_child = b%child()
    ispc = b_child%ispecies()
    nper = b_child%b(ispc)%g + 1
    nsym = rot(ispc)%nsym()
    allocate (res%nodes(nper * (nsym + 1)))
!
    do concurrent(iper=1:nper, isym=0:nsym)
      block
        real(RK) :: z(dmn), w(dmn)
        call copy(dmn, x, z)
        call centering(b%d, mn, z)
        call pack_child(rot(ispc), b_child, iper, isym, y, w)
        res%nodes(isym * nper + iper) = node(b_child, z, w)
      end block
    end do
!
  end function node_generate_breadth
!
  pure subroutine pack_child(r, b, iper, isym, x, w)
    integer(IK), intent(in)               :: iper, isym
    class(molecular_rotation), intent(in) :: r
    type(mol_block_list), intent(in)      :: b
    real(RK), intent(in)                  :: x(*)
    real(RK), intent(inout)               :: w(*)
    integer(IK)                           :: i, s, ispc, mn
!
    ispc = b%ispecies()
    s = b%nspecies()
    do concurrent(i=1:s)
      if (i == ispc) then
        call pack_x(r, b%d, b%b(i), iper, isym, x(b%b(i)%p), w(b%b(i)%p))
      else
        call copy(b%d * b%b(i)%m * b%b(i)%n, x(b%b(i)%p), w(b%b(i)%p))
      end if
    end do
!
    mn = SUM(b%b%m * b%b%n)
    call centering(b%d, mn, w)
!
  end subroutine pack_child
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
  pure subroutine pack_x(rot, d, b, iper, isym, y, w)
    class(molecular_rotation), intent(in) :: rot
    integer(IK), intent(in)               :: d, iper, isym
    type(mol_block), intent(in)           :: b
    real(RK), intent(in)                  :: y(d, b%m, *)
    real(RK), intent(inout)               :: w(d, b%m, *)
    integer(IK)                           :: i, j, k
    do concurrent(i=1:d, j=1:b%m, k=1:b%g)
      w(i, j, k) = y(i, j, k)
    end do
    do concurrent(i=1:d, j=1:b%m)
      w(i, j, b%g + 1) = y(i, j, b%g + iper)
    end do
    call rot%swap(d, w(1, 1, b%g + 1), isym)
    do concurrent(i=1:d, j=1:b%m, k=b%g + 1:b%g + iper - 1)
      w(i, j, k + 1) = y(i, j, k)
    end do
    do concurrent(i=1:d, j=1:b%m, k=b%g + iper + 1:b%n)
      w(i, j, k) = y(i, j, k)
    end do
  end subroutine pack_x
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
  pure elemental subroutine node_destroy(this)
    type(node), intent(inout) :: this
    this%lower = -HUGE(ZERO)
  end subroutine node_destroy
!
  pure elemental subroutine breadth_destroy(this)
    type(breadth), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
  end subroutine breadth_destroy
!
end module mod_tree

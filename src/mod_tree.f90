module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_lower_bound
  implicit none
  private
  public :: node, childs
!
  type node
!
    private
!
    integer(IK), allocatable :: fixed_indices(:)
    !! Molecules indices for fix
    integer(IK), allocatable :: free_indices(:)
    !! Molecules indices that can free rotation
    real(RK), public         :: lower = -HUGE(ZERO)
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
  pure function node_new(d, m, n, fixed, free, x, y) result(res)
    integer(IK), intent(in) :: d
    !! space dimension, d > 0
    integer(IK), intent(in) :: m
    !! number of atom per molecule, m > 0
    integer(IK), intent(in) :: n
    !! number of molecule, n > 0
    integer(IK), intent(in) :: fixed(:)
    !! Molecules indices for fix, starting from 1.
    integer(IK), intent(in) :: free(:)
    !! Molecules indices that can free rotation, starting from 1.
    real(RK), intent(in)    :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)    :: y(*)
    !! target molecular coordinate, y(d,m,n)
    real(RK)                :: w(node_worksize(d, m, n, free))
    type(node)              :: res
    integer(IK)             :: mn, dmn, ix, iy, iw
!
    res%fixed_indices = fixed
    res%free_indices = free
!
    mn  = m * n
    dmn = d * mn
!
    ix = 1
    iy = ix + dmn
    iw = iy + dmn
!
    w(ix:ix + dmn - 1) = x(:dmn)
    call copy_y(d, m, n, fixed, free, y, w(iy))
!
    call lower_bound(d, mn, nlist(m, free), w(ix), w(iy), w(iw))
!
    res%lower = w(iw)
!
  contains
!
    pure subroutine copy_y(d, m, n, fixed, free, y, z)
    integer(IK), intent(in) :: d, m, n, fixed(:), free(:)
    real(RK), intent(in)    :: y(d, m, n)
    real(RK), intent(inout) :: z(d, m, n)
    integer(IK)             :: i, j, k, f, g
!
      f = SIZE(fixed)
      g = SIZE(free)
!
      do concurrent(i=1:d, j=1:m, k=1:f)
        z(i, j, k) = y(i, j, fixed(k))
      enddo
!
      do concurrent(i=1:d, j=1:m, k=1:g)
        z(i, j, f + k) = y(i, j, free(k))
      enddo
!
    end subroutine copy_y
!
  end function node_new
!
!| generate childe nodes instance
  pure function node_generate_childs(this, d, m, n, x, y) result(res)
    class(node), intent(in) :: this
    integer(IK), intent(in) :: d
    !! space dimension, d > 0
    integer(IK), intent(in) :: m
    !! number of atom per molecule, m > 0
    integer(IK), intent(in) :: n
    !! number of molecule, n > 0
    real(RK), intent(in)    :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)    :: y(*)
    !! target molecular coordinate, y(d,m,n)
    type(childs)            :: res
    integer(IK)             :: f, g, i
!
    f = SIZE(this%free_indices)
    g = SIZE(this%fixed_indices)
!
    ALLOCATE(res%nodes(f))
!
    do concurrent(i=1:f)
      block
        integer(IK) :: fixed(g + 1), free(f - 1)
        fixed = [this%fixed_indices, this%free_indices(i)]
        free = [this%free_indices(:i - 1), this%free_indices(i + 1:)]
        res%nodes(i) = node(d, m, n, fixed, free, x, y)
      end block
    end do
!
  end function node_generate_childs
!
  pure function node_worksize(d, m, n, free) result(res)
  integer(IK), intent(in) :: d, m, n, free(:)
  integer(IK)             :: res
!
    res = 2 * (d * m * n) + lower_bound_worksize(d, m * n, nlist(m, free))
!
  end function node_worksize
!
  pure function nlist(m, free) result(res)
  integer(IK), intent(in) :: m, free(:)
  integer(IK)             :: res(SIZE(free) * m)
  integer(IK)             :: i, j
!
    do concurrent(j=1:SIZE(free), i=1:m)
      res((j - 1) * m + i) = (free(j) - 1) * m + i
    end do
!
  end function nlist
!
  pure elemental subroutine node_destroy(this)
    type(node), intent(inout) :: this
    if(ALLOCATED(this%fixed_indices)) deallocate(this%fixed_indices)
    if(ALLOCATED(this%free_indices))  deallocate(this%free_indices)
    this%lower = -HUGE(ZERO)
  end subroutine node_destroy
!
  pure elemental subroutine childs_destroy(this)
    type(childs), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
  end subroutine childs_destroy
!
end module mod_tree

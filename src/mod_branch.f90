module mod_branch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_lower_bound
  implicit none
!
  type node
!
    integer(IK), allocatable :: fixed_indices(:)
    !! Molecules indices for fix
    integer(IK), allocatable :: free_indices(:)
    !! Molecules indices that can free rotation
    real(RK)                 :: lower = -HUGE(ZERO)
    !! the lower bound
!
  contains
!
    final             :: node_destroy
!
  end type node
!
contains
!
!| generate instance
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
    integer(IK)             :: dnm, ix, iy, iw
!
    res%fixed_indices = fixed
    res%free_indices = free
!
    dmn = d * n * m
!
    ix = 1
    iy = ix + dmn
    iw = iy + dmn
!
    w(ix:ix + dmn - 1) = x(:dmn)
    call copy_y(d, m, n, fixed, free, y, w(iy))
!
    call lower_bound(d, n, nlist(m, free), w(ix), w(iy), w(iw))
!
    res%lower = w(iw)
!
  contains
!
    pure function node_worksize(d, m, n, free) result(res)
    integer(IK), intent(in) :: d, m, n, free(:)
    integer(IK)             :: i, j, res
!
      res = 2 * (d * m * n) + lower_bound_worksize(d, m * n, nlist(m, free))
!
    end function node_worksize
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
  end function node_new
!
  pure elemental subroutine node_destroy(this)
    type(node), intent(inout) :: this
    if(ALLOCATED(this%fixed_indices)) deallocate(this%fixed_indices)
    if(ALLOCATED(this%free_indices))  deallocate(this%free_indices)
    this%lower = -HUGE(ZERO)
  end subroutine node_destroy
!
end module mod_branch

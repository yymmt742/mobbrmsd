!| Define spatial dimension, \(D\),
!  and provide optimized functions for dimension.
module mod_dimspec_functions
  use mod_kinds, only: RK, IK
  implicit none
  private
  public  :: D, DD, ND
  public  :: setup_dimension
  public  :: compute_com
  public  :: covdot
  public  :: covcopy
  !| Spatial dimension, \(D\).
  integer(IK), protected, save :: D = 3
  !| Square spatial dimension, \(D^2\).
  integer(IK), protected, save :: DD = 9
  !| Node memory size, defined by \(1 + 1 + D^2\).
  !  Let \([L, G, \mathbf{C}]\) be a node,
  !  where \(L, G\in\mathbb{R}\) and \(\mathbf{C}\in\mathbb{R}^{D\times D}\).
  integer(IK), protected, save :: ND = 9 + 2
!
  real(RK), parameter :: ZERO = 0.0_RK
  real(RK), parameter :: ONE = 1.0_RK
!
contains
!| Sets the dimensions of the space. <br>
!  Caution, this routine affects global.
  subroutine setup_dimension(d_)
    integer(IK), intent(in) :: d_
    D = MAX(1, d_)
    DD = D * D
    ND = DD + 2
  end subroutine setup_dimension
!
!| Calculate center of mass.
  pure subroutine compute_com(d, n, X, C)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *)
    real(RK), intent(inout) :: C(d)
    real(RK)                :: rn
    integer(IK)             :: i, j
    do concurrent(i=1:d)
      C(i) = ZERO
    end do
    do j = 1, n
      do concurrent(i=1:d)
        C(i) = C(i) + X(i, j)
      end do
    end do
    rn = ONE / real(n, RK)
    do concurrent(i=1:d)
      C(i) = C(i) * rn
    end do
  end subroutine compute_com
!
!| Calculate \(\text{tr}[(\mathbf Y-\bar{\mathbf Y})(\mathbf X-\bar{\mathbf X})^\top]\).
  pure function covdot(d, n, X, Y, CX, CY) result(res)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), Y(d, *)
    real(RK), intent(in)    :: CX(d), CY(d)
    real(RK)                :: res, su(d)
    integer(IK)             :: i, j
    do concurrent(i=1:d)
      su(i) = ZERO
    end do
    do j = 1, n
      do concurrent(i=1:d)
        su(i) = su(i) + (X(i, j) - CX(i)) * (Y(i, j) - CY(i))
      end do
    end do
    res = ZERO
    do i = 1, d
      res = res + su(i)
    end do
  end function covdot
!
!| Compute \(\mathbf{Y} \gets \mathbf{X}-\bar{\mathbf{X}}\).
  pure subroutine covcopy(d, n, X, CX, Y)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), CX(d)
    real(RK), intent(inout) :: Y(d, *)
    integer(IK)             :: i, j
    do concurrent(i=1:d, j=1:n)
      Y(i, j) = X(i, j) - CX(i)
    end do
  end subroutine covcopy
end module mod_dimspec_functions


!| Calculate the covariance matrix. <br>
module mod_cov
  use mod_kinds, only: IK, RK
  implicit none
  private
  public :: covdot
  public :: covcopy
!
  real(RK), parameter    :: ZERO = 0.0_RK
!
contains
!| Calculate \(\text{tr}[(\mathbf Y-\bar{\mathbf Y})(\mathbf X-\bar{\mathbf X})^\top]\)
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
!| Compute \(\mathbf{Y} \gets \mathbf{X}-\bar{\mathbf{X}} \)
  pure subroutine covcopy(d, n, X, CX, Y)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), CX(d)
    real(RK), intent(inout) :: Y(d, *)
    integer(IK)             :: i, j
    do concurrent(i=1:d, j=1:n)
      Y(i, j) = X(i, j) - CX(i)
    end do
  end subroutine covcopy
!
end module mod_cov


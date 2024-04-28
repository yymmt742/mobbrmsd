!| Calculate the covariance matrix for \(D=3\). <br>
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
!  for \(D=3\)
  pure function covdot(d, n, X, Y, CX, CY) result(res)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), Y(d, *)
    real(RK), intent(in)    :: CX(d), CY(d)
    real(RK)                :: res, su1, su2
    integer(IK)             :: i
    res = ZERO
    do i = 1, n
      res = res + (X(1, i) - CX(1)) * (Y(1, i) - CY(1))
      su1 = su1 + (X(2, i) - CX(2)) * (Y(2, i) - CY(2))
      su2 = su2 + (X(3, i) - CX(3)) * (Y(3, i) - CY(3))
    end do
    res = res + su1 + su2
  end function covdot
!
!| Compute \(\mathbf{Y} \gets \mathbf{X}-\bar{\mathbf{X}} \)
!  for \(D=3\)
  pure subroutine covcopy(d, n, X, CX, Y)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), CX(d)
    real(RK), intent(inout) :: Y(d, *)
    integer(IK)             :: i
    do concurrent(i=1:n)
      Y(1, i) = X(1, i) - CX(1)
      Y(2, i) = X(2, i) - CX(2)
      Y(3, i) = X(3, i) - CX(3)
    end do
  end subroutine covcopy
!
end module mod_cov


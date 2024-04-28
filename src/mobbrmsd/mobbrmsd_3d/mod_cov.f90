!| Calculate the covariance matrix for \(D=3\). <br>
module mod_cov
  use mod_kinds, only: IK, RK
  implicit none
  private
  public :: com
  public :: covdot
  public :: covcopy
!
  real(RK), parameter :: ZERO = 0.0_RK
  real(RK), parameter :: ONE = 1.0_RK
!
contains
!| Calculate center of mass for \(D=3\)
  pure subroutine com(d, n, X, C)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *)
    real(RK), intent(inout) :: C(d)
    real(RK)                :: rn
    integer(IK)             :: i
    C(1) = ZERO
    C(2) = ZERO
    C(3) = ZERO
    do i = 1, n, 2
      C(1) = C(1) + X(1, i + 0) + X(1, i + 1)
      C(2) = C(2) + X(2, i + 0) + X(2, i + 1)
      C(3) = C(3) + X(3, i + 0) + X(3, i + 1)
    end do
    if (MODULO(n, 2) == 1) then
      C(1) = C(1) + X(1, n)
      C(2) = C(2) + X(2, n)
      C(3) = C(3) + X(3, n)
    end if
    rn = ONE / real(n, RK)
    C(1) = C(1) * rn
    C(2) = C(2) * rn
    C(3) = C(3) * rn
  end subroutine com
!
!| Calculate \(\text{tr}[(\mathbf Y-\bar{\mathbf Y})(\mathbf X-\bar{\mathbf X})^\top]\)
!  for \(D=3\)
  pure function covdot(d, n, X, Y, CX, CY) result(res)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), Y(d, *)
    real(RK), intent(in)    :: CX(d), CY(d)
    real(RK)                :: res, su1, su2
    integer(IK)             :: i
    res = ZERO
    su1 = ZERO
    su2 = ZERO
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


!| Calculate the rmsd.
module mod_rmsd
  use mod_params, only: IK, RK, ZERO => RZERO
  implicit none
  private
  public :: sd, msd, rmsd
!
contains
!
!| Calculate the root mean squared displacement.
  pure function rmsd(d, n, x, y) result(res)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: n
    !! matrix row dimension.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK)                :: res
!
    res = SQRT(msd(d, n, x, y))
!
  end function rmsd
!
!| Calculate the mean squared displacement.
  pure function msd(d, n, x, y) result(res)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: n
    !! matrix row dimension.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK)                :: res
!
    res = sd(d, n, x, y) / n
!
  end function msd
!
!| Calculate the sum of squared displacement.
  pure function sd(d, n, x, y) result(res)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: n
    !! matrix row dimension.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK)                :: res
    integer(IK)             :: i, nd
!
    nd = d * n
    res = ZERO
    do i = 1, nd
      res = res + (x(i) - y(i))**2
    enddo
!
  end function sd
!
end module mod_rmsd

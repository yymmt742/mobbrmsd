module mod_det
  use mod_params
  implicit none
  private
  public :: det
!
contains
!
   subroutine det(d, x, res)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: x(d, d)
    real(RK), intent(inout) :: res
!
    if (d < 1) then
!
      res = RZERO
!
    elseif (d == 1) then
!
      res = x(1, 1)
!
    elseif (d == 2) then
!
      res = x(1, 1) * x(2, 2) - x(1, 2) * x(2, 1)
!
    elseif (d == 3) then
!
      res = x(1, 1) * x(2, 2) * x(3, 3) &
        & + x(1, 2) * x(2, 3) * x(3, 1) &
        & + x(1, 3) * x(2, 1) * x(3, 2) &
        & - x(1, 3) * x(2, 2) * x(3, 1) &
        & - x(1, 2) * x(2, 1) * x(3, 3) &
        & - x(1, 1) * x(2, 3) * x(3, 2)
!
    else
!
      block
        integer(IK) :: i, ipiv(0:d)
        call DGETRF(d, d, x, d, ipiv(1:d), ipiv(0))
        res = RONE
        do i = 1, d
          res = res * x(i, i)
          if (ipiv(i) /= i) res = -res
        end do
      end block

    end if
!
  end subroutine det
!
end module mod_det

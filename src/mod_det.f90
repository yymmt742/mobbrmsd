module mod_det
  use mod_params, only : IK, RK, ONE => RONE, ZERO => RZERO
  implicit none
  private
  public :: det, det_sign
!
  interface
    include 'dgetrf.h'
  end interface
!
contains
!
   pure subroutine det(d, x, res)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: x(d, d)
    real(RK), intent(inout) :: res
!
    if (d < 1) then
!
      res = ZERO
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
      res = x(1, 1) * ( x(2, 2) * x(3, 3) - x(2, 3) * x(3, 2) ) +&
          & x(1, 2) * ( x(2, 3) * x(3, 1) - x(2, 1) * x(3, 3) ) +&
          & x(1, 3) * ( x(2, 1) * x(3, 2) - x(2, 2) * x(3, 1) )
!
    else
!
      block
        integer(IK) :: i, ipiv(0:d)
        call DGETRF(d, d, x, d, ipiv(1:d), ipiv(0))
        res = ONE
        do i = 1, d
          res = res * x(i, i)
          if (ipiv(i) /= i) res = -res
        end do
      end block

    end if
!
  end subroutine det
!
  pure subroutine det_sign(d, x, res)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: x(d, d)
    real(RK), intent(inout) :: res
!
    if (d < 1) then
!
      res = ONE
!
    elseif (d == 1) then
!
      res = SIGN(ONE, x(1, 1))
!
    elseif (d == 2) then
!
      res = SIGN(ONE, x(1, 1) * x(2, 2) - x(1, 2) * x(2, 1))
!
    elseif (d == 3) then
!
      res = SIGN(ONE, x(1, 1) * (x(2, 2) * x(3, 3) - x(2, 3) * x(3, 2))&
        &           + x(1, 2) * (x(2, 3) * x(3, 1) - x(2, 1) * x(3, 3))&
        &           + x(1, 3) * (x(2, 1) * x(3, 2) - x(2, 2) * x(3, 1)))
!
    else
!
      block
        integer(IK) :: i, ipiv(0:d)
        call DGETRF(d, d, x, d, ipiv(1:d), ipiv(0))
        res = MERGE( ONE, -ONE, &
                     MODULO(COUNT([(ipiv(i) == i, i=1, d)]) + &
            &        COUNT([(x(i, i) > ZERO, i=1, d)]),&
            &        2) == 0 )
      end block

    end if
!
  end subroutine det_sign
!
end module mod_det

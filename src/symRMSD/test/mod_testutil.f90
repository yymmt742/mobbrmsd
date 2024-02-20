!
!| Utility functions for testing.
module mod_testutil
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: sample
!
contains
!
  function sample(d, n) result(res)
    integer(IK), intent(in) :: d, n
    real(RK)                :: cnt(d)
    real(RK)                :: res(d, n)
    integer(IK)             :: i
    call RANDOM_NUMBER(res)
    cnt = SUM(res, 2) / n
    do concurrent(i=1:n)
      res(:, i) = res(:, i) - cnt
    end do
  end function sample
!
end module mod_testutil

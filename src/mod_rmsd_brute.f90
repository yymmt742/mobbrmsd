!| Calculate the rmsd.
module mod_rmsd_brute
  use mod_params, only: IK, RK, ZERO => RZERO
  use mod_permutation
  implicit none
  private
  public :: rmsd_brute
!
contains
!
!| Calculate the root mean squared displacement.
   function rmsd_brute(d, n, x, y) result(res)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: n
    !! matrix row dimension.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK)                :: res, tmp
    type(permutation)       :: perm
!
    perm = permutation(n)
    res = HUGE(ZERO)
!
    do while(.not.perm%endl())
      call sd_perm(d, n, perm%id, x, y, tmp)
      if (tmp < res) res = tmp
      call perm%next()
    enddo
    res = SQRT(res / n)
!
  end function rmsd_brute
!
  pure subroutine sd_perm(d, n, id, x, y, res)
    integer(IK), intent(in) :: d, n, id(n)
    real(RK), intent(in)    :: x(d, n), y(d, n)
    real(RK), intent(inout) :: res
    integer(IK)             :: i, j
    res = ZERO
    do j = 1, n
      do i = 1, d
        res = res + (x(i, j) - y(i, id(j)))**2
      end do
    end do
  end subroutine sd_perm
!
end module mod_rmsd_brute

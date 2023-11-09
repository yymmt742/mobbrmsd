!| Calculate the covariance matrix.
module mod_cov
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  implicit none
  private
  public :: cov, cov_row_major
!
  interface
    include 'dgemm.h'
  end interface
!
contains
!
  pure elemental function optarg(l) result(res)
    logical, intent(in), optional :: l
    logical                       :: res
    if (PRESENT(l)) then
      res = l
    else
      res = .false.
    end if
  end function optarg
!
!| Calculate the covariance matrix, X@YT.
  pure subroutine cov(d, n, x, y, res)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: n
    !! matrix row dimension.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK), intent(inout) :: res(*)
    !! returns d*d covariance matrix
!
    call DGEMM('N', 'T', d, d, n, ONE, x, d, y, d, ZERO, res, d)
!
  end subroutine cov
!
!| Calculate the covariance matrix, XT@Y.
  pure subroutine cov_row_major(d, n, x, y, res)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: n
    !! matrix row dimension.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK), intent(inout) :: res(*)
    !! returns n*n covariance matrix
!
    call DGEMM('T', 'N', n, n, d, ONE, x, d, y, d, ZERO, res, n)
!
  end subroutine cov_row_major
!
end module mod_cov

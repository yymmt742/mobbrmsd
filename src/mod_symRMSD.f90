module mod_symRMSD
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_lower_bound
  implicit none
!
contains
!
! subroutine X()
! type(unittest)               :: u
! logical,parameter            :: T=.TRUE., F=.FALSE.
! call u%init( 'test_unittest' )
! call u%assert(                T,           'assert               bool    0  ' )
! call u%finish_and_terminate()
! end subroutine
!
!| Calculate the rotation matrix from covariance matrix.
  subroutine branch_and_plune(d, n, free_indices, X, Y)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)       :: n
    !! matrix row dimension.
    integer(IK), intent(in)       :: free_indices(:)
    !! matrix row dimension.
    real(RK), intent(in)          :: X(*)
    !! target d*n array
    real(RK), intent(inout)       :: Y(*)
    !! rotation d*d matrix
!
!
  end subroutine branch_and_plune
!
end module mod_symRMSD

module mod_symRMSD
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_tree
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
  subroutine branch_and_plune(d, m, n, X, Y)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)       :: m
    !! matrix collumn dimension.
    integer(IK), intent(in)       :: n
    !! matrix row dimension.
    real(RK), intent(in)          :: X(*)
    !! reference d*m*n array
    real(RK), intent(inout)       :: Y(*)
    !! target d*m*n array
!
!
  end subroutine branch_and_plune
!
end module mod_symRMSD

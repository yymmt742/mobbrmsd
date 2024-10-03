pure subroutine DLASQ1( N, D, E, WORK, INFO )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: N
integer,intent(out)    :: INFO
real(DP),intent(inout) :: D( * ), E( * )
real(DP),intent(out)   :: WORK( * )
end subroutine DLASQ1

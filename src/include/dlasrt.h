pure subroutine DLASRT( ID, N, D, INFO )
use LA_CONSTANTS, only: DP
character(*),intent(in) :: ID
integer,intent(in)      :: N
integer,intent(out)     :: INFO
real(DP),intent(inout)  :: D( * )
end subroutine DLASRT

pure subroutine DLANGE( NORM, M, N, A, LDA, RES, WORK )
use LA_CONSTANTS, only: DP
character,intent(in)   :: NORM
integer,intent(in)     :: LDA, M, N
real(DP),intent(inout) :: A( LDA, * )
real(DP),intent(out)   :: RES, WORK( * )
end subroutine DLANGE

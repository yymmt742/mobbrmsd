pure subroutine DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: K, LDA, M, N
integer,intent(out)    :: INFO
real(DP),intent(in)    :: TAU( * )
real(DP),intent(inout) :: A( LDA, * )
real(DP),intent(out)   :: WORK( * )
end subroutine DORG2R


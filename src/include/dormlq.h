pure subroutine DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, &
               &        LDC, WORK, LWORK, INFO )
use LA_CONSTANTS, only: DP
character,intent(in)   :: SIDE, TRANS
integer,intent(in)     :: K, LDA, LDC, LWORK, M, N
integer,intent(out)    :: INFO
real(DP),intent(in)    :: TAU( * )
real(DP),intent(inout) :: A( LDA, * ), C( LDC, * )
real(DP),intent(out)   :: WORK( * )
end subroutine DORMLQ


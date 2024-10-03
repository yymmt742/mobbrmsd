pure subroutine DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
use LA_CONSTANTS, only: DP
character,intent(in)   :: VECT
INTEGER,intent(in)     :: K, LDA, LWORK, M, N
INTEGER,intent(out)    :: INFO
real(DP),intent(in)    :: TAU( * )
real(DP),intent(inout) :: A( LDA, * )
real(DP),intent(out)   :: WORK( * )
end subroutine DORGBR

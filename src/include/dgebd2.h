pure subroutine DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: LDA, M, N
integer,intent(out)    :: INFO
real(DP),intent(inout) :: A( LDA, * )
real(DP),intent(out)   :: D( * ), E( * ), TAUP( * ), TAUQ( * ), WORK( * )
end subroutine DGEBD2

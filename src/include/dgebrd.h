pure subroutine DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: LDA, LWORK, M, N
integer,intent(out)    :: INFO
real(DP),intent(inout) :: A( LDA, * )
real(DP),intent(out)   :: D( * ), E( * ), TAUP( * ),TAUQ( * ), WORK( * )
end subroutine DGEBRD

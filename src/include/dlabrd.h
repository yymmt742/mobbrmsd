pure subroutine DLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: LDA, LDX, LDY, M, N, NB
real(DP),intent(inout) :: A( LDA, * )
real(DP),intent(out)   :: D( * ), E( * ), TAUP( * ), TAUQ( * ), X( LDX, * ), Y( LDY, * )
end subroutine DLABRD

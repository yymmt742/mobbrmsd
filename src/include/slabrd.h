pure subroutine SLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY )
integer,intent(in) :: LDA, LDX, LDY, M, N, NB
real,intent(inout) :: A( LDA, * )
real,intent(out)   :: D( * ), E( * ), TAUP( * ), TAUQ( * ), X( LDX, * ), Y( LDY, * )
end subroutine SLABRD

pure subroutine SGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
integer,intent(in)  :: LDA, M, N
integer,intent(out) :: INFO
real,intent(inout)  :: A( LDA, * )
real,intent(out)    :: D( * ), E( * ), TAUP( * ), TAUQ( * ), WORK( * )
end subroutine SGEBD2

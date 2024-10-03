pure subroutine SGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO )
integer,intent(in)  :: LDA, LWORK, M, N
integer,intent(out) :: INFO
real,intent(inout)  :: A( LDA, * )
real,intent(out)    :: D( * ), E( * ), TAUP( * ),TAUQ( * ), WORK( * )
end subroutine SGEBRD

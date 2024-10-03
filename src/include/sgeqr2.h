pure subroutine SGEQR2( M, N, A, LDA, TAU, WORK, INFO )
integer,intent(in)  :: LDA, M, N
integer,intent(out) :: INFO
real,intent(inout)  :: A( LDA, * )
real,intent(out)    :: TAU( * ), WORK( * )
end subroutine SGEQR2

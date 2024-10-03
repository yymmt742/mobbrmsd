pure subroutine SGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
integer,intent(in)  :: LDA, LWORK, M, N
integer,intent(out) :: INFO
real,intent(inout)  :: A( LDA, * )
real,intent(out)    :: TAU( * ), WORK( * )
end subroutine SGELQF

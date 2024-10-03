pure subroutine SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
integer,intent(in)  :: K, LDA, LWORK, M, N
integer,intent(out) :: INFO
real,intent(in)     :: TAU( * )
real,intent(inout)  :: A( LDA, * )
real,intent(out)    :: WORK( * )
end subroutine SORGQR

pure subroutine SORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
character,intent(in) :: VECT
INTEGER,intent(in)   :: K, LDA, LWORK, M, N
INTEGER,intent(out)  :: INFO
real,intent(in)      :: TAU( * )
real,intent(inout)   :: A( LDA, * )
real,intent(out)     :: WORK( * )
end subroutine SORGBR

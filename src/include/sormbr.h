pure subroutine SORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, &
               &        C, LDC, WORK, LWORK, INFO )
character,intent(in) :: SIDE, TRANS, VECT
integer,intent(in)   :: K, LDA, LDC, LWORK, M, N
integer,intent(out)  :: INFO
real,intent(in)      :: TAU( * )
real,intent(inout)   :: A( LDA, * ), C( LDC, * )
real,intent(out)     :: WORK( * )
end subroutine SORMBR

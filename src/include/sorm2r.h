pure subroutine SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, &
               &        LDC, WORK, INFO )
character,intent(in) :: SIDE, TRANS
integer,intent(in)   :: K, LDA, LDC, M, N
integer,intent(out)  :: INFO
real,intent(in)      :: TAU( * )
real,intent(inout)   :: A( LDA, * ), C( LDC, * )
real,intent(out)     :: WORK( * )
end subroutine SORM2R


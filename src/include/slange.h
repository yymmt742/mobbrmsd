pure subroutine SLANGE( NORM, M, N, A, LDA, RES, WORK )
character,intent(in) :: NORM
integer,intent(in)   :: LDA, M, N
real,intent(inout)   :: A( LDA, * )
real,intent(out)     :: RES, WORK( * )
end subroutine SLANGE

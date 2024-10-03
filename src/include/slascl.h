pure subroutine SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
character,intent(in) :: TYPE
integer,intent(in)   :: KL, KU, LDA, M, N
integer,intent(out)  :: INFO
real,intent(in)      :: CFROM, CTO
real,intent(inout)   :: A( LDA, * )
end subroutine SLASCL

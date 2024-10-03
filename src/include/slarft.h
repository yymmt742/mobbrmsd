pure subroutine SLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
character,intent(in) :: DIRECT, STOREV
integer,intent(in)   :: K, LDT, LDV, N
real,intent(in)      :: TAU( * ), V( LDV, * )
real,intent(out)     :: T( LDT, * )
end subroutine SLARFT

pure subroutine DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
use LA_CONSTANTS, only: DP
character,intent(in) :: DIRECT, STOREV
integer,intent(in)   :: K, LDT, LDV, N
real(DP),intent(in)  :: TAU( * ), V( LDV, * )
real(DP),intent(out) :: T( LDT, * )
end subroutine DLARFT

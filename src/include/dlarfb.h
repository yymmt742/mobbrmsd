pure subroutine DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, &
               &        LDV, T, LDT, C, LDC, WORK, LDWORK )
use LA_CONSTANTS, only: DP
character,intent(in)   :: DIRECT, SIDE, STOREV, TRANS
integer,intent(in)     :: K, LDC, LDT, LDV, LDWORK, M, N
real(DP),intent(inout) :: C( LDC, * )
real(DP),intent(in)    :: T( LDT, * ), V( LDV, * )
real(DP),intent(out)   :: WORK( LDWORK, * )
end subroutine DLARFB

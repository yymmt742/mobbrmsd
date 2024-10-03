pure subroutine DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
use LA_CONSTANTS, only: DP
character,intent(in)   :: TYPE
integer,intent(in)     :: KL, KU, LDA, M, N
integer,intent(out)    :: INFO
real(DP),intent(in)    :: CFROM, CTO
real(DP),intent(inout) :: A( LDA, * )
end subroutine DLASCL

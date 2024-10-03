pure subroutine DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
use LA_CONSTANTS, only: DP
character,intent(in) :: UPLO
integer,intent(in)   :: LDA, M, N
real(DP),intent(in)  :: ALPHA, BETA
real(DP),intent(out) :: A( LDA, * )
end subroutine DLASET

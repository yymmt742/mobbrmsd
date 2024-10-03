pure subroutine DLACPY( UPLO, M, N, A, LDA, B, LDB )
use LA_CONSTANTS, only: DP
character(*),intent(in) :: UPLO
integer,intent(in)      :: LDA, LDB, M, N
real(DP),intent(in)     :: A( LDA, * )
real(DP),intent(out)    :: B( LDB, * )
end subroutine DLACPY

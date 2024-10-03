pure recursive subroutine DGETRF2( M, N, A, LDA, IPIV, INFO )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: LDA, M, N
integer,intent(out)    :: INFO
integer,intent(out)    :: IPIV( * )
real(DP),intent(inout) :: A( LDA, * )
end subroutine DGETRF2


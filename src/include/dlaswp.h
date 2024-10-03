pure subroutine DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
use LA_CONSTANTS, only: dp
integer,intent(in)     :: INCX, K1, K2, LDA, N
integer,intent(in)     :: IPIV( * )
real(DP),intent(inout) :: A( LDA, * )
end subroutine DLASWP

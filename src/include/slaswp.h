pure subroutine SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)     :: INCX, K1, K2, LDA, N
integer,intent(in)     :: IPIV( * )
real(SP),intent(inout) :: A( LDA, * )
end subroutine SLASWP

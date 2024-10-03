pure recursive subroutine SGETRF2( M, N, A, LDA, IPIV, INFO )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)     :: LDA, M, N
integer,intent(out)    :: INFO
integer,intent(out)    :: IPIV( * )
real(SP),intent(inout) :: A( LDA, * )
end subroutine SGETRF2


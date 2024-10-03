pure subroutine SORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)     :: K, LDA, M, N
integer,intent(out)    :: INFO
real(SP),intent(in)    :: TAU( * )
real(SP),intent(inout) :: A( LDA, * )
real(SP),intent(out)   :: WORK( * )
end subroutine SORG2R


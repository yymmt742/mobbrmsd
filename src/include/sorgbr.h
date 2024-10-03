pure subroutine SORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
character,intent(in)   :: VECT
INTEGER,intent(in)     :: K, LDA, LWORK, M, N
INTEGER,intent(out)    :: INFO
real(SP),intent(in)    :: TAU( * )
real(SP),intent(inout) :: A( LDA, * )
real(SP),intent(out)   :: WORK( * )
end subroutine SORGBR

pure subroutine SORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, &
               &        LDC, WORK, INFO )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
character(*),intent(in) :: SIDE, TRANS
integer,intent(in)      :: K, LDA, LDC, M, N
integer,intent(out)     :: INFO
real(SP),intent(in)     :: TAU( * )
real(SP),intent(inout)  :: A( LDA, * ), C( LDC, * )
real(SP),intent(out)    :: WORK( * )
end subroutine SORML2


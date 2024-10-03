pure subroutine SBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, &
               &        LDU, C, LDC, WORK, INFO )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
CHARACTER,intent(in)   :: UPLO
INTEGER,intent(in)     :: LDC, LDU, LDVT, N, NCC, NCVT, NRU
INTEGER,intent(out)    :: INFO
real(SP),intent(inout) :: C( LDC, * ), D( * ), E( * ), U( LDU, * ), VT( LDVT, * )
real(SP),intent(out)   :: WORK( * )
end subroutine SBDSQR

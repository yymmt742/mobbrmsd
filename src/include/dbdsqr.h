pure subroutine DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, &
               &        LDU, C, LDC, WORK, INFO )
use LA_CONSTANTS, only: DP
CHARACTER,intent(in)    :: UPLO
INTEGER,intent(in)      :: LDC, LDU, LDVT, N, NCC, NCVT, NRU
INTEGER,intent(out)     :: INFO
real(DP),intent(inout)  :: C( LDC, * ), D( * ), E( * ), U( LDU, * ), VT( LDVT, * )
real(DP),intent(out)    :: WORK( * )
end subroutine DBDSQR

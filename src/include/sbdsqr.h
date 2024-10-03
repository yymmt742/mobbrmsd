pure subroutine SBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, &
               &        LDU, C, LDC, WORK, INFO )
CHARACTER,intent(in) :: UPLO
INTEGER,intent(in)   :: LDC, LDU, LDVT, N, NCC, NCVT, NRU
INTEGER,intent(out)  :: INFO
real,intent(inout)   :: C( LDC, * ), D( * ), E( * ), U( LDU, * ), VT( LDVT, * )
real,intent(out)     :: WORK( * )
end subroutine SBDSQR

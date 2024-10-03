pure subroutine SLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, &
               &        LDV, T, LDT, C, LDC, WORK, LDWORK )
character,intent(in) :: DIRECT, SIDE, STOREV, TRANS
integer,intent(in)   :: K, LDC, LDT, LDV, LDWORK, M, N
real,intent(inout)   :: C( LDC, * )
real,intent(in)      :: T( LDT, * ), V( LDV, * )
real,intent(out)     :: WORK( LDWORK, * )
end subroutine SLARFB

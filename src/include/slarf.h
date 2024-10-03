pure subroutine SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
character(*),intent(in) :: SIDE
integer,intent(in)      :: INCV, LDC, M, N
real,intent(in)         :: TAU
real,intent(in)         :: V( * )
real,intent(inout)      :: C( LDC, * )
real,intent(out)        :: WORK( LDC, * )
end subroutine SLARF

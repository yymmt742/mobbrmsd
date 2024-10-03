pure subroutine DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
use LA_CONSTANTS, only: DP
character(*),intent(in) :: SIDE
integer,intent(in)      :: INCV, LDC, M, N
real(DP),intent(in)     :: TAU
real(DP),intent(in)     :: V( * )
real(DP),intent(inout)  :: C( LDC, * )
real(DP),intent(out)    :: WORK( LDC, * )
end subroutine DLARF

pure subroutine DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
use LA_CONSTANTS, only: DP
character(*),intent(in) :: DIRECT, PIVOT, SIDE
integer,intent(in)      :: LDA, M, N
real(DP),intent(in)     :: C( * ), S( * )
real(DP),intent(inout)  :: A( LDA, * )
end subroutine DLASR

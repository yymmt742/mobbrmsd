pure subroutine DGELQ2( M, N, A, LDA, TAU, WORK, INFO )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: LDA, M, N
integer,intent(out)    :: INFO
real(DP),intent(inout) :: A( LDA, * )
real(DP),intent(out)   :: TAU( * ), WORK( * )
end subroutine DGELQ2

pure subroutine DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: LDA, LWORK, M, N
integer,intent(out)    :: INFO
real(DP),intent(inout) :: A( LDA, * )
real(DP),intent(out)   :: TAU( * ), WORK( * )
end subroutine DGELQF

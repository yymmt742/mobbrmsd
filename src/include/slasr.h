pure subroutine SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
character(*),intent(in) :: DIRECT, PIVOT, SIDE
integer,intent(in)      :: LDA, M, N
real,intent(in)         :: C( * ), S( * )
real,intent(inout)      :: A( LDA, * )
end subroutine SLASR

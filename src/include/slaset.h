pure subroutine SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
character, intent(in) :: UPLO
integer, intent(in)   :: LDA, M, N
real, intent(in)      :: ALPHA, BETA
real, intent(out)     :: A( LDA, * )
end subroutine SLASET

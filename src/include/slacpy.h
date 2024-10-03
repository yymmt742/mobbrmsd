pure subroutine SLACPY( UPLO, M, N, A, LDA, B, LDB )
character(*),intent(in) :: UPLO
integer,intent(in)      :: LDA, LDB, M, N
real,intent(in)         :: A( LDA, * )
real,intent(out)        :: B( LDB, * )
end subroutine SLACPY

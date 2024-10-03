pure subroutine STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
character,intent(in) :: DIAG,SIDE,TRANSA,UPLO
integer,intent(in)   :: LDA,LDB,M,N
real,intent(in)      :: ALPHA
real,intent(in)      :: A(LDA,*)
real,intent(inout)   :: B(LDB,*)
end subroutine STRMM

pure subroutine STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
character,intent(in) :: DIAG,SIDE,TRANSA,UPLO
real,intent(in)      :: ALPHA
integer,intent(in)   :: LDA,LDB,M,N
real,intent(in)      :: A(LDA,*)
real,intent(inout)   :: B(LDB,*)
end subroutine STRSM

pure subroutine STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
integer,intent(in)   :: INCX,LDA,N
character,intent(in) :: DIAG,TRANS,UPLO
real,intent(in)      :: A(LDA,*)
real,intent(inout)   :: X(*)
end subroutine STRMV

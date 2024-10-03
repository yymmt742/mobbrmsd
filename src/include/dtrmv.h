pure subroutine DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
use LA_CONSTANTS, only: DP
integer,intent(in)     :: INCX,LDA,N
character,intent(in)   :: DIAG,TRANS,UPLO
real(DP),intent(in)    :: A(LDA,*)
real(DP),intent(inout) :: X(*)
end subroutine DTRMV

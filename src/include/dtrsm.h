pure subroutine DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
use LA_CONSTANTS, only: DP
character,intent(in)   :: DIAG,SIDE,TRANSA,UPLO
real(DP),intent(in)    :: ALPHA
integer,intent(in)     :: LDA,LDB,M,N
real(DP),intent(in)    :: A(LDA,*)
real(DP),intent(inout) :: B(LDB,*)
end subroutine DTRSM

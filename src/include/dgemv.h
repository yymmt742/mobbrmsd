pure subroutine DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
use LA_CONSTANTS, only: DP
real(DP),intent(in)    :: ALPHA,BETA
integer,intent(in)     :: INCX,INCY,LDA,M,N
CHARACTER,intent(in)   :: TRANS
real(DP),intent(in)    :: A(LDA,*),X(*)
real(DP),intent(inout) :: Y(*)
end subroutine DGEMV

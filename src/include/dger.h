pure subroutine DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
use LA_CONSTANTS, only: DP
real(DP),intent(in)    :: ALPHA
integer,intent(in)     :: INCX,INCY,LDA,M,N
real(DP),intent(in)    :: X(*),Y(*)
real(DP),intent(inout) :: A(LDA,*)
end subroutine DGER

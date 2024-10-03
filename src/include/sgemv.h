pure subroutine SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
real,intent(in)      :: ALPHA, BETA
integer,intent(in)   :: INCX, INCY, LDA, M, N
CHARACTER,intent(in) :: TRANS
real,intent(in)      :: A(LDA,*), X(*)
real,intent(inout)   :: Y(*)
end subroutine SGEMV

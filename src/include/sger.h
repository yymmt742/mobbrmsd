pure subroutine SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
real,intent(in)    :: ALPHA
integer,intent(in) :: INCX,INCY,LDA,M,N
real,intent(in)    :: X(*),Y(*)
real,intent(inout) :: A(LDA,*)
end subroutine SGER

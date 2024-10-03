pure subroutine SROT(N,DX,INCX,DY,INCY,C,S)
real,intent(in)    :: C,S
integer,intent(in) :: INCX,INCY,N
real,intent(inout) :: DX(*),DY(*)
end subroutine SROT

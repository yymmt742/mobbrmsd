pure subroutine DROT(N,DX,INCX,DY,INCY,C,S)
use LA_CONSTANTS, only: DP
real(DP),intent(in)    :: C,S
integer,intent(in)     :: INCX,INCY,N
real(DP),intent(inout) :: DX(*),DY(*)
end subroutine DROT

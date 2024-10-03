pure subroutine DSWAP(N,DX,INCX,DY,INCY)
use LA_CONSTANTS, only: DP
integer,intent(in)     :: INCX,INCY,N
real(DP),intent(inout) :: DX(*),DY(*)
end subroutine DSWAP

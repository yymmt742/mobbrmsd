pure subroutine DCOPY(N,DX,INCX,DY,INCY)
use LA_CONSTANTS, only: DP
integer,intent(in)   :: INCX,INCY,N
real(DP),intent(in)  :: DX(*)
real(DP),intent(out) :: DY(*)
end subroutine DCOPY

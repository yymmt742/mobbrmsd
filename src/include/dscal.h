pure subroutine DSCAL(N,DA,DX,INCX)
use LA_CONSTANTS, only: DP
integer,intent(in)     :: INCX,N
real(DP),intent(in)    :: DA
real(DP),intent(inout) :: DX(*)
end subroutine DSCAL

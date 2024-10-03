pure function IDAMAX(N,DX,INCX)
use LA_CONSTANTS, only: DP
integer,intent(in)  :: INCX,N
real(DP),intent(in) :: DX(*)
integer             :: IDAMAX
end function IDAMAX


pure function ISAMAX(N, DX, INCX)
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)  :: INCX, N
real(SP),intent(in) :: DX(*)
integer             :: ISAMAX
end function ISAMAX


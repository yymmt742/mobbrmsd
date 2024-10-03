pure subroutine SLASSQ( N, X, INCX, SCL, SUMSQ )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)     :: INCX, N
real(SP),intent(inout) :: SCL, SUMSQ
real(SP),intent(in)    :: X(*)
end subroutine SLASSQ

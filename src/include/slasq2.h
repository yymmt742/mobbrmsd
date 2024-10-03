pure subroutine SLASQ2( N, Z, INFO )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)  ::  N
integer,intent(out) ::  INFO
real,intent(inout)  :: Z( * )
end subroutine SLASQ2

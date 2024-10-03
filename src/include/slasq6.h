pure subroutine SLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, &
               &        DNM1, DNM2 )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)     :: I0, N0, PP
real(SP),intent(out)   :: DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
real(SP),intent(inout) :: Z( * )
end subroutine SLASQ6


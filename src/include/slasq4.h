pure subroutine SLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, &
               &        DN1, DN2, TAU, TTYPE, G )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)     :: I0, N0, PP, N0IN
integer,intent(out)    :: TTYPE
real(SP),intent(in)    :: DMIN, DMIN1, DMIN2, DN, DN1, DN2
real(SP),intent(inout) :: G
real(SP),intent(out)   :: TAU
real(SP),intent(in)    :: Z( * )
end subroutine SLASQ4


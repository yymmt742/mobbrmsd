pure subroutine SLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, &
               &        DN, DNM1, DNM2, IEEE, EPS )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
logical,intent(in)     :: IEEE
integer,intent(in)     :: I0, N0, PP
real(SP),intent(in)    :: SIGMA, EPS
real(SP),intent(out)   :: TAU, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
real(SP),intent(inout) :: Z( * )
end subroutine SLASQ5


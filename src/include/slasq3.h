pure subroutine SLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, &
               &        ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, &
               &        DN2, G, TAU )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)     :: I0
integer,intent(inout)  :: N0, PP, NFAIL, ITER, NDIV, TTYPE
real(SP),intent(out)   :: DMIN,  SIGMA
real(SP),intent(inout) :: DESIG, QMAX, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
logical,intent(in)     :: IEEE
real(SP),intent(inout) :: Z( * )
end subroutine SLASQ3


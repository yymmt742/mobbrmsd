pure subroutine DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, &
               &        ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, &
               &        DN2, G, TAU )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: I0
integer,intent(inout)  :: N0, PP, NFAIL, ITER, NDIV,TTYPE
real(DP),intent(out)   :: DMIN,  SIGMA
real(DP),intent(inout) :: DESIG, QMAX, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
logical,intent(in)     :: IEEE
real(DP),intent(inout) :: Z( * )
end subroutine DLASQ3


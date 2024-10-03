pure subroutine DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, &
               &        DN, DNM1, DNM2, IEEE, EPS )
use LA_CONSTANTS, only: DP
logical,intent(in)     :: IEEE
integer,intent(in)     :: I0, N0, PP
real(DP),intent(in)    :: SIGMA, EPS
real(DP),intent(out)   :: TAU, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
real(DP),intent(inout) :: Z( * )
end subroutine DLASQ5


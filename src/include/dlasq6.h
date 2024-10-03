pure subroutine DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, &
               &        DNM1, DNM2 )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: I0, N0, PP
real(DP),intent(out)   :: DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
real(DP),intent(inout) :: Z( * )
end subroutine DLASQ6


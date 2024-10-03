pure function DNRM2( N, X, INCX )
use LA_CONSTANTS, only: DP
integer,intent(in)  :: INCX, N
real(DP),intent(in) :: X(*)
real(DP)            :: DNRM2
end function DNRM2


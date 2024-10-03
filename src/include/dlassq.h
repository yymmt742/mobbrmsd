pure subroutine DLASSQ( n, x, incx, scl, sumsq )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: incx, n
real(DP),intent(inout) :: scl, sumsq
real(DP),intent(in)    :: x(*)
end subroutine DLASSQ

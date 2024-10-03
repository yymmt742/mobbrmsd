pure subroutine DLARFG( N, ALPHA, X, INCX, TAU )
use LA_CONSTANTS, only: DP
integer,intent(in)     :: INCX, N
real(DP),intent(inout) :: ALPHA
real(DP),intent(out)   :: TAU
real(DP),intent(inout) :: X( * )
end subroutine DLARFG

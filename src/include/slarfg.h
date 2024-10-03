pure subroutine SLARFG( N, ALPHA, X, INCX, TAU )
integer,intent(in) :: INCX, N
real,intent(inout) :: ALPHA
real,intent(out)   :: TAU
real,intent(inout) :: X( * )
end subroutine SLARFG

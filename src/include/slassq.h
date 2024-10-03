pure subroutine SLASSQ( N, X, INCX, SCL, SUMSQ )
integer,intent(in) :: INCX, N
real,intent(inout) :: SCL, SUMSQ
real,intent(in)    :: X(*)
end subroutine SLASSQ

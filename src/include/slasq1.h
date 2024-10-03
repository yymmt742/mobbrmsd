pure subroutine SLASQ1( N, D, E, WORK, INFO )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
integer,intent(in)  :: N
integer,intent(out) :: INFO
real,intent(inout)  :: D( * ), E( * )
real,intent(out)    :: WORK( * )
end subroutine SLASQ1

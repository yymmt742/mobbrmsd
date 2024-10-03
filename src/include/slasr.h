pure subroutine SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32
character(*),intent(in) :: DIRECT, PIVOT, SIDE
integer,intent(in)      :: LDA, M, N
real(SP),intent(in)     :: C( * ), S( * )
real(SP),intent(inout)  :: A( LDA, * )
end subroutine SLASR

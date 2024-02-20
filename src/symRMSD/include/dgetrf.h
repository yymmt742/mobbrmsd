pure subroutine DGETRF( M, N, A, LDA, IPIV, INFO )
  use, intrinsic :: ISO_FORTRAN_ENV, only: wp=>REAL64
  integer,intent(in)      :: LDA, M, N
  integer,intent(out)     :: INFO
  integer,intent(out)     :: IPIV( * )
  real(wp),intent(inout)  :: A( LDA, * )
end subroutine DGETRF
